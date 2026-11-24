!-------------------------------------------------------------------------------
!
!+  Program driver : NSW6 cloud-microphysics mini-app
!
!-------------------------------------------------------------------------------
program prg_mp_miniapp
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       Stand-alone driver for ONE call of the NICAM cloud-microphysics
  !       package, mod_mp_driver/mp_driver (negative filter -> NSW6 core ->
  !       precipitation transport -> saturation adjustment -> cloud-fraction
  !       diagnosis).
  !
  !       This is the mini-app counterpart of nhm/driver/prg_driver.f90: the
  !       same setup chain, the same profiler, but every field that NICAM would
  !       have produced by reading a restart file and by building the grid is
  !       instead read from the Phase-A dump under nicam_mp/dataset/<size>/,
  !       raw Fortran-unformatted big-endian records (see nicam_mp/README.md).
  !
  !       The whole NSW6 chain is column-local: no halo exchange, no horizontal
  !       operator, no pole (_pl) handling, so a single process is enough and
  !       the region index l is nothing but the slowest array dimension.
  !
  !       Flow:
  !          ADM_openparam -> ADM_proc_init -> ADM_setup -> PROF_setup
  !          -> CNST_setup -> GRD_setup -> runconf_setup -> SATURATION_setup
  !          -> mp_init(MP_TYPE)
  !          -> load the 16 static + 15 per-step input records
  !          -> repeat: restore the pristine state, call mp_driver per region
  !          -> validate the 23 mp_driver outputs against the dump
  !          -> (optionally) write the computed outputs back as big-endian records
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      1.00      2026-08    Phase B: extracted from NICAM (mini-app)
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio
  use mod_prof
  use mod_adm, only: &
     GLOBAL_prc_me,   &
     GLOBAL_prc_master, &
     ADM_openparam,   &
     ADM_proc_init,   &
     ADM_proc_stop,   &
     ADM_proc_finish, &
     ADM_setup,       &
     ADM_KNONE,       &
     ADM_gall_in,     &
     ADM_kall,        &
     ADM_lall
  use mod_cnst, only: &
     CNST_EPS_DP, & !--- MINIAPP: validator zero-threshold
     CNST_setup
  use mod_grd, only: &
     GRD_setup, & ! provides GRD_rdgz / rdgzh / afact..dfact used by the MP chain
     GRD_gz,    & ! [OpenACC]
     GRD_gzh,   & ! [OpenACC]
     GRD_dgz,   & ! [OpenACC]
     GRD_dgzh,  & ! [OpenACC]
     GRD_rdgz,  & ! [OpenACC]
     GRD_rdgzh, & ! [OpenACC]
     GRD_afact, & ! [OpenACC]
     GRD_bfact, & ! [OpenACC]
     GRD_cfact, & ! [OpenACC]
     GRD_dfact    ! [OpenACC]
  use mod_runconf, only: &
     runconf_setup, &
     MP_TYPE,       &
     TRC_VMAX,      &
     TRC_name,      & !--- MINIAPP: per-category validation row labels
     CPW,           & ! [OpenACC]
     CVW              ! [OpenACC]
  use mod_bndcnd, only: &
     BNDCND_setup
  use mod_satadjust, only: &
     SATURATION_setup
  use mod_vadv1d, only: &
     vadv1d_init
  use mod_mp_driver, only: &
     mp_init, &
     mp_driver
  !-----------------------------------------------------------------------------
  implicit none

  character(len=H_FILE) :: param_fname = "mp_miniapp.cnf" ! name of parameter file

  !--- &MP_MINIAPPPARAM : the mini-app's own knobs -----------------------------
  character(len=H_LONG) :: dataset_dir  = './dataset/region0' !< where the mpF_*.bin records live
  character(len=H_LONG) :: output_dir   = './output'          !< where to write our own records
  integer               :: num_repeat   = 1                   !< how many times to run mp_driver (timing)
  real(DP)              :: tolerance    = 1.E-13_DP           !< relative tolerance of the validation
  real(DP)              :: tolerance_abs= 0.0_DP              !< absolute tolerance of the validation
  logical               :: validate     = .true.              !< compare against the dumped *_out records
  logical               :: write_output = .false.             !< write our own big-endian records

  namelist / MP_MINIAPPPARAM / &
     dataset_dir,   &
     output_dir,    &
     num_repeat,    &
     tolerance,     &
     tolerance_abs, &
     validate,      &
     write_output

  !--- dimensions (from ADM_setup / runconf_setup) -----------------------------
  integer :: ijdim, kdim, lall, nqmax

  !--- static input (16 records) -----------------------------------------------
  real(DP), allocatable :: z      (:,:,:), zh     (:,:,:), phi    (:,:,:)
  real(DP), allocatable :: gam2   (:,:,:), gam2h  (:,:,:)
  real(DP), allocatable :: gsgam2 (:,:,:), gsgam2h(:,:,:)
  real(DP), allocatable :: ix(:,:), iy(:,:), iz(:,:)
  real(DP), allocatable :: jx(:,:), jy(:,:), jz(:,:)
  real(DP), allocatable :: unccn(:,:,:)
  real(DP), allocatable :: AOPclim_dNovrN  (:,:,:)
  real(DP)              :: mp_params(2)
  real(DP)              :: dt, ct

  !--- prognostic / diagnostic state (working copy) ----------------------------
  real(DP), allocatable :: rhog  (:,:,:), rhogvx(:,:,:), rhogvy(:,:,:)
  real(DP), allocatable :: rhogvz(:,:,:), rhogw (:,:,:), rhoge (:,:,:)
  real(DP), allocatable :: rhogq (:,:,:,:)
  real(DP), allocatable :: rho   (:,:,:), pre   (:,:,:), tem   (:,:,:)
  real(DP), allocatable :: vx    (:,:,:), vy    (:,:,:), vz    (:,:,:)
  real(DP), allocatable :: w     (:,:,:), q     (:,:,:,:)

  !--- pristine copy of the same state (restored before every repeat) ----------
  !--- MINIAPP: pristine copies of the per-step state. No counterpart in NICAM --
  !> the model advances once per timestep; the mini-app re-runs the SAME step
  !> num_repeat times for timing, so the state is restored from these before
  !> every iteration (timed separately as MAIN_Restore_state).
  real(DP), allocatable :: p_rhog  (:,:,:), p_rhogvx(:,:,:), p_rhogvy(:,:,:)
  real(DP), allocatable :: p_rhogvz(:,:,:), p_rhogw (:,:,:), p_rhoge (:,:,:)
  real(DP), allocatable :: p_rhogq (:,:,:,:)
  real(DP), allocatable :: p_rho   (:,:,:), p_pre   (:,:,:), p_tem   (:,:,:)
  real(DP), allocatable :: p_vx    (:,:,:), p_vy    (:,:,:), p_vz    (:,:,:)
  real(DP), allocatable :: p_w     (:,:,:), p_q     (:,:,:,:)

  !--- forcing that physics_step hands over as zero (CP/TB = 'NONE') -----------
  real(DP), allocatable :: CBMFX(:,:,:), qke(:,:,:)

  !--- mp_driver output --------------------------------------------------------
  real(DP), allocatable :: precip_mp(:,:,:), precip_ISO1_mp(:,:,:), precip_ISO2_mp(:,:,:)
  real(DP), allocatable :: precip_lh_mp(:,:), precip_rhoein_mp(:,:), precip_rhophi_mp(:,:), precip_rhokin_mp(:,:), precip_energy_mp(:,:)
  real(DP), allocatable :: cldfrac_mp(:,:,:), qhyd_mp(:,:,:)

  !--- validation scratch ------------------------------------------------------
  real(DP), allocatable :: refbuf(:)
  integer               :: nfail_total, ncheck
  real(DP)              :: worst_rel, worst_abs
  character(len=H_MID)  :: worst_name, worst_aname

  integer :: l, iter
  !=============================================================================

  !---< open parameter file >---
  call ADM_openparam(param_fname)

  !---< start process >---
  call ADM_proc_init

  !---< admin module setup >---
  call ADM_setup(param_fname,'msg_mp')

  !---< profiler module setup >---
  call PROF_setup

  !#############################################################################
  call PROF_setprefx('INIT')
  call PROF_rapstart('Initialize',0)

  write(IO_FID_LOG,*)
  write(IO_FID_LOG,*) '##### [NICAM-MP] start  setup     #####'
  if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
     write(*,*) '##### [NICAM-MP] start  setup     #####'
  endif

  call read_miniapp_param

  !---< constant module setup >---
  call CNST_setup

  !---< grid module setup (vertical axis only, see mod_grd) >---
  call GRD_setup

  !---< global configulation setup >---
  call runconf_setup

  !---< boundary condition setup >---
  ! NICAM calls this from dynamics_setup (mod_dynamics.f90:149).  The mini-app
  ! has no dynamics, but CNVVAR_prg2diag_in -> BNDCND_thermo needs it: without
  ! it rho/tem at k=kmin-1 and k=kmax+1 are left untouched.
  call BNDCND_setup

  !---< saturation process setup >---
  call SATURATION_setup

  !---< microphysics module setup >---
  ! physics_setup order (mod_physics.f90:319 and :349)
  call mp_init( MP_TYPE )

  !--- initialize vertical advection for OpenACC (mod_physics.f90:349).
  ! No-op on the CPU build; on the GPU build it creates the device-side
  ! kcell / dw work arrays that vadv1d_getflux dereferences.
  call vadv1d_init

  ijdim = ADM_gall_in
  kdim  = ADM_kall
  lall  = ADM_lall
  nqmax = TRC_VMAX

  call allocate_all

  call PROF_rapstart('Load_input',1)
  call load_input
  call PROF_rapend  ('Load_input',1)

  write(IO_FID_LOG,*)
  write(IO_FID_LOG,*) '##### [NICAM-MP] finish setup     #####'
  if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
     write(*,*) '##### [NICAM-MP] finish setup     #####'
  endif

  call flush(IO_FID_LOG)

  call PROF_rapend  ('Initialize',0)
  !#############################################################################
  call PROF_setprefx('MAIN')
  call PROF_rapstart('Main_Loop',0)

  write(IO_FID_LOG,*) '##### [NICAM-MP] start  main loop #####'
  if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
     write(*,*) '##### [NICAM-MP] start  main loop #####'
  endif

  !--- the persistent device-data region: mirrors the one prg_driver.f90 wraps
  !--- around its own time loop, restricted to what the MP chain reads.
  !$acc data &
  !$acc pcopyin(GRD_gz,GRD_gzh,GRD_dgz,GRD_dgzh,GRD_rdgz,GRD_rdgzh) &
  !$acc pcopyin(GRD_afact,GRD_bfact,GRD_cfact,GRD_dfact) &
  !$acc pcopyin(CPW,CVW) &
  !$acc pcopyin(z,zh,phi,gam2,gam2h,gsgam2,gsgam2h) &
  !$acc pcopyin(ix,iy,iz,jx,jy,jz,unccn,AOPclim_dNovrN,CBMFX,qke) &
  !$acc pcopy(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq) &
  !$acc pcopy(rho,pre,tem,vx,vy,vz,w,q) &
  !$acc pcopy(precip_mp,precip_ISO1_mp,precip_ISO2_mp) &
  !$acc pcopy(precip_lh_mp,precip_rhoein_mp,precip_rhophi_mp,precip_rhokin_mp,precip_energy_mp,cldfrac_mp,qhyd_mp)

  do iter = 1, num_repeat

     !--- every iteration must be the SAME work: mp_driver is INOUT on the whole
     !--- state, so the pristine copy is restored first.  The restore is counted
     !--- separately ('Restore_state') and is not part of 'MP_driver'.
     call PROF_rapstart('Restore_state',1)
     call restore_state
     call PROF_rapend  ('Restore_state',1)

     call PROF_rapstart('MP_driver',1)
     do l = 1, lall
        call mp_driver( ijdim,                           & ! [IN]
                        kdim,                            & ! [IN]
                        l,                               & ! [IN]
                        rhog            (:,:,l),         & ! [INOUT]
                        rhogvx          (:,:,l),         & ! [INOUT]
                        rhogvy          (:,:,l),         & ! [INOUT]
                        rhogvz          (:,:,l),         & ! [INOUT]
                        rhogw           (:,:,l),         & ! [INOUT]
                        rhoge           (:,:,l),         & ! [INOUT]
                        rhogq           (:,:,:,l),       & ! [INOUT]
                        rho             (:,:,l),         & ! [INOUT]
                        pre             (:,:,l),         & ! [INOUT]
                        tem             (:,:,l),         & ! [INOUT]
                        vx              (:,:,l),         & ! [INOUT]
                        vy              (:,:,l),         & ! [INOUT]
                        vz              (:,:,l),         & ! [INOUT]
                        w               (:,:,l),         & ! [INOUT]
                        q               (:,:,:,l),       & ! [INOUT]
                        precip_mp       (:,:,l),         & ! [OUT]
                        precip_ISO1_mp  (:,:,l),         & ! [OUT]
                        precip_ISO2_mp  (:,:,l),         & ! [OUT]
                        precip_lh_mp    (:,l),           & ! [OUT]
                        precip_rhoein_mp(:,l),           & ! [OUT]
                        precip_rhophi_mp(:,l),           & ! [OUT]
                        precip_rhokin_mp(:,l),           & ! [OUT]
                        precip_energy_mp(:,l),           & ! [OUT]
                        cldfrac_mp      (:,:,l),         & ! [INOUT]
                        qhyd_mp         (:,:,l),         & ! [OUT]
                        CBMFX           (:,:,l),         & ! [IN]
                        unccn           (:,:,l),         & ! [IN]
                        AOPclim_dNovrN  (:,ADM_KNONE,l), & ! [IN] factor_n_aie
                        qke             (:,:,l),         & ! [IN] rhogqke
                        z               (:,:,l),         & ! [IN]
                        zh              (:,:,l),         & ! [IN]
                        ix              (:,l),           & ! [IN]
                        iy              (:,l),           & ! [IN]
                        iz              (:,l),           & ! [IN]
                        jx              (:,l),           & ! [IN]
                        jy              (:,l),           & ! [IN]
                        jz              (:,l),           & ! [IN]
                        phi             (:,:,l),         & ! [IN]
                        gam2            (:,:,l),         & ! [IN]
                        gam2h           (:,:,l),         & ! [IN]
                        gsgam2          (:,:,l),         & ! [IN]
                        gsgam2h         (:,:,l),         & ! [IN]
                        dt,                              & ! [IN]
                        ct                               ) ! [IN]
     enddo
     !$acc wait !--- MINIAPP: drain async(1) so MP_driver includes ALL device work
     call PROF_rapend  ('MP_driver',1)

     write(IO_FID_LOG,'(1x,A,I6,A,I6,A)') '### mp_driver iteration ', iter, ' /', num_repeat, ' done'
     call flush(IO_FID_LOG)
  enddo

  !$acc update host(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq)
  !$acc update host(rho,pre,tem,vx,vy,vz,w,q)
  !$acc update host(precip_mp,precip_lh_mp,precip_rhoein_mp,precip_rhophi_mp,precip_rhokin_mp,precip_energy_mp,cldfrac_mp,qhyd_mp)

  !$acc end data
  !$acc wait

  write(IO_FID_LOG,*)
  write(IO_FID_LOG,*) '##### [NICAM-MP] finish main loop #####'
  if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
     write(*,*) '##### [NICAM-MP] finish main loop #####'
  endif

  call PROF_rapend('Main_Loop',0)
  !#############################################################################
  call PROF_setprefx('CHCK')

  if ( validate ) then
     call PROF_rapstart('Validate',0)
     call validate_all
     call PROF_rapend  ('Validate',0)
  endif

  if ( write_output ) then
     call PROF_rapstart('Write_output',0)
     call write_all
     call PROF_rapend  ('Write_output',0)
  endif

  call PROF_rapreport

  !--- finalize all process
  call ADM_proc_finish

  stop

contains

  !-----------------------------------------------------------------------------
  !> read the mini-app's own namelist
  subroutine read_miniapp_param
    implicit none

    integer :: ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[mp_miniapp]/Category[nhm driver]'
    rewind(IO_FID_CONF)
    read(IO_FID_CONF,nml=MP_MINIAPPPARAM,iostat=ierr)
    if ( ierr < 0 ) then
       write(IO_FID_LOG,*) '*** MP_MINIAPPPARAM is not specified. use default.'
    elseif( ierr > 0 ) then
       write(*         ,*) 'xxx Not appropriate names in namelist MP_MINIAPPPARAM. STOP.'
       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist MP_MINIAPPPARAM. STOP.'
       call ADM_proc_stop
    endif
    if( IO_FID_NML > 0 ) write(IO_FID_NML,nml=MP_MINIAPPPARAM)

    write(IO_FID_LOG,'(1x,A,A)'     ) '--- dataset directory   : ', trim(dataset_dir)
    write(IO_FID_LOG,'(1x,A,I10)'   ) '--- number of repeat    : ', num_repeat
    write(IO_FID_LOG,'(1x,A,L10)'   ) '--- validate?           : ', validate
    write(IO_FID_LOG,'(1x,A,ES15.8)') '--- relative tolerance  : ', tolerance
    write(IO_FID_LOG,'(1x,A,ES15.8)') '--- absolute tolerance  : ', tolerance_abs
    write(IO_FID_LOG,'(1x,A,L10)'   ) '--- write own output?   : ', write_output
    if ( write_output ) then
       write(IO_FID_LOG,'(1x,A,A)'  ) '--- output directory    : ', trim(output_dir)
    endif

    return
  end subroutine read_miniapp_param

  !-----------------------------------------------------------------------------
  subroutine allocate_all
    implicit none
    !---------------------------------------------------------------------------

    allocate( z      (ijdim,kdim,lall) )
    allocate( zh     (ijdim,kdim,lall) )
    allocate( phi    (ijdim,kdim,lall) )
    allocate( gam2   (ijdim,kdim,lall) )
    allocate( gam2h  (ijdim,kdim,lall) )
    allocate( gsgam2 (ijdim,kdim,lall) )
    allocate( gsgam2h(ijdim,kdim,lall) )
    allocate( ix(ijdim,lall), iy(ijdim,lall), iz(ijdim,lall) )
    allocate( jx(ijdim,lall), jy(ijdim,lall), jz(ijdim,lall) )
    allocate( unccn(ijdim,kdim,lall) )
    allocate( AOPclim_dNovrN  (ijdim,ADM_KNONE,lall) )

    allocate( rhog  (ijdim,kdim,lall) ); allocate( p_rhog  (ijdim,kdim,lall) )
    allocate( rhogvx(ijdim,kdim,lall) ); allocate( p_rhogvx(ijdim,kdim,lall) )
    allocate( rhogvy(ijdim,kdim,lall) ); allocate( p_rhogvy(ijdim,kdim,lall) )
    allocate( rhogvz(ijdim,kdim,lall) ); allocate( p_rhogvz(ijdim,kdim,lall) )
    allocate( rhogw (ijdim,kdim,lall) ); allocate( p_rhogw (ijdim,kdim,lall) )
    allocate( rhoge (ijdim,kdim,lall) ); allocate( p_rhoge (ijdim,kdim,lall) )
    allocate( rho   (ijdim,kdim,lall) ); allocate( p_rho   (ijdim,kdim,lall) )
    allocate( pre   (ijdim,kdim,lall) ); allocate( p_pre   (ijdim,kdim,lall) )
    allocate( tem   (ijdim,kdim,lall) ); allocate( p_tem   (ijdim,kdim,lall) )
    allocate( vx    (ijdim,kdim,lall) ); allocate( p_vx    (ijdim,kdim,lall) )
    allocate( vy    (ijdim,kdim,lall) ); allocate( p_vy    (ijdim,kdim,lall) )
    allocate( vz    (ijdim,kdim,lall) ); allocate( p_vz    (ijdim,kdim,lall) )
    allocate( w     (ijdim,kdim,lall) ); allocate( p_w     (ijdim,kdim,lall) )
    allocate( rhogq (ijdim,kdim,nqmax,lall) ); allocate( p_rhogq(ijdim,kdim,nqmax,lall) )
    allocate( q     (ijdim,kdim,nqmax,lall) ); allocate( p_q    (ijdim,kdim,nqmax,lall) )

    allocate( CBMFX(ijdim,kdim,lall) )
    allocate( qke  (ijdim,kdim,lall) )

    allocate( precip_mp   (ijdim,2,lall) )
    allocate( precip_ISO1_mp(ijdim,2,lall) )
    allocate( precip_ISO2_mp(ijdim,2,lall) )
    allocate( precip_lh_mp (ijdim,lall), precip_rhoein_mp(ijdim,lall) )
    allocate( precip_rhophi_mp(ijdim,lall), precip_rhokin_mp(ijdim,lall) )
    allocate( precip_energy_mp (ijdim,lall) )
    allocate( cldfrac_mp(ijdim,kdim,lall) )
    allocate( qhyd_mp (ijdim,kdim,lall) )

    allocate( refbuf(ijdim*kdim*nqmax*lall) )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'    ) '====== mini-app problem size ======'
    write(IO_FID_LOG,'(1x,A,I10)') '--- columns per region   (ijdim) : ', ijdim
    write(IO_FID_LOG,'(1x,A,I10)') '--- vertical levels      (kdim)  : ', kdim
    write(IO_FID_LOG,'(1x,A,I10)') '--- regions              (lall)  : ', lall
    write(IO_FID_LOG,'(1x,A,I10)') '--- tracers              (nqmax) : ', nqmax
    write(IO_FID_LOG,'(1x,A,I10)') '--- total columns                : ', ijdim*lall

    return
  end subroutine allocate_all

  !-----------------------------------------------------------------------------
  !> read one Fortran-sequential-unformatted big-endian record into a(n)
  subroutine read_record( name, a, n )
    implicit none

    character(len=*), intent(in)  :: name
    integer,          intent(in)  :: n
    real(DP),         intent(out) :: a(n)

    character(len=H_LONG) :: fname
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    fname = trim(dataset_dir)//'/'//trim(name)//'.bin'

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          status = 'old',         &
          iostat = ierr           )
    if ( ierr /= 0 ) then
       write(*         ,*) 'xxx Cannot open dump record: ', trim(fname)
       write(IO_FID_LOG,*) 'xxx Cannot open dump record: ', trim(fname)
       call ADM_proc_stop
    endif

    read(fid,iostat=ierr) a(:)
    if ( ierr /= 0 ) then
       write(*         ,*) 'xxx Record size mismatch: ', trim(fname), ' expected ', n
       write(IO_FID_LOG,*) 'xxx Record size mismatch: ', trim(fname), ' expected ', n
       call ADM_proc_stop
    endif
    close(fid)

    return
  end subroutine read_record

  !-----------------------------------------------------------------------------
  !> write one Fortran-sequential-unformatted big-endian record (same format)
  subroutine write_record( name, a, n )
    implicit none

    character(len=*), intent(in) :: name
    integer,          intent(in) :: n
    real(DP),         intent(in) :: a(n)

    character(len=H_LONG) :: fname
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    fname = trim(output_dir)//'/'//trim(name)//'.bin'

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          form   = 'unformatted', &
          status = 'replace',     &
          iostat = ierr           )
    if ( ierr /= 0 ) then
       write(*         ,*) 'xxx Cannot open output record: ', trim(fname)
       write(IO_FID_LOG,*) 'xxx Cannot open output record: ', trim(fname)
       call ADM_proc_stop
    endif

    write(fid) a(:)
    close(fid)

    return
  end subroutine write_record

  !-----------------------------------------------------------------------------
  !> the 16 static + 15 per-step records of INPUTS.txt sections 1-2
  subroutine load_input
    implicit none

    integer :: n3, n2, nq3
    !---------------------------------------------------------------------------

    n3  = ijdim*kdim*lall
    n2  = ijdim*lall
    nq3 = ijdim*kdim*nqmax*lall

    !--- section 1 : static (read once, identical every step)
    call read_record( 'mpF_params'     , mp_params    , 2   )
    dt = mp_params(1)
    ct = mp_params(2)

    call read_record( 'mpF_z'          , z      , n3 )
    call read_record( 'mpF_zh'         , zh     , n3 )
    call read_record( 'mpF_phi'        , phi    , n3 )
    call read_record( 'mpF_gam2'       , gam2   , n3 )
    call read_record( 'mpF_gam2h'      , gam2h  , n3 )
    call read_record( 'mpF_gsgam2'     , gsgam2 , n3 )
    call read_record( 'mpF_gsgam2h'    , gsgam2h, n3 )
    call read_record( 'mpF_ix'         , ix     , n2 )
    call read_record( 'mpF_iy'         , iy     , n2 )
    call read_record( 'mpF_iz'         , iz     , n2 )
    call read_record( 'mpF_jx'         , jx     , n2 )
    call read_record( 'mpF_jy'         , jy     , n2 )
    call read_record( 'mpF_jz'         , jz     , n2 )
    call read_record( 'mpF_unccn'      , unccn  , n3 )
    call read_record( 'mpF_aop_dnovrn' , AOPclim_dNovrN    , n2 )  ! (gall,KNONE,lall)

    !--- section 2 : per-step state, kept as the pristine copy
    call read_record( 'mpF_drv_rhog_in'  , p_rhog  , n3  )
    call read_record( 'mpF_drv_rhogvx_in', p_rhogvx, n3  )
    call read_record( 'mpF_drv_rhogvy_in', p_rhogvy, n3  )
    call read_record( 'mpF_drv_rhogvz_in', p_rhogvz, n3  )
    call read_record( 'mpF_drv_rhogw_in' , p_rhogw , n3  )
    call read_record( 'mpF_drv_rhoge_in' , p_rhoge , n3  )
    call read_record( 'mpF_drv_rhogq_in' , p_rhogq , nq3 )
    call read_record( 'mpF_drv_rho_in'   , p_rho   , n3  )
    call read_record( 'mpF_drv_pre_in'   , p_pre   , n3  )
    call read_record( 'mpF_drv_tem_in'   , p_tem   , n3  )
    call read_record( 'mpF_drv_vx_in'    , p_vx    , n3  )
    call read_record( 'mpF_drv_vy_in'    , p_vy    , n3  )
    call read_record( 'mpF_drv_vz_in'    , p_vz    , n3  )
    call read_record( 'mpF_drv_w_in'     , p_w     , n3  )
    call read_record( 'mpF_drv_q_in'     , p_q     , nq3 )

    !--- physics_step hands these over as zero when CP/TB = 'NONE'
    !--- (mod_physics.f90:1489-1503); GDCFRC enters as zero, GDCLW is [OUT].
    CBMFX(:,:,:) = 0.0_DP
    qke  (:,:,:) = 0.0_DP

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)'      ) '====== dumped step ======'
    write(IO_FID_LOG,'(1x,A,ES20.10)') '--- dt [s]                       : ', dt
    write(IO_FID_LOG,'(1x,A,ES20.10)') '--- ct [s]                       : ', ct

    return
  end subroutine load_input

  !-----------------------------------------------------------------------------
  !> put the state back to what the dump says it was before mp_driver
  subroutine restore_state
    implicit none
    !---------------------------------------------------------------------------

    rhog  (:,:,:)   = p_rhog  (:,:,:)
    rhogvx(:,:,:)   = p_rhogvx(:,:,:)
    rhogvy(:,:,:)   = p_rhogvy(:,:,:)
    rhogvz(:,:,:)   = p_rhogvz(:,:,:)
    rhogw (:,:,:)   = p_rhogw (:,:,:)
    rhoge (:,:,:)   = p_rhoge (:,:,:)
    rho   (:,:,:)   = p_rho   (:,:,:)
    pre   (:,:,:)   = p_pre   (:,:,:)
    tem   (:,:,:)   = p_tem   (:,:,:)
    vx    (:,:,:)   = p_vx    (:,:,:)
    vy    (:,:,:)   = p_vy    (:,:,:)
    vz    (:,:,:)   = p_vz    (:,:,:)
    w     (:,:,:)   = p_w     (:,:,:)
    rhogq (:,:,:,:) = p_rhogq (:,:,:,:)
    q     (:,:,:,:) = p_q     (:,:,:,:)

    !--- cldfrac_mp after cp_driver('NONE'): GDCFRC/GDCLW enter as zero
    cldfrac_mp(:,:,:) = 0.0_DP
    qhyd_mp (:,:,:) = 0.0_DP

    !$acc update device(rhog,rhogvx,rhogvy,rhogvz,rhogw,rhoge,rhogq)
    !$acc update device(rho,pre,tem,vx,vy,vz,w,q,cldfrac_mp,qhyd_mp)

    return
  end subroutine restore_state

  !-----------------------------------------------------------------------------
  !> compare one computed array with its dumped counterpart
  !> PASS metric identical to nicam-kks/tools/compare_dump.py :
  !>    pass  <=>  |d| <= atol  .OR.  |d| <= rtol*|ref|
  !> MAX_REL is |d|/|ref| over the elements with ref /= 0; the elements with
  !> ref == 0 and d /= 0 have no relative measure and are counted in N_ZREF
  !> (their |d| is already in MAX|DIFF|).
  !-----------------------------------------------------------------------------
  !--- MINIAPP: per-category rows for a (ijdim,kdim,nqmax,lall) tracer array.
  !> The Fortran carries every water category in ONE array indexed by the
  !> mod_runconf tracer index (I_QV..I_QG = 1..6 on the NSW6 path), so a single
  !> row hides WHICH category differs. Report one row per category, labelled
  !> with the model's own short name from TRC_name (qv, qc, qr, qi, qs, qg).
  !> The reference record is read once and sliced, not re-read per category.
  subroutine check_trc( name, label, a )
    implicit none

    character(len=*), intent(in) :: name   ! dump record name (see OUTPUTS.txt)
    character(len=*), intent(in) :: label  !--- MINIAPP: mp_driver dummy name
    real(DP),         intent(in) :: a(ijdim,kdim,nqmax,lall)

    real(DP) :: d, r, maxabs, maxrel
    integer  :: i, k, l, nq, idx, nbad, nelem
    character(len=6)  :: verdict
    character(len=24) :: rowname
    !---------------------------------------------------------------------------

    call read_record( name, refbuf, ijdim*kdim*nqmax*lall )

    nelem = ijdim*kdim*lall

    do nq = 1, nqmax
       maxabs = 0.0_DP
       maxrel = 0.0_DP
       nbad   = 0

       do l  = 1, lall
       do k  = 1, kdim
       do i  = 1, ijdim
          idx = i + ijdim*( (k-1) + kdim*( (nq-1) + nqmax*(l-1) ) )

          d = abs( a(i,k,nq,l) - refbuf(idx) )
          if ( d > maxabs ) maxabs = d

          if ( abs(refbuf(idx)) > CNST_EPS_DP ) then !--- MINIAPP: NICAM convention — |ref|<=EPS counts as zero (N_ZREF), no denormal division
             r = d / abs(refbuf(idx))
             if ( r > maxrel ) maxrel = r
          endif

          if ( .NOT. ( d <= tolerance_abs .OR. d <= tolerance*abs(refbuf(idx)) ) ) then
             nbad = nbad + 1
          endif
       enddo
       enddo
       enddo

       if ( nbad == 0 ) then
          verdict = 'OK    '
       else
          verdict = 'DIFFER'
       endif

       rowname = trim(label)//':'//trim(TRC_name(nq))

       write(IO_FID_LOG,'(1x,A24,1x,A6,2ES13.4,I12,I14)') &
            rowname, verdict, maxabs, maxrel, nbad, nelem
       if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
          write(*,'(1x,A24,1x,A6,2ES13.4,I12,I14)') &
               rowname, verdict, maxabs, maxrel, nbad, nelem
       endif

       ncheck      = ncheck + 1
       nfail_total = nfail_total + nbad
       if ( maxrel > worst_rel ) then
          worst_rel  = maxrel
          worst_name = rowname
       endif
       if ( maxabs > worst_abs ) then
          worst_abs   = maxabs
          worst_aname = rowname
       endif
    enddo

    return
  end subroutine check_trc

  !-----------------------------------------------------------------------------
  !--- MINIAPP: per-phase rows for the (ijdim,2,lall) surface precipitation
  !> flux. Index 1 = liquid (rain), 2 = solid (snow+graupel+ice) -- see the
  !> mp_driver [OUT] argument at mod_mp_driver.f90:830/:909.
  subroutine check_phase( name, label, a )
    implicit none

    character(len=*), intent(in) :: name   ! dump record name (see OUTPUTS.txt)
    character(len=*), intent(in) :: label  !--- MINIAPP: mp_driver dummy name
    real(DP),         intent(in) :: a(ijdim,2,lall)

    real(DP) :: d, r, maxabs, maxrel
    integer  :: i, l, ip, idx, nbad, nelem
    character(len=6)  :: verdict
    character(len=24) :: rowname
    character(len=3)  :: plabel(2)
    !---------------------------------------------------------------------------

    plabel(1) = 'liq'
    plabel(2) = 'sol'

    call read_record( name, refbuf, ijdim*2*lall )

    nelem = ijdim*lall

    do ip = 1, 2
       maxabs = 0.0_DP
       maxrel = 0.0_DP
       nbad   = 0

       do l = 1, lall
       do i = 1, ijdim
          idx = i + ijdim*( (ip-1) + 2*(l-1) )

          d = abs( a(i,ip,l) - refbuf(idx) )
          if ( d > maxabs ) maxabs = d

          if ( abs(refbuf(idx)) > CNST_EPS_DP ) then !--- MINIAPP: NICAM convention — |ref|<=EPS counts as zero (N_ZREF), no denormal division
             r = d / abs(refbuf(idx))
             if ( r > maxrel ) maxrel = r
          endif

          if ( .NOT. ( d <= tolerance_abs .OR. d <= tolerance*abs(refbuf(idx)) ) ) then
             nbad = nbad + 1
          endif
       enddo
       enddo

       if ( nbad == 0 ) then
          verdict = 'OK    '
       else
          verdict = 'DIFFER'
       endif

       rowname = trim(label)//':'//plabel(ip)

       write(IO_FID_LOG,'(1x,A24,1x,A6,2ES13.4,I12,I14)') &
            rowname, verdict, maxabs, maxrel, nbad, nelem
       if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
          write(*,'(1x,A24,1x,A6,2ES13.4,I12,I14)') &
               rowname, verdict, maxabs, maxrel, nbad, nelem
       endif

       ncheck      = ncheck + 1
       nfail_total = nfail_total + nbad
       if ( maxrel > worst_rel ) then
          worst_rel  = maxrel
          worst_name = rowname
       endif
       if ( maxabs > worst_abs ) then
          worst_abs   = maxabs
          worst_aname = rowname
       endif
    enddo

    return
  end subroutine check_phase

  subroutine check( name, label, a, n )
    implicit none

    character(len=*), intent(in) :: name   ! dump record name (see OUTPUTS.txt)
    character(len=*), intent(in) :: label  !--- MINIAPP: mp_driver dummy name
    integer,          intent(in) :: n
    real(DP),         intent(in) :: a(n)

    real(DP) :: d, r, maxabs, maxrel
    integer  :: i, nbad
    character(len=6)  :: verdict
    character(len=24) :: rowname
    !---------------------------------------------------------------------------

    call read_record( name, refbuf, n )

    maxabs = 0.0_DP
    maxrel = 0.0_DP
    nbad   = 0

    do i = 1, n
       d = abs( a(i) - refbuf(i) )
       if ( d > maxabs ) maxabs = d

       if ( abs(refbuf(i)) > CNST_EPS_DP ) then !--- MINIAPP: NICAM convention — |ref|<=EPS counts as zero (N_ZREF), no denormal division
          r = d / abs(refbuf(i))
          if ( r > maxrel ) maxrel = r
       endif

       if ( .NOT. ( d <= tolerance_abs .OR. d <= tolerance*abs(refbuf(i)) ) ) then
          nbad = nbad + 1
       endif
    enddo

    if ( nbad == 0 ) then
       verdict = 'OK    '
    else
       verdict = 'DIFFER'
    endif

    rowname = label

    write(IO_FID_LOG,'(1x,A24,1x,A6,2ES13.4,I12,I14)') &
         rowname, verdict, maxabs, maxrel, nbad, n
    if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
       write(*,'(1x,A24,1x,A6,2ES13.4,I12,I14)') &
            rowname, verdict, maxabs, maxrel, nbad, n
    endif

    ncheck      = ncheck + 1
    nfail_total = nfail_total + nbad
    if ( maxrel > worst_rel ) then
       worst_rel  = maxrel
       worst_name = rowname
    endif
    if ( maxabs > worst_abs ) then
       worst_abs   = maxabs
       worst_aname = rowname
    endif

    return
  end subroutine check

  !-----------------------------------------------------------------------------
  !> the 23 mp_driver entries of OUTPUTS.txt
  subroutine validate_all
    implicit none

    integer :: n3, n2, nq3, n2p
    !---------------------------------------------------------------------------

    n3  = ijdim*kdim*lall
    n2  = ijdim*lall
    n2p = ijdim*2*lall
    nq3 = ijdim*kdim*nqmax*lall

    ncheck      = 0
    nfail_total = 0
    worst_rel   = 0.0_DP
    worst_abs   = 0.0_DP
    worst_name  = 'none'
    worst_aname = 'none'

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A)') '====== VALIDATION : mp_driver (23 arrays, water categories listed separately) ======'
    write(IO_FID_LOG,'(1x,A)') '   ARRAY = the caller-side variable name used by' // &
                               ' mod_physics.f90 at its mp_driver call site;'
    write(IO_FID_LOG,'(1x,A)') '           dump record names are listed in OUTPUTS.txt'
    write(IO_FID_LOG,'(1x,A24,1x,A6,A13,A13,A12,A14)') &
         'ARRAY','STATUS','MAX|DIFF|','MAX_REL','N_FAIL','N_ELEM'
    if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
       write(*,*)
       write(*,'(1x,A)') '====== VALIDATION : mp_driver (23 arrays, water categories listed separately) ======'
       write(*,'(1x,A24,1x,A6,A13,A13,A12,A14)') &
            'ARRAY','STATUS','MAX|DIFF|','MAX_REL','N_FAIL','N_ELEM'
    endif

    call check( 'mpF_drv_rhog_out'     , 'rhog'         , rhog  , n3  )
    call check( 'mpF_drv_rhogvx_out'   , 'rhogvx'       , rhogvx, n3  )
    call check( 'mpF_drv_rhogvy_out'   , 'rhogvy'       , rhogvy, n3  )
    call check( 'mpF_drv_rhogvz_out'   , 'rhogvz'       , rhogvz, n3  )
    call check( 'mpF_drv_rhogw_out'    , 'rhogw'        , rhogw , n3  )
    call check( 'mpF_drv_rhoge_out'    , 'rhoge'        , rhoge , n3  )
    call check_trc( 'mpF_drv_rhogq_out'  , 'rhogq' , rhogq ) !--- MINIAPP: per category
    call check( 'mpF_drv_rho_out'      , 'rho'          , rho   , n3  )
    call check( 'mpF_drv_pre_out'      , 'pre'          , pre   , n3  )
    call check( 'mpF_drv_tem_out'      , 'tem'          , tem   , n3  )
    call check( 'mpF_drv_vx_out'       , 'vx'           , vx    , n3  )
    call check( 'mpF_drv_vy_out'       , 'vy'           , vy    , n3  )
    call check( 'mpF_drv_vz_out'       , 'vz'           , vz    , n3  )
    call check( 'mpF_drv_w_out'        , 'w'            , w     , n3  )
    call check_trc( 'mpF_drv_q_out'      , 'q'     , q     ) !--- MINIAPP: per category
    call check_phase( 'mpF_drv_precip'   , 'precip_mp', precip_mp ) !--- MINIAPP: liq / sol
    call check( 'mpF_drv_precip_lh'    , 'precip_lh_mp' , precip_lh_mp   , n2  )
    call check( 'mpF_drv_precip_rhoein', 'precip_rhoein_mp', precip_rhoein_mp  , n2  )
    call check( 'mpF_drv_precip_rhophi', 'precip_rhophi_mp', precip_rhophi_mp  , n2  )
    call check( 'mpF_drv_precip_rhokin', 'precip_rhokin_mp', precip_rhokin_mp  , n2  )
    call check( 'mpF_drv_precip_energy', 'precip_energy_mp', precip_energy_mp   , n2  )
    call check( 'mpF_drv_cldfrac'      , 'cldfrac_mp'   , cldfrac_mp, n3  )
    call check( 'mpF_drv_qhyd'         , 'qhyd_mp'      , qhyd_mp , n3  )

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,'(1x,A,I6)'        ) '--- rows checked              : ', ncheck
    write(IO_FID_LOG,'(1x,A,I14)'       ) '--- failing elements (total)  : ', nfail_total
    write(IO_FID_LOG,'(1x,A,ES13.4,A,A)') '--- worst relative difference : ', worst_rel, &
                                          '  in ', trim(worst_name)
    write(IO_FID_LOG,'(1x,A,ES13.4,A,A)') '--- worst absolute difference : ', worst_abs, &
                                          '  in ', trim(worst_aname)
    write(IO_FID_LOG,'(1x,A,ES13.4)'    ) '--- tolerance (relative)      : ', tolerance
    write(IO_FID_LOG,'(1x,A,ES13.4)'    ) '--- tolerance (absolute)      : ', tolerance_abs
    if ( nfail_total == 0 ) then
       write(IO_FID_LOG,'(1x,A)') '##### [NICAM-MP] VALIDATION PASS #####'
    else
       write(IO_FID_LOG,'(1x,A)') '##### [NICAM-MP] VALIDATION FAIL #####'
    endif

    if ( GLOBAL_prc_me == GLOBAL_prc_master ) then
       write(*,*)
       write(*,'(1x,A,I6)'        ) '--- rows checked              : ', ncheck
       write(*,'(1x,A,I14)'       ) '--- failing elements (total)  : ', nfail_total
       write(*,'(1x,A,ES13.4,A,A)') '--- worst relative difference : ', worst_rel, &
                                    '  in ', trim(worst_name)
       write(*,'(1x,A,ES13.4,A,A)') '--- worst absolute difference : ', worst_abs, &
                                    '  in ', trim(worst_aname)
       write(*,'(1x,A,ES13.4)'    ) '--- tolerance (relative)      : ', tolerance
       write(*,'(1x,A,ES13.4)'    ) '--- tolerance (absolute)      : ', tolerance_abs
       if ( nfail_total == 0 ) then
          write(*,'(1x,A)') '##### [NICAM-MP] VALIDATION PASS #####'
       else
          write(*,'(1x,A)') '##### [NICAM-MP] VALIDATION FAIL #####'
       endif
    endif

    return
  end subroutine validate_all

  !-----------------------------------------------------------------------------
  !> write our own copy of the 23 mp_driver outputs, same record format
  subroutine write_all
    implicit none

    integer :: n3, n2, nq3, n2p
    !---------------------------------------------------------------------------

    n3  = ijdim*kdim*lall
    n2  = ijdim*lall
    n2p = ijdim*2*lall
    nq3 = ijdim*kdim*nqmax*lall

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** write mini-app output records to ', trim(output_dir)

    call write_record( 'mpF_drv_rhog_out'      , rhog     , n3  )
    call write_record( 'mpF_drv_rhogvx_out'    , rhogvx   , n3  )
    call write_record( 'mpF_drv_rhogvy_out'    , rhogvy   , n3  )
    call write_record( 'mpF_drv_rhogvz_out'    , rhogvz   , n3  )
    call write_record( 'mpF_drv_rhogw_out'     , rhogw    , n3  )
    call write_record( 'mpF_drv_rhoge_out'     , rhoge    , n3  )
    call write_record( 'mpF_drv_rhogq_out'     , rhogq    , nq3 )
    call write_record( 'mpF_drv_rho_out'       , rho      , n3  )
    call write_record( 'mpF_drv_pre_out'       , pre      , n3  )
    call write_record( 'mpF_drv_tem_out'       , tem      , n3  )
    call write_record( 'mpF_drv_vx_out'        , vx       , n3  )
    call write_record( 'mpF_drv_vy_out'        , vy       , n3  )
    call write_record( 'mpF_drv_vz_out'        , vz       , n3  )
    call write_record( 'mpF_drv_w_out'         , w        , n3  )
    call write_record( 'mpF_drv_q_out'         , q        , nq3 )
    call write_record( 'mpF_drv_precip'        , precip_mp   , n2p )
    call write_record( 'mpF_drv_precip_lh'     , precip_lh_mp      , n2  )
    call write_record( 'mpF_drv_precip_rhoein' , precip_rhoein_mp     , n2  )
    call write_record( 'mpF_drv_precip_rhophi' , precip_rhophi_mp     , n2  )
    call write_record( 'mpF_drv_precip_rhokin' , precip_rhokin_mp     , n2  )
    call write_record( 'mpF_drv_precip_energy' , precip_energy_mp      , n2  )
    call write_record( 'mpF_drv_cldfrac'       , cldfrac_mp   , n3  )
    call write_record( 'mpF_drv_qhyd'          , qhyd_mp    , n3  )

    return
  end subroutine write_all

end program prg_mp_miniapp
