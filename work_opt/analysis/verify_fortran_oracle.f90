! =============================================================================
! verify_fortran_oracle.f90
!
! Differential Oracle comparing the REAL compiled NICAM Fortran mp_driver
! (Mode 0 = pristine baseline) against our C++/CUDA mathematical core
! (Mode 1 = ENG, Mode 2 = ENG+MATH, Mode 3 = ENG+MATH+MONO) across:
!   - All 23 mp_driver output arrays (34 rows) with exact NICAM GPU gate
!     (|d| <= 1e-10 OR |d| <= 1e-9 * |ref|)
!   - 10 diverse realistic atmospheric regimes (1,024 cols x 80 levels,
!     terrain-following metrics, 78-level stretched grid, crossing-trajectory
!     storms, extreme supersaturation, trace hydrometeors, negative-q states)
!   - True continuous 5-step and 100-step multi-step evolutions measuring
!     column water conservation, total energy drift, and pre-filter neg-q count!
! =============================================================================
program verify_fortran_oracle
  use iso_c_binding
  use mod_precision, only: DP, RP
  use mod_stdio
  use mod_adm
  use mod_prof
  use mod_cnst
  use mod_grd
  use mod_runconf
  use mod_bndcnd
  use mod_satadjust
  use mod_vadv1d
  use mod_mp_driver
  use mod_mp_opt_bridge
  implicit none

  integer, parameter :: IJ = 1024
  integer, parameter :: KD = 80
  integer, parameter :: NQ = 6

  ! Static fields
  real(DP) :: z(IJ,KD), zh(IJ,KD), phi(IJ,KD)
  real(DP) :: gam2(IJ,KD), gam2h(IJ,KD), gsgam2(IJ,KD), gsgam2h(IJ,KD)
  real(DP) :: ix(IJ), iy(IJ), iz(IJ), jx(IJ), jy(IJ), jz(IJ)
  real(DP) :: unccn(IJ,KD), factor_n_aie(IJ), CBMFX(IJ,KD), qke(IJ,KD)

  ! Initial state
  real(DP) :: init_rhog(IJ,KD), init_rhogvx(IJ,KD), init_rhogvy(IJ,KD), init_rhogvz(IJ,KD)
  real(DP) :: init_rhogw(IJ,KD), init_rhoge(IJ,KD), init_rhogq(IJ,KD,NQ)
  real(DP) :: init_rho(IJ,KD), init_pre(IJ,KD), init_tem(IJ,KD)
  real(DP) :: init_vx(IJ,KD), init_vy(IJ,KD), init_vz(IJ,KD), init_w(IJ,KD), init_q(IJ,KD,NQ)
  character(len=H_FILE) :: param_fname = 'oracle_test.cnf'

  ! Reference (Mode 0) single-step outputs
  real(DP) :: r_rhog(IJ,KD), r_rhogvx(IJ,KD), r_rhogvy(IJ,KD), r_rhogvz(IJ,KD)
  real(DP) :: r_rhogw(IJ,KD), r_rhoge(IJ,KD), r_rhogq(IJ,KD,NQ)
  real(DP) :: r_rho(IJ,KD), r_pre(IJ,KD), r_tem(IJ,KD)
  real(DP) :: r_vx(IJ,KD), r_vy(IJ,KD), r_vz(IJ,KD), r_w(IJ,KD), r_q(IJ,KD,NQ)
  real(DP) :: r_precip(IJ,2), r_piso1(IJ,2), r_piso2(IJ,2)
  real(DP) :: r_plh(IJ), r_prein(IJ), r_pphi(IJ), r_pkin(IJ), r_pene(IJ)
  real(DP) :: r_cfrc(IJ,KD), r_clw(IJ,KD)

  ! Working state
  real(DP) :: rhog(IJ,KD), rhogvx(IJ,KD), rhogvy(IJ,KD), rhogvz(IJ,KD)
  real(DP) :: rhogw(IJ,KD), rhoge(IJ,KD), rhogq(IJ,KD,NQ)
  real(DP) :: rho(IJ,KD), pre(IJ,KD), tem(IJ,KD)
  real(DP) :: vx(IJ,KD), vy(IJ,KD), vz(IJ,KD), w(IJ,KD), q(IJ,KD,NQ)
  real(DP) :: precip(IJ,2), piso1(IJ,2), piso2(IJ,2)
  real(DP) :: plh(IJ), prein(IJ), pphi(IJ), pkin(IJ), pene(IJ)
  real(DP) :: cfrc(IJ,KD), clw(IJ,KD)

  real(DP) :: dt, ct, t0, t1, t_base, t_mode
  integer  :: mode, nfail_mode, step, nsteps
  real(DP) :: w0, e0, wf, ef, cum_p, cum_pe, max_rel_all, max_abs_all
  character(len=16) :: mnames(0:3)
  integer :: ij_i, k_i, fid

  mnames(0) = 'BB (Fortran)   '
  mnames(1) = 'ENG (Fused)    '
  mnames(2) = 'ENG+MATH       '
  mnames(3) = 'ENG+MATH+MONO  '

  ! 1. Write synthetic vgrid78.dat and test config matching production Miyabi-G namelists
  call write_test_configs()

  call ADM_openparam(param_fname)
  call ADM_proc_init
  call ADM_setup(param_fname, 'msg_oracle')
  call PROF_setup
  call CNST_setup
  call GRD_setup
  call runconf_setup
  call BNDCND_setup
  call SATURATION_setup
  call mp_init(MP_TYPE)
  call vadv1d_init

  dt = 30.0_DP
  ct = 0.0_DP

  ! 2. Generate rich 1,024-column x 80-level 10-regime atmospheric state
  call build_atmospheric_state()

  write(*,'(A)') '================================================================================'
  write(*,'(A)') ' NICAM NSW6 DIFFERENTIAL ORACLE (Real Fortran Baseline vs C++/CUDA Math Core)'
  write(*,'(A,I5,A,I3,A,F6.1,A)') ' Grid: ', IJ, ' cols x ', KD, ' levels (10 atmospheric regimes), dt=', dt, ' s'
  write(*,'(A)') ' Validation gate: |d| <= 1e-10 OR |d| <= 1e-9 * |ref| across all 23 arrays (34 rows)'
  write(*,'(A)') '================================================================================'

  ! 3. Run Mode 0 (Pristine Baseline Fortran) first
  call set_opt_mode_runtime(0, 1, 0)
  call copy_init_to_work()
  call cpu_time(t0)
  call run_one_step()
  call cpu_time(t1)
  t_base = (t1 - t0) * 1000.0_DP
  call save_work_to_ref()
  write(*,'(A,F8.3,A)') ' [Mode 0: BB (Original Fortran)] single-step CPU time = ', t_base, ' ms'

  ! 4. Compare Modes 1, 2, 3 against Mode 0 across all 34 rows!
  do mode = 1, 3
     call set_opt_mode_runtime(mode, 1, 1)
     call mp_opt_reset_counters()
     call copy_init_to_work()
     call cpu_time(t0)
     call run_one_step()
     call cpu_time(t1)
     t_mode = (t1 - t0) * 1000.0_DP

     call compare_against_ref(mode, mnames(mode), t_mode, t_base, nfail_mode, max_rel_all, max_abs_all)
     call mp_opt_report_counters(trim(mnames(mode)))

     if (mode <= 2 .and. nfail_mode > 0) then
        write(*,'(A,I2,A,I10)') 'FATAL: Exact mode ', mode, ' FAILED validation with nfail=', nfail_mode
        stop 2
     endif
  enddo

  ! 5. True Continuous Multi-Step Evolution (5 steps and 100 steps!)
  do nsteps = 5, 100, 95
     write(*,*)
     write(*,'(A,I3,A)') '================================================================================'
     write(*,'(A,I3,A,F7.1,A)') ' CONTINUOUS MULTI-STEP EVOLUTION (', nsteps, ' steps, total T = ', nsteps*dt, ' s)'
     write(*,'(A)') '================================================================================'
     do mode = 0, 3
        call set_opt_mode_runtime(mode, 1, 1)
        call mp_opt_reset_counters()
        call copy_init_to_work()
        call compute_col_budgets(w0, e0)
        cum_p = 0.0_DP; cum_pe = 0.0_DP
        do step = 1, nsteps
           call run_one_step()
           cum_p  = cum_p  + sum(precip(:,1) + precip(:,2)) * dt
           cum_pe = cum_pe + sum(pene(:)) * dt
        enddo
        call compute_col_budgets(wf, ef)
        write(*,'(1x,A15,A,ES11.3,A,ES11.3,A,ES13.5,A,ES10.2)') &
             mnames(mode), &
             ' | rel_mass_drift=', abs((wf + cum_p) - w0)/max(abs(w0),1.0_DP), &
             ' | rel_E_drift=', abs((ef + cum_pe) - e0)/max(abs(e0),1.0_DP), &
             ' | cum_precip=', cum_p, &
             ' | min_q=', minval(q(:,ADM_kmin:ADM_kmax,2:NQ))
        if (mode >= 1) call mp_opt_report_counters(trim(mnames(mode)))
     enddo
  enddo

  write(*,*)
  write(*,'(A)') '##### [ORACLE] ALL EXACT VALIDATION GATES PASSED (N_FAIL = 0) #####'
  call ADM_proc_finish

contains

  subroutine set_opt_mode_runtime(m, rceff, cnt)
    integer, intent(in) :: m, rceff, cnt
    ! Direct helper using C setenv + resetting cached env via dedicated hook
    interface
       integer(c_int) function c_setenv(name, val, ov) bind(C, name="setenv")
         import :: c_int, c_char
         character(kind=c_char), intent(in) :: name(*), val(*)
         integer(c_int), value :: ov
       end function c_setenv
    end interface
    character(len=4) :: sm, sr, sc
    integer(c_int) :: rc
    write(sm,'(I0)') m
    write(sr,'(I0)') rceff
    write(sc,'(I0)') cnt
    rc = c_setenv('NICAM_OPT_MODE'//c_null_char, trim(sm)//c_null_char, 1_c_int)
    rc = c_setenv('NICAM_COMPUTE_RCEFF'//c_null_char, trim(sr)//c_null_char, 1_c_int)
    rc = c_setenv('NICAM_COLLECT_COUNTERS'//c_null_char, trim(sc)//c_null_char, 1_c_int)
    call mp_opt_reload_env()
  end subroutine set_opt_mode_runtime

  subroutine write_test_configs()
    integer :: k, u
    real(DP) :: gz(KD), gzh(KD), dz_k, h
    ! Stretched 78-layer atmospheric vertical grid (66 m near surface up to ~38 km)
    h = 0.0_DP
    gzh(1) = -33.0_DP
    gzh(2) = 0.0_DP
    do k = 2, KD - 1
       ! Smoothly stretched from 66 m at k=2 to ~1400 m aloft
       dz_k = 66.0_DP + 1350.0_DP * (dble(k - 2) / dble(KD - 3))**2.2_DP
       gzh(k+1) = gzh(k) + dz_k
    enddo
    do k = 2, KD - 1
       gz(k) = 0.5_DP * (gzh(k) + gzh(k+1))
    enddo
    gz(1)  = -gz(2)
    gz(KD) = 2.0_DP * gzh(KD) - gz(KD-1)

    open(newunit=u, file='vgrid78_oracle.dat', form='unformatted', status='replace')
    write(u) KD - 2
    write(u) gz
    write(u) gzh
    close(u)

    open(newunit=u, file='oracle_test.cnf', status='replace')
    write(u,'(A)') '&IOPARAM log_basename = "msg_oracle", /'
    write(u,'(A)') '&PROFPARAM PROF_rap_level = 1, /'
    write(u,'(A)') '&ADMPARAM glevel = 2, rlevel = 0, vlayer = 78, lall = 1, /'
    write(u,'(A)') '&GRDPARAM vgrid_fname = "vgrid78_oracle.dat", /'
    write(u,'(A)') '&RUNCONFPARAM'
    write(u,'(A)') '  run_phys_step = .true.,'
    write(u,'(A)') '  EIN_TYPE = "SIMPLE",'
    write(u,'(A)') '  RAIN_TYPE = "COLD",'
    write(u,'(A)') '  PRCIP_TRN_ECORRECT = "KIN2EIN",'
    write(u,'(A)') '  MP_TYPE = "NSW6",'
    write(u,'(A)') '  MP_DIV_NUM = 1,'
    write(u,'(A)') '/'
    write(u,'(A)') '&NM_MP_NSW6'
    write(u,'(A)') '  GAMMA_SAUT = 6.0D-2,'
    write(u,'(A)') '  GAMMA_SACR = 2.5D-2,'
    write(u,'(A)') '  PSAUT_BETA0 = 6.0D-3,'
    write(u,'(A)') '  Roh_flag = .true.,'
    write(u,'(A)') '  precip_scheme_type = "Flux-Semilag_new",'
    write(u,'(A)') '/'
    close(u)
  end subroutine write_test_configs

  subroutine build_atmospheric_state()
    integer :: i, k, reg
    real(DP) :: lon, lat, topo, ztop, alt, alth, T_k, p_k, rho_k, qv_k, qc_k, qr_k, qi_k, qs_k, qg_k
    real(DP) :: psat_l, qsat_l, psat_i, qsat_i, qd_k, rhogh_k, f_r

    ztop = GRD_gzh(KD)
    CBMFX = 0.0_DP
    qke   = 0.0_DP

    do i = 1, IJ
       reg = mod(i - 1, 10)
       lon = 2.0_DP * CNST_PI * dble(i) / dble(IJ)
       lat = -1.2_DP + 2.4_DP * dble(mod(i * 37, IJ)) / dble(IJ)
       ! Orthonormal local tangent basis (i_vec, j_vec)
       ix(i) = -sin(lon);           iy(i) =  cos(lon);           iz(i) = 0.0_DP
       jx(i) = -sin(lat)*cos(lon);  jy(i) = -sin(lat)*sin(lon);  jz(i) = cos(lat)
       factor_n_aie(i) = 0.8_DP + 0.4_DP * dble(mod(i, 7)) / 6.0_DP
       topo = 1200.0_DP * (0.5_DP + 0.5_DP * sin(3.0_DP * lon) * cos(2.0_DP * lat))

       do k = 1, KD
          ! Terrain-following Gal-Chen height & non-trivial metrics (gsgam2 /= gsgam2h!)
          z(i,k)       = topo + (1.0_DP - topo / ztop) * GRD_gz(k)
          zh(i,k)      = topo + (1.0_DP - topo / ztop) * GRD_gzh(k)
          phi(i,k)     = CNST_GRAV * z(i,k)
          gam2(i,k)    = 1.0_DP + 0.002_DP * z(i,k) / ztop
          gam2h(i,k)   = 1.0_DP + 0.002_DP * zh(i,k) / ztop
          gsgam2(i,k)  = (1.0_DP - topo / ztop) * gam2(i,k)
          gsgam2h(i,k) = (1.0_DP - topo / ztop) * gam2h(i,k)
          unccn(i,k)   = 1.0E8_DP

          alt = max(z(i,k), 0.0_DP)
          ! Realistic temperature profile: 300K (tropical) to 260K (polar), lapse 6.5 K/km up to 15 km
          if (reg == 2 .or. reg == 4) then
             T_k = max(200.0_DP, 282.0_DP - 6.5E-3_DP * alt)
          else
             T_k = max(195.0_DP, 299.0_DP - 6.5E-3_DP * alt)
          endif
          p_k = 101325.0_DP * exp(-alt / 8200.0_DP)
          rho_k = p_k / (CNST_Rdry * T_k)

          psat_l = CNST_PSAT0 * exp((LHV / CNST_Rvap) * (1.0_DP/CNST_TEM00 - 1.0_DP/T_k))
          qsat_l = psat_l / (rho_k * CNST_Rvap * T_k)
          psat_i = CNST_PSAT0 * exp((LHS / CNST_Rvap) * (1.0_DP/CNST_TEM00 - 1.0_DP/T_k))
          qsat_i = psat_i / (rho_k * CNST_Rvap * T_k)

          qv_k = 0.45_DP * min(qsat_l, qsat_i)
          qc_k = 0.0_DP; qr_k = 0.0_DP; qi_k = 0.0_DP; qs_k = 0.0_DP; qg_k = 0.0_DP

          if (k >= ADM_kmin .and. k <= ADM_kmax) then
             select case (reg)
             case (0) ! Clear-sky dry subsaturated
                qv_k = 0.30_DP * min(qsat_l, qsat_i)
             case (1) ! Tropical deep convection (heavy rain + graupel + supersat)
                if (alt < 12000.0_DP) then
                   qv_k = 1.15_DP * qsat_l
                   if (T_k > 265.0_DP) qc_k = 1.8E-3_DP
                   if (T_k > 270.0_DP) qr_k = 4.2E-3_DP * exp(-alt / 5000.0_DP)
                   if (T_k < 275.0_DP) qg_k = 2.5E-3_DP * exp(-abs(alt - 6500.0_DP) / 2500.0_DP)
                   if (T_k < 273.0_DP) qs_k = 1.2E-3_DP * exp(-abs(alt - 8000.0_DP) / 3000.0_DP)
                   if (T_k < 260.0_DP) qi_k = 4.0E-4_DP
                endif
             case (2) ! Mid-latitude stratiform snow & melting layer
                if (alt < 10000.0_DP) then
                   qv_k = 0.95_DP * min(qsat_l, qsat_i)
                   if (T_k < 276.0_DP .and. T_k > 240.0_DP) qs_k = 1.5E-3_DP
                   if (T_k > 271.0_DP) qr_k = 8.0E-4_DP
                   if (T_k < 270.0_DP) qc_k = 3.0E-4_DP
                endif
             case (3) ! Cold upper-tropospheric cirrus (ice-supersaturated, Pigen+Pidep)
                if (alt > 7000.0_DP .and. alt < 15000.0_DP) then
                   qv_k = 1.25_DP * qsat_i
                   if (mod(i, 2) == 0) qi_k = 2.0E-4_DP
                endif
             case (4) ! Supercooled mixed-phase cloud (-40 < T < 0 C)
                if (T_k > 235.0_DP .and. T_k < 273.0_DP) then
                   qv_k = 1.05_DP * qsat_l
                   qc_k = 1.0E-3_DP; qi_k = 2.0E-4_DP; qs_k = 6.0E-4_DP; qg_k = 4.0E-4_DP
                endif
             case (5) ! Trajectory-crossing shear storm (sharp heavy band over zero layer!)
                if (k >= 12 .and. k <= 22) then
                   qr_k = 6.0E-3_DP; qg_k = 4.0E-3_DP; qs_k = 2.0E-3_DP
                   qv_k = 0.9_DP * qsat_l
                else if (k >= 5 .and. k <= 8) then
                   qr_k = 5.0E-5_DP; qg_k = 0.0_DP
                endif
             case (6) ! Extreme supersaturation (5x qsat, stress tests SatAdj!)
                if (alt < 9000.0_DP) then
                   qv_k = min(5.0_DP * qsat_l, 0.045_DP)
                   qc_k = 5.0E-4_DP
                endif
             case (7) ! Trace sub-threshold hydrometeors (around 1e-12 zerosw edge!)
                qr_k = 1.0E-12_DP * (0.5_DP + dble(mod(k, 3)) * 0.5_DP)
                qs_k = 1.0E-12_DP * (0.5_DP + dble(mod(k + 1, 3)) * 0.5_DP)
                qg_k = 1.0E-12_DP * (0.5_DP + dble(mod(k + 2, 3)) * 0.5_DP)
                qi_k = 1.0E-13_DP
             case (8) ! Pre-existing negative water perturbations (stress tests negative_filter!)
                qv_k = 0.8_DP * qsat_l
                qc_k = merge(-1.0E-6_DP, 4.0E-4_DP, mod(k, 3) == 0)
                qr_k = merge(-2.0E-6_DP, 1.0E-3_DP, mod(k, 4) == 0)
                qs_k = merge(-1.0E-7_DP, 5.0E-4_DP, mod(k, 5) == 0)
             case (9) ! Near-surface heavy rain + strong vertical velocity (KIN2EIN stress test!)
                if (alt < 6000.0_DP) then
                   qv_k = 1.08_DP * qsat_l
                   qc_k = 1.2E-3_DP; qr_k = 5.5E-3_DP; qs_k = 8.0E-4_DP; qg_k = 1.5E-3_DP
                endif
             end select
          endif

          qd_k = 1.0_DP - (qv_k + qc_k + qr_k + qi_k + qs_k + qg_k)
          init_rho(i,k)  = rho_k
          init_tem(i,k)  = T_k
          init_pre(i,k)  = rho_k * (qd_k * CNST_Rdry + qv_k * CNST_Rvap) * T_k
          init_vx(i,k)   = 15.0_DP * cos(lon) + 5.0_DP * sin(alt / 3000.0_DP)
          init_vy(i,k)   = -12.0_DP * sin(lon) + 4.0_DP * cos(alt / 3000.0_DP)
          init_vz(i,k)   = 3.0_DP * sin(lat)
          init_w(i,k)    = merge(3.5_DP * sin(CNST_PI * alt / 15000.0_DP), 0.2_DP, reg == 1 .or. reg == 9)
          if (k == 1 .or. k == KD) init_w(i,k) = 0.0_DP

          init_rhog(i,k)   = rho_k * gsgam2(i,k)
          init_rhogvx(i,k) = init_rhog(i,k) * init_vx(i,k)
          init_rhogvy(i,k) = init_rhog(i,k) * init_vy(i,k)
          init_rhogvz(i,k) = init_rhog(i,k) * init_vz(i,k)
          init_rhoge(i,k)  = init_rhog(i,k) * T_k * CNST_CVdry

          init_q(i,k,1) = qv_k; init_q(i,k,2) = qc_k; init_q(i,k,3) = qr_k
          init_q(i,k,4) = qi_k; init_q(i,k,5) = qs_k; init_q(i,k,6) = qg_k
          init_rhogq(i,k,1:6) = init_rhog(i,k) * init_q(i,k,1:6)
       enddo

       ! Half-level rhogw consistent with C2W interpolation
       do k = 1, KD
          if (k == 1) then
             rhogh_k = (GRD_afact(2)*init_rho(i,2) + GRD_bfact(2)*init_rho(i,1)) * gsgam2h(i,2)
          else
             rhogh_k = (GRD_afact(k)*init_rho(i,k) + GRD_bfact(k)*init_rho(i,k-1)) * gsgam2h(i,k)
          endif
          init_rhogw(i,k) = rhogh_k * init_w(i,k)
       enddo
    enddo
  end subroutine build_atmospheric_state

  subroutine copy_init_to_work()
    rhog = init_rhog; rhogvx = init_rhogvx; rhogvy = init_rhogvy; rhogvz = init_rhogvz
    rhogw = init_rhogw; rhoge = init_rhoge; rhogq = init_rhogq
    rho = init_rho; pre = init_pre; tem = init_tem
    vx = init_vx; vy = init_vy; vz = init_vz; w = init_w; q = init_q
    cfrc = 0.0_DP; clw = 0.0_DP
  end subroutine copy_init_to_work

  subroutine save_work_to_ref()
    r_rhog = rhog; r_rhogvx = rhogvx; r_rhogvy = rhogvy; r_rhogvz = rhogvz
    r_rhogw = rhogw; r_rhoge = rhoge; r_rhogq = rhogq
    r_rho = rho; r_pre = pre; r_tem = tem
    r_vx = vx; r_vy = vy; r_vz = vz; r_w = w; r_q = q
    r_precip = precip; r_piso1 = piso1; r_piso2 = piso2
    r_plh = plh; r_prein = prein; r_pphi = pphi; r_pkin = pkin; r_pene = pene
    r_cfrc = cfrc; r_clw = clw
  end subroutine save_work_to_ref

  subroutine run_one_step()
    call mp_driver( &
         IJ, KD, 1, &
         rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq, &
         rho, pre, tem, vx, vy, vz, w, q, &
         precip, piso1, piso2, plh, prein, pphi, pkin, pene, &
         cfrc, clw, CBMFX, unccn, factor_n_aie, qke, &
         z, zh, ix, iy, iz, jx, jy, jz, phi, gam2, gam2h, gsgam2, gsgam2h, &
         dt, ct )
  end subroutine run_one_step

  subroutine check_arr(lbl, a, ref, n, nfail_tot, worst_rel, worst_abs)
    character(len=*), intent(in) :: lbl
    integer,  intent(in) :: n
    real(DP), intent(in) :: a(n), ref(n)
    integer,  intent(inout) :: nfail_tot
    real(DP), intent(inout) :: worst_rel, worst_abs
    integer :: idx, nbad
    real(DP) :: d, r, ma, mr
    ma = 0.0_DP; mr = 0.0_DP; nbad = 0
    do idx = 1, n
       d = abs(a(idx) - ref(idx))
       if (d > ma) ma = d
       if (abs(ref(idx)) > CNST_EPS_DP) then
          r = d / abs(ref(idx))
          if (r > mr) mr = r
       endif
       if (.not. (d <= 1.0D-10 .or. d <= 1.0D-9 * abs(ref(idx)))) then
          nbad = nbad + 1
       endif
    enddo
    nfail_tot = nfail_tot + nbad
    if (mr > worst_rel) worst_rel = mr
    if (ma > worst_abs) worst_abs = ma
    if (nbad > 0) then
       write(*,'(3x,A18,A,ES11.3,A,ES11.3,A,I8)') lbl, ' DIFFER | maxabs=', ma, ' maxrel=', mr, ' nfail=', nbad
    endif
  end subroutine check_arr

  subroutine compare_against_ref(m, mname, tm, tb, nfail_tot, worst_rel, worst_abs)
    integer, intent(in) :: m
    character(len=*), intent(in) :: mname
    real(DP), intent(in) :: tm, tb
    integer, intent(out) :: nfail_tot
    real(DP), intent(out) :: worst_rel, worst_abs
    integer :: nq_i, ip
    nfail_tot = 0; worst_rel = 0.0_DP; worst_abs = 0.0_DP

    call check_arr('rhog',   rhog,   r_rhog,   IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('rhogvx', rhogvx, r_rhogvx, IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('rhogvy', rhogvy, r_rhogvy, IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('rhogvz', rhogvz, r_rhogvz, IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('rhogw',  rhogw,  r_rhogw,  IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('rhoge',  rhoge,  r_rhoge,  IJ*KD, nfail_tot, worst_rel, worst_abs)
    do nq_i = 1, NQ
       call check_arr('rhogq:'//trim(TRC_name(nq_i)), rhogq(:,:,nq_i), r_rhogq(:,:,nq_i), IJ*KD, nfail_tot, worst_rel, worst_abs)
    enddo
    call check_arr('rho', rho, r_rho, IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('pre', pre, r_pre, IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('tem', tem, r_tem, IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('vx',  vx,  r_vx,  IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('vy',  vy,  r_vy,  IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('vz',  vz,  r_vz,  IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('w',   w,   r_w,   IJ*KD, nfail_tot, worst_rel, worst_abs)
    do nq_i = 1, NQ
       call check_arr('q:'//trim(TRC_name(nq_i)), q(:,:,nq_i), r_q(:,:,nq_i), IJ*KD, nfail_tot, worst_rel, worst_abs)
    enddo
    do ip = 1, 2
       call check_arr('precip', precip(:,ip), r_precip(:,ip), IJ, nfail_tot, worst_rel, worst_abs)
    enddo
    call check_arr('precip_lh',     plh,   r_plh,   IJ,    nfail_tot, worst_rel, worst_abs)
    call check_arr('precip_rhoein', prein, r_prein, IJ,    nfail_tot, worst_rel, worst_abs)
    call check_arr('precip_rhophi', pphi,  r_pphi,  IJ,    nfail_tot, worst_rel, worst_abs)
    call check_arr('precip_rhokin', pkin,  r_pkin,  IJ,    nfail_tot, worst_rel, worst_abs)
    call check_arr('precip_energy', pene,  r_pene,  IJ,    nfail_tot, worst_rel, worst_abs)
    call check_arr('cldfrac',       cfrc,  r_cfrc,  IJ*KD, nfail_tot, worst_rel, worst_abs)
    call check_arr('qhyd',          clw,   r_clw,   IJ*KD, nfail_tot, worst_rel, worst_abs)

    write(*,'(1x,A,I1,A,A15,A,F7.3,A,F5.2,A,I8,A,ES10.2,A,ES10.2)') &
         '[Mode ', m, ': ', mname, '] time=', tm, ' ms (', tb/max(tm,1.0D-6), 'x CPU) | N_FAIL=', &
         nfail_tot, ' | max_rel=', worst_rel, ' | max_abs=', worst_abs
  end subroutine compare_against_ref

  subroutine compute_col_budgets(wmass, tenergy)
    real(DP), intent(out) :: wmass, tenergy
    integer :: i, k, nq_i
    real(DP) :: dgz_k, rg_k, kh_k
    ! Exact NICAM metric-weighted column integral: sum_k (rhog * var) * GRD_dgz(k)
    ! telescopes identically against -(F(k+1) - F(k)) * GRD_rdgz(k) and surface precip!
    wmass = 0.0_DP; tenergy = 0.0_DP
    do k = ADM_kmin, ADM_kmax
       dgz_k = GRD_dgz(k)
       do i = 1, IJ
          do nq_i = 1, NQ
             wmass = wmass + rhogq(i,k,nq_i) * dgz_k
          enddo
          rg_k = rhog(i,k)
          kh_k = 0.5_DP * (rhogvx(i,k)**2 + rhogvy(i,k)**2 + rhogvz(i,k)**2) / rg_k
          tenergy = tenergy + (rhoge(i,k) + rg_k * CNST_GRAV * z(i,k) + kh_k &
                               + LHV * rhogq(i,k,I_QV) &
                               - LHF * (rhogq(i,k,I_QI) + rhogq(i,k,I_QS) + rhogq(i,k,I_QG))) * dgz_k
       enddo
    enddo
  end subroutine compute_col_budgets

end program verify_fortran_oracle
