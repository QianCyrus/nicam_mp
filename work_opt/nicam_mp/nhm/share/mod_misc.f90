!-------------------------------------------------------------------------------
!> Module MISCELLANEOUS
!!
!! @par Description
!!         Miscellaneous functions
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_misc
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_stdio, only: &
     MISC_get_available_fid => IO_get_available_fid
  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: MISC_make_idstr
  public :: MISC_make_idstr2
  public :: MISC_get_available_fid
  public :: MISC_sec2ymdhms
  public :: MISC_ymdhms2sec
  public :: MISC_gammafunc
  public :: MISC_igammafunc
  public :: MISC_iigammafunc
  public :: MISC_scale_offset
  public :: MISC_make_fname_sepdir

  interface MISC_gammafunc
     module procedure MISC_gammafunc_SP
     module procedure MISC_gammafunc_DP
  end interface MISC_gammafunc

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters, variables & subroutines
  !
  !--- minimum available fid
  integer,parameter, private :: min_fid = 7
  !
  !--- maximum available fid.
  integer,parameter, private :: max_fid = 9999

  !--- number of separated file starts from 0 ?
  logical, private, parameter :: NSTR_ZERO_START = .true.

  !--- digit of separated file
  integer, private :: NSTR_MAX_DIGIT = 5
  !
  !<----  if running on ES,
  !<----      NSTR_ZERO_START = .TRUE.
  !<----      NSTR_MAX_DIGIT  = 5
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine MISC_make_idstr( &
       str,    &
       prefix, &
       ext,    &
       numID,  &
       digit   )
    implicit none

    character(len=*), intent(out) :: str    ! combined string (file name)
    character(len=*), intent(in)  :: prefix ! prefix
    character(len=*), intent(in)  :: ext    ! extention( e.g. .rgn )
    integer,          intent(in)  :: numID  ! number

    integer, optional, intent(in) :: digit  ! digit

    character(len=128) :: rankstr
    integer            :: setdigit
    !---------------------------------------------------------------------------

    if ( NSTR_ZERO_START ) then
       write(rankstr,'(I128.128)') numID-1
    else
       write(rankstr,'(I128.128)') numID
    endif

    if ( present(digit) ) then
       setdigit = digit
    else
       setdigit = NSTR_MAX_DIGIT
    endif

    rankstr(1:setdigit) = rankstr(128-(setdigit-1):128)
    rankstr(setdigit+1:128) = ' '

    str = trim(prefix)//'.'//trim(ext)//trim(rankstr) ! -> prefix.ext00000

    return
  end subroutine MISC_make_idstr

  !-----------------------------------------------------------------------------
  subroutine MISC_make_idstr2( & ! separate directory
       str,    &
       prefix, &
       ext,    &
       numID,  &
       digit   )
    implicit none

    character(len=*), intent(out) :: str    ! combined string (file name)
    character(len=*), intent(in)  :: prefix ! prefix
    character(len=*), intent(in)  :: ext    ! extention( e.g. .rgn )
    integer,          intent(in)  :: numID  ! number

    integer, optional, intent(in) :: digit  ! digit

    character(len=128) :: rankstr
    integer            :: setdigit
    character(len=5)   :: dirname
    integer :: i, idir
    !---------------------------------------------------------------------------

    if ( NSTR_ZERO_START ) then
       write(rankstr,'(I128.128)') numID-1
       idir = (numID-1)/1000
    else
       write(rankstr,'(I128.128)') numID
       idir = (numID)/1000
    endif

    if ( present(digit) ) then
       setdigit = digit
    else
       setdigit = NSTR_MAX_DIGIT
    endif

    rankstr(1:setdigit) = rankstr(128-(setdigit-1):128)
    rankstr(setdigit+1:128) = ' '

    write(dirname,'(I5.5)') idir ! can be allow 8 diget per 1000 file
    i = index(prefix,'/',back=.true.)

    if (i == 0 ) then
       str = './'//trim(dirname)//'/'//trim(prefix(i+1:))//'.'//trim(ext)//trim(rankstr) ! -> prefix.ext00000
    else
       str = trim(prefix(1:i-1))//'/'//trim(dirname)//'/'//trim(prefix(i+1:))//'.'//trim(ext)//trim(rankstr)  ! -> prefix.ext00000
    end if

    return
  end subroutine MISC_make_idstr2

  !-----------------------------------------------------------------------------
  subroutine MISC_sec2ymdhms( &
       ctime_second,     &
       year, month, day, &
       hour, min,   sec  )
    implicit none

    real(DP), intent(in)  :: ctime_second
    integer,  intent(out) :: year
    integer,  intent(out) :: month
    integer,  intent(out) :: day
    integer,  intent(out) :: hour
    integer,  intent(out) :: min
    integer,  intent(out) :: sec

    integer :: second
    !---------------------------------------------------------------------------

    second = idnint(ctime_second)
    year   = second / (60*60*24*30*12)

    second = mod(second,60*60*24*30*12)
    month  = second / (60*60*24*30)

    second = mod(second,60*60*24*30)
    day    = second / (60*60*24)

    second = mod(second,60*60*24)
    hour   = second / (60*60)

    second = mod(second,60*60)
    min    = second / (60)

    second = mod(second,60)
    sec    = second

    return
  end subroutine MISC_sec2ymdhms

  !-----------------------------------------------------------------------------
  subroutine MISC_ymdhms2sec( &
       ctime_second,     &
       year, month, day, &
       hour, min,   sec  )
    implicit none

    real(DP), intent(out) :: ctime_second
    integer,  intent(in)  :: sec
    integer,  intent(in)  :: min
    integer,  intent(in)  :: hour
    integer,  intent(in)  :: day
    integer,  intent(in)  :: month
    integer,  intent(in)  :: year
    !---------------------------------------------------------------------------

    ctime_second = real( year * 60*60*24*30*12 &
                       + month* 60*60*24*30    &
                       + day  * 60*60*24       &
                       + hour * 60*60          &
                       + min  * 60             &
                       + sec, kind=DP          )

    return
  end subroutine MISC_ymdhms2sec

  !-----------------------------------------------------------------------------
  !> Gamma function
  ! Lanczos approximation
  ! g = 7.0, N = 8
  function MISC_gammafunc_SP(x) result(f)
    implicit none

    real(SP), intent(in) :: x        !< input
    real(SP)             :: f

    real(DP), parameter :: LRTP = 0.9189385332046727_DP ! log(sqrt(2*PI))
    real(DP)            :: Ag, dx, work
    !---------------------------------------------------------------------------

    dx = real(x-1.0_SP,kind=DP)

    Ag = 0.999999999999809932276847004735_DP                     &
       + 676.520368121885098567009190444_DP    / ( dx + 1.0_DP ) &
       - 1259.13921672240287047156078755_DP    / ( dx + 2.0_DP ) &
       + 771.323428777653078848652825889_DP    / ( dx + 3.0_DP ) &
       - 176.615029162140599065845513540_DP    / ( dx + 4.0_DP ) &
       + 12.5073432786869048144589368533_DP    / ( dx + 5.0_DP ) &
       - 0.138571095265720116895547069851_DP   / ( dx + 6.0_DP ) &
       + 9.98436957801957085956266899569E-6_DP / ( dx + 7.0_DP ) &
       + 1.50563273514931155833835577539E-7_DP / ( dx + 8.0_DP )

    work = exp( LRTP + (dx+0.5_DP)*log(dx+7.5_DP) - (dx+7.5_DP) + log(Ag) )

    f = real(work,kind=SP)

    return
  end function MISC_gammafunc_SP

  !-----------------------------------------------------------------------------
  !> Gamma function
  ! Lanczos approximation
  ! g = 7.0, N = 8
  function MISC_gammafunc_DP(x) result(f)
    implicit none

    real(DP), intent(in) :: x        !< input
    real(DP)             :: f

    real(DP), parameter :: LRTP = 0.9189385332046727_DP ! log(sqrt(2*PI))
    real(DP)            :: Ag, dx, work
    !---------------------------------------------------------------------------

    dx = real(x-1.0_DP,kind=DP)

    Ag = 0.999999999999809932276847004735_DP                     &
       + 676.520368121885098567009190444_DP    / ( dx + 1.0_DP ) &
       - 1259.13921672240287047156078755_DP    / ( dx + 2.0_DP ) &
       + 771.323428777653078848652825889_DP    / ( dx + 3.0_DP ) &
       - 176.615029162140599065845513540_DP    / ( dx + 4.0_DP ) &
       + 12.5073432786869048144589368533_DP    / ( dx + 5.0_DP ) &
       - 0.138571095265720116895547069851_DP   / ( dx + 6.0_DP ) &
       + 9.98436957801957085956266899569E-6_DP / ( dx + 7.0_DP ) &
       + 1.50563273514931155833835577539E-7_DP / ( dx + 8.0_DP )

    work = exp( LRTP + (dx+0.5_DP)*log(dx+7.5_DP) - (dx+7.5_DP) + log(Ag) )

    f = real(work,kind=DP)

    return
  end function MISC_gammafunc_DP

  !-----------------------------------------------------------------------------
  function MISC_igammafunc( x, alpha, nm ) result(g)
    ! Section 6.2 in 'Numerical Recipes in C'
    ! function is represented by(6.2.1)
    ! g(x,alpha)=1/gamma(alpha)*integral_0^x { exp(-t) * t^(alpha-1) }dt
    ! g(0)=0 and g(+infinity)=1
    implicit none

    real(DP), intent(in) :: x
    real(DP), intent(in) :: alpha
    integer, intent(in), optional :: nm ! nmax
    real(DP) :: g
    ! work
    real(DP), parameter :: eps=1.d-30
    real(DP) :: a,b,c,d,h,an,del
    real(DP) :: gm, lgm, rn
    integer  :: nmax
    integer  :: n
    !---------------------------------------------------------------------------

    gm  = MISC_gammafunc(alpha)
    lgm = log(gm)
    if (x<alpha+1.d0) then ! series expansion (6.2.5)
       if (present(nm)) then
          nmax=nm
       else
          ! nmax=10 is enough to represent Droplet Size Distribution
          ! relative error is below 1% arround local maximum(x=alpha-1)
          nmax=10
       endif

       a = 1.d0/alpha
       g = a
       do n=1,nmax
          rn = real(n,kind=DP)
          a = a*x/(alpha+rn)
          g = g+a
       enddo
       g = g*exp( -x + alpha*log(x) - lgm )
    else ! continued fraction expansion (6.2.6)
       if (present(nm)) then
          nmax=nm
       else
          ! nmax=2 is enough to represent Droplet Size Distribution
          ! relative error is below 1% arround local maximum(x=alpha-1)
          nmax=2
       endif

       b = x+1.d0-alpha
       c = 1.d0/eps
       d = 1.d0/b
       h=d
       do n=1, nmax
          rn = real(n,kind=DP)
          an = -rn*(rn-alpha)
          b  = b + 2.d0
          d  = an*d+b
          if(abs(d)<eps) d=eps
          c  = b+an/c
          if(abs(c)<eps) c=eps
          d  = 1.d0/d
          del= d*c
          h  = h*del
       enddo
       g = 1.d0 - exp( -x+alpha*log(x) - lgm )*h
    endif

    return
  end function MISC_igammafunc

  !-----------------------------------------------------------------------------
  function MISC_iigammafunc( p, alpha ) result(x)
    ! p = MISC_igammafunc( x, alpha, nm )
    ! Didonato and Morris, 1986,
    ! Computation of the Incomplete Gamma Function Ratios and their Inverse,
    ! ACM Transactions on Mathmatical Software, Vol.12, No.4,
    implicit none

    real(DP), intent(in) :: alpha  !
    real(DP), intent(in) :: p
    real(DP) :: x

    real(DP) :: s,s2,s3,s4,s5
    real(DP) :: t,t2,t3,t4
    real(DP) :: tau
    real(DP) :: m
    real(DP) :: a_n, a_d
    real(DP) :: sa ! sqrt(alpha)

    real(DP), parameter :: a0=3.31125922108741d0
    real(DP), parameter :: a1=11.6616720288968d0
    real(DP), parameter :: a2=4.28342155967104d0
    real(DP), parameter :: a3=0.213623493715853d0
    real(DP), parameter :: b1=6.61053765625462d0
    real(DP), parameter :: b2=6.40691597760039d0
    real(DP), parameter :: b3=1.27364489782223d0
    real(DP), parameter :: b4=0.361170810188420d-1

    real(DP), parameter :: eps = 1.d-10
    !---------------------------------------------------------------------------

    if (alpha < 1.d0) then
       write(*,*) 'Invalid parameter, alpha=',alpha
       write(*,*) 'This routine does not support.'
       return
    endif
    if ( (alpha > 1.d0-eps) .AND. (alpha < 1.d0+eps) ) then ! alpha = 1
       x = -log(1.d0 - p)
       return
    endif

    if (p>=0.5d0) then
       m  = 1.d0
       tau= 1.d0-p
    else
       m  =-1.d0
       tau= p
    endif
    t    = sqrt(-2.d0*log(tau))
    t2   = t*t
    t3   = t2*t
    t4   = t2*t2
    a_n  = a0   + a1*t + a2*t2 + a3*t3
    a_d  = 1.d0 + b1*t + b2*t2 + b3*t3 + b4*t4
    s    = m*(t-a_n/a_d) !(32)

    s2   = s*s
    s3   = s*s2
    s4   = s2*s2
    s5   = s3*s2
    sa   = sqrt(alpha)

    x    = alpha + s*sa + (s2-1.d0)/3.d0 + (s3 - 7.d0*s)/(36.d0*sa) &
         - (3.d0*s4 + 7.d0*s2 - 16.d0)/(810.d0*alpha)               &
         + (9.d0*s5 + 256.d0*s3 - 433.d0*s)/(38880.d0*alpha*sa)

    return
  end function MISC_iigammafunc

  !-----------------------------------------------------------------------------
  subroutine MISC_scale_offset( &
       imax,   &
       jmax,   &
       var_r4, &
       scale,  &
       offset, &
       var_i2  )
    use mod_cnst, only: &
       CNST_UNDEF4, &
       CNST_UNDEF2
    implicit none

    integer,    intent(in)  :: imax
    integer,    intent(in)  :: jmax
    real(SP),   intent(in)  :: var_r4(imax,jmax)
    real(DP),   intent(out) :: scale
    real(DP),   intent(out) :: offset
    integer(2), intent(out) :: var_i2(imax,jmax)

    real(DP), parameter :: i2range = 65535.0_DP

    real(SP) :: var_min, var_max
    real(DP) :: range
    real(SP) :: minlimit, maxlimit

    integer  :: i, j
    !---------------------------------------------------------------------------

    ! search max/min
    var_min =  1.E+20_SP
    var_max = -1.E+20_SP
    do j = 1, jmax
    do i = 1, imax
       if ( var_r4(i,j) > CNST_UNDEF4 / 10.0_SP ) then ! neglect undefined value
          var_min = min( var_min, var_r4(i,j) )
          var_max = max( var_max, var_r4(i,j) )
       endif
    enddo
    enddo

    ! calculate scale factor and offset value
    offset   = real(var_min,kind=DP)
    range    = real(var_max,kind=DP) - offset
    scale    = range / i2range

    minlimit = real( offset - scale              , kind=SP )
    maxlimit = real( offset + scale * (i2range+1), kind=SP )

    write(*,*) '*** [MISC_scale_offset] scale-offset compression parameter'
    write(*,*) '*** Original      : min   = ', var_min,  ', max    = ', var_max
    write(*,*) '***               : undef = ', CNST_UNDEF4

    write(*,*) '*** Compression   : scale = ', scale,    ', offset = ', offset
    write(*,*) '*** + value limit : min     ', minlimit, ', max    = ', maxlimit
    write(*,*) '***               : undef = ', CNST_UNDEF2

    ! compress 4-byte floating point into 2-byte integer

    do j = 1, jmax
    do i = 1, imax

       if ( var_r4(i,j) > CNST_UNDEF4 / 10.0_SP ) then ! neglect undefined value

          if    ( var_r4(i,j) >= maxlimit ) then
             write(*,*) '*** [MISC_scale_offset] exceed max value limit. ', var_r4(i,j)
             stop
          elseif( var_r4(i,j) <= minlimit ) then
             write(*,*) '*** [MISC_scale_offset] exceed min value limit. ', var_r4(i,j)
             stop
          endif

          if ( scale == 0.0_DP ) then
             var_i2(i,j) = 0
          else
             var_i2(i,j) = nint( ( real(var_r4(i,j),kind=DP)-offset ) / scale, kind=2 )
          endif

       else
          var_i2(i,j) = CNST_UNDEF2
       endif

    enddo
    enddo

    return
  end subroutine MISC_scale_offset

  !-----------------------------------------------------------------------------
  subroutine MISC_make_fname_sepdir( & ! change filename with separate directory
       str,    &
       filename, &
       numID )
    implicit none
    ! from: filename = XXX/XXX.{rgn,pe}*
    ! to: filename = XXX/00000/XXX.{rgn,pe}*
    character(len=*), intent(out) :: str    ! combined string (file name)
    character(len=*), intent(in)  :: filename ! original filename
    integer,          intent(in)  :: numID  ! number

    character(len=128) :: rankstr
    character(len=5)   :: dirname
    integer :: i, idir
    !---------------------------------------------------------------------------

    if ( NSTR_ZERO_START ) then
       write(rankstr,'(I128.128)') numID-1
       idir = (numID-1)/1000
    else
       write(rankstr,'(I128.128)') numID
       idir = (numID)/1000
    endif

    write(dirname,'(I5.5)') idir ! can be allow 8 diget per 1000 file
    i = index(filename,'/',back=.true.)

    if (i == 0 ) then
       str = './'//trim(dirname)//'/'//trim(filename(i+1:))
    else
       str = trim(filename(1:i-1))//'/'//trim(dirname)//'/'//trim(filename(i+1:))
    end if

    return
  end subroutine MISC_make_fname_sepdir

end module mod_misc
