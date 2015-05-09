!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -*- Mode: F90 -*- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! module_hydro_commun.f90 --- 
!!!!
!! module hydro_precision
!! module hydro_commons
!! module hydro_parameters 
!! module hydro_const
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module hydro_precision
  integer, parameter :: prec_real=kind(1.d0)
  integer, parameter :: prec_int=4
  integer, parameter :: prec_output=4
end module hydro_precision

module hydro_commons
  use hydro_precision
  integer(kind=prec_int) :: imin,imax,jmin,jmax
  real(kind=prec_real),allocatable,dimension(:,:,:) :: uold
  real(kind=prec_real)   :: t=0.
  integer(kind=prec_int) :: nstep=0
end module hydro_commons

module hydro_parameters
  use hydro_precision
  integer(kind=prec_int) :: nx=2
  integer(kind=prec_int) :: ny=2
  integer(kind=prec_int) :: idimbloc=1
  integer(kind=prec_int) :: jdimbloc=1
  integer(kind=prec_int) :: nvar=4
  real(kind=prec_real)   :: dx=1.0
  real(kind=prec_real)   :: tend=0.0
  real(kind=prec_real)   :: gamma=1.4d0
  real(kind=prec_real)   :: courant_factor=0.5d0
  real(kind=prec_real)   :: smallc=1.d-10
  real(kind=prec_real)   :: smallr=1.d-10
  integer(kind=prec_int) :: niter_riemann=10
  integer(kind=prec_int) :: iorder=2
  real(kind=prec_real)   :: slope_type=1.
  character(LEN=20)      :: scheme='muscl'
  integer(kind=prec_int) :: boundary_right=1
  integer(kind=prec_int) :: boundary_left =1
  integer(kind=prec_int) :: boundary_down =1
  integer(kind=prec_int) :: boundary_up   =1
  integer(kind=prec_int) :: noutput=100
  integer(kind=prec_int) :: nstepmax=1000000
  logical                :: on_output=.true.
end module hydro_parameters

module hydro_const
  use hydro_precision
  real(kind=prec_real)   :: zero = 0.0
  real(kind=prec_real)   :: one = 1.0
  real(kind=prec_real)   :: two = 2.0
  real(kind=prec_real)   :: three = 3.0
  real(kind=prec_real)   :: four = 4.0
  real(kind=prec_real)   :: two3rd = 0.6666666666666667d0
  real(kind=prec_real)   :: half = 0.5
  real(kind=prec_real)   :: third = 0.33333333333333333d0
  real(kind=prec_real)   :: forth = 0.25
  real(kind=prec_real)   :: sixth = 0.16666666666666667d0
  integer(kind=prec_int) :: ID=1
  integer(kind=prec_int) :: IU=2
  integer(kind=prec_int) :: IV=3
  integer(kind=prec_int) :: IP=4
end module hydro_const

module hydro_work_space
  use hydro_precision

  ! Work arrays
  real(kind=prec_real),dimension(:,:),allocatable :: u,q,qxm,qxp,dq,qleft,qright,qgdnv,flux
  real(kind=prec_real),dimension(:)  ,allocatable :: c
  real(kind=prec_real),dimension(:)  ,allocatable :: rl,ul,pl,cl,wl
  real(kind=prec_real),dimension(:)  ,allocatable :: rr,ur,pr,cr,wr
  real(kind=prec_real),dimension(:)  ,allocatable :: ro,uo,po,co,wo
  real(kind=prec_real),dimension(:)  ,allocatable :: rstar,ustar,pstar,cstar
  real(kind=prec_real),dimension(:)  ,allocatable :: sgnm,spin,spout,ushock
  real(kind=prec_real),dimension(:)  ,allocatable :: frac,scr,delp,pold
  integer(kind=prec_int),dimension(:),allocatable :: ind,ind2
end module hydro_work_space
