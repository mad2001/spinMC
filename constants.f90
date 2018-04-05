module constants

  implicit none

  private

  integer, parameter, public          :: dp = selected_real_kind(15,200)
  real(kind=dp), parameter, public    :: pi=3.141592653589793238462643383279_dp
  real(kind=dp), parameter, public    :: twopi = 2*pi
  complex(kind=dp), parameter, public :: cmplx_i = (0.0_dp,1.0_dp)
  complex(kind=dp), parameter, public :: cmplx_0 = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter, public :: cmplx_1 = (1.0_dp,0.0_dp)

  real(kind=dp), parameter, public :: bohr = 0.5291772108_dp

  real(kind=dp), parameter, public    :: eps2  = 1.0e-2_dp
  real(kind=dp), parameter, public    :: eps5  = 1.0e-5_dp
  real(kind=dp), parameter, public    :: eps6  = 1.0e-6_dp
  real(kind=dp), parameter, public    :: eps7  = 1.0e-7_dp
  real(kind=dp), parameter, public    :: eps8  = 1.0e-8_dp
  real(kind=dp), parameter, public    :: eps10 = 1.0e-10_dp

  integer, parameter, public          :: io_unit = 101
  integer, parameter, public          :: in_unit = 201
  integer, parameter, public          :: states_unit = 401
  integer, parameter, public          :: en_unit = 501
  integer, parameter, public          :: chg_unit = 601
  integer, parameter, public          :: susz_unit = 701
  integer, parameter, public          :: pol_unit = 801
  integer, parameter, public          :: conv_unit = 1001
end module constants
