module utility
use constants, only : dp
  implicit none

contains

subroutine cross(v1,v2,vr)
  real(kind=dp),dimension(3)::v1,v2,vr
  vr(1) = v1(2) * v2(3) - v1(3) * v2(2)
  vr(2) = v1(3) * v2(1) - v1(1) * v2(3)
  vr(3) = v1(1) * v2(2) - v1(2) * v2(1)
end subroutine cross

subroutine vecNorm(vector)
  real(kind=dp),dimension(3),intent(inout) :: vector
  real(kind=dp)                            :: v0
  v0=dot_product(vector,vector)
  vector(:)=vector(:)/sqrt(v0)
  return
end subroutine vecNorm

!function get_rnds(spin)
!  real(kind=dp),intent(out) :: get_rnds(3)
!  real(kind=dp),intent(in)  :: spin(3)
!  real(kind=dp)             :: ri,rk,rx,ry,r1,r2,r3
!  ri=2.0_dp
!  do while (ri>=1.0_dp)
!    call random_number(rnd_tmp)
!    rx=1.0_dp-2.0_dp*rnd_tmp
!    call random_number(rnd_tmp)
!    ry=1.0_dp-2.0_dp*rnd_tmp
!    ri=rx*rx+ry*ry
!  end do
!  rk=sqrt(1.0_dp-ri)
!  r1=2.0_dp*rx*rk
!  r2=2.0_dp*ry*rk
!  r3=1.0_dp-2.0_dp*ri
!
!  get_rnds = 0.9682458365518542213_dp*spin+0.25_dp*(/r1,r2,r3/)
!  return
!end function get_rnds

subroutine init_random_seed()

    use iso_fortran_env, only: int64

    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                     + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                     + dt(3) * 24_int64 * 60 * 60 * 1000 &
                     + dt(5) * 60 * 60 * 1000 &
                     + dt(6) * 60 * 1000 + dt(7) * 1000 &
                     + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
end subroutine init_random_seed

end module utility
