module sstates

use constants,  only : dp
use lattice,    only : nnn, nnn2
implicit none

!==== spin states
real(kind=dp), allocatable :: s_states(:,:) !spin matrix (3,n_s)
integer,       allocatable :: n_nns(:,:) !nearest neighor site index (6,n_s)
real(kind=dp), allocatable :: vec_nD(:,:)
integer,       save        :: n_s

real(kind=dp), save        :: beta !* 1/kT
! private variables
integer,       private     :: ierr
real(kind=dp), private     :: rnd_tmp
contains

!========
!output spin states
!========
subroutine out_states(iounit)
  use lattice, only : out_states_l
  integer :: iounit
  call out_states_l(iounit,s_states)
end subroutine out_states
!========
!read spin states
!========
subroutine read_states(iounit)
  use lattice, only : read_states_l
  integer :: iounit
  call read_states_l(iounit,s_states)
end subroutine read_states

!========
!initialize the neighbor and the direction of DM
!========
subroutine ini_nns()
  use lattice, only : ini_nns_l
  integer          :: it,inn
  allocate(n_nns(nnn,n_s),stat=ierr)
  allocate(vec_nD(3,nnn2),stat=ierr)
  call ini_nns_l(n_nns,vec_nD)
end subroutine ini_nns

subroutine out_nns(iounit)
  integer,intent(in):: iounit
  integer           :: it,inn
  do it=1,n_s
    write(iounit,*) (n_nns(inn,it),inn=1,nnn)
  end do
  do it=1,nnn2
    write(iounit,'(3f16.8)') (vec_nD(inn,it),inn=1,3)
  end do
end subroutine out_nns
!========
!initialize the spin states
!========
subroutine ini_sstates()
  use lattice, only : get_n_s
  integer          :: it
  real(kind=dp)    :: ri,rx,ry,rk
  n_s = get_n_s()
  allocate(s_states(3,n_s+1),stat=ierr)
  s_states = 0.0_dp
  do it=1,n_s
    ri=2.0_dp
    do while (ri>=1.0_dp)
      call random_number(rnd_tmp)
      rx=1.0_dp-2.0_dp*rnd_tmp
      call random_number(rnd_tmp)
      ry=1.0_dp-2.0_dp*rnd_tmp
      ri=rx*rx+ry*ry
    end do
    rk=sqrt(1.0_dp-ri)
    s_states(1,it)=2.0_dp*rx*rk
    s_states(2,it)=2.0_dp*ry*rk
    s_states(3,it)=1.0_dp-2.0_dp*ri
  end do
end subroutine ini_sstates

!========
!transfer temperature to beta
!========
subroutine update_beta(temperature)
  real(kind=dp), intent(in) :: temperature
  beta = 1.0_dp/temperature
end subroutine

subroutine deallocate_all()
  deallocate(s_states,n_nns)
end subroutine deallocate_all

!========
!update one spin
!========
subroutine rndlat(sr,s0)
  use utility, only : vecNorm
  use parameters, only : T_0
  real(kind=dp),intent(out) :: sr(3)
  real(kind=dp),intent(in)  :: s0(3)
  real(kind=dp)             :: ri,rx,ry,rk
  real(kind=dp)             :: state0(3)
  ri=2.0_dp
  do while (ri>=1.0_dp)
    call random_number(rnd_tmp)
    rx=1.0_dp-2.0_dp*rnd_tmp
    call random_number(rnd_tmp)
    ry=1.0_dp-2.0_dp*rnd_tmp
    ri=rx*rx+ry*ry
  end do
  rk=sqrt(1.0_dp-ri)
  state0(1)=2.0_dp*rx*rk
  state0(2)=2.0_dp*ry*rk
  state0(3)=1.0_dp-2.0_dp*ri

  sr = (beta*T_0)*s0+0.25_dp*state0
  call vecNorm(sr)
end subroutine rndlat

!========
! effective field on one site
!========
subroutine field(eff,it,is_int)
  use utility,    only : cross
  use parameters, only : r_J,r_D,r_B,r_K,n_l
  real(kind=dp),intent(out) :: eff(3)
  integer, intent(in)       :: it
!====
!is intergral(total energy) or not
!Magnetic field is special
!when we calculate single particle field, it's B
!when we calculate total energy, it's 2B because we need half later.
!====
  logical                   :: is_int
  integer                   :: ii
  real(kind=dp)             :: tmpf(3),tmpcross(3),state0(3)
  eff=0.0_dp
  do ii = 1, nnn !neighbor
    eff = eff + s_states(:,n_nns(ii,it))
  end do
  eff = r_J*eff

  tmpf=0.0_dp
  do ii = 1, nnn2 !half of neighbor
    state0 = s_states(:,n_nns(ii*2-1,it))+s_states(:,n_nns(ii*2,it))
    call cross(vec_nD(:,ii),state0,tmpcross)
    if (mod(it+it/n_l, 2) == 0) then
      tmpf = tmpf - tmpcross
    else
      tmpf = tmpf + tmpcross
    endif
  end do
  eff = eff + r_D*tmpf

  eff(3) = eff(3) + r_K*s_states(3,it) !anisotropic potential
  if(is_int) then
    eff(3) = eff(3) + 2.0_dp*r_B
  else
    eff(3) = eff(3) + r_B
  endif
end subroutine field

!========
!Metroplis MC update
!========
subroutine mcupdate()
  integer                :: it
  real(kind=dp)          :: diffEn
  real(kind=dp)          :: state1(3),state2(3),effField(3)
  do it = 1,n_s
    call rndlat(state1, s_states(:,it))
    call field(effField,it,.false.)
    state2 = s_states(:,it) - state1
    diffEn = dot_product(state2,effField)
    if (diffEn < 0.0_dp) then
      s_states(:,it)=state1
    else
      call random_number(rnd_tmp)
      if (rnd_tmp<exp(-diffEn*beta)) s_states(:,it)=state1
    end if
  end do
end subroutine mcupdate

!========
!which only can be 0 or 1
!========
subroutine overrelaxation()
  integer                    :: it
  real(kind=dp)              :: effField(3),state0(3),state1(3)
!  integer, intent(in) :: which
  do it = 1,n_s
!    if (mod((ix+iy),2) == which) cycle
    state0 = s_states(:,it)
    call field(effField,it,.false.)
    state1=2.0_dp*effField(:)*dot_product(state0,effField)/dot_product(effField,effField)-state0
    s_states(:,it) = state1
  end do
end subroutine overrelaxation

subroutine statesNorm()
  use utility, only:vecNorm
  integer  :: it
  do it = 1,n_s
    call vecNorm(s_states(:,it))
  end do
end subroutine statesNorm

subroutine update_onestep()
  use parameters, only : n_over
  integer        :: ii
  do ii = 1, n_over
    call overrelaxation()
  end do
  call mcupdate()
end subroutine update_onestep

!========
!total energy
!========
real(kind=dp) function fullEnergy()
  real(kind=dp)              :: e
  real(kind=dp)              :: de
  integer                    :: it
  real(kind=dp)              :: effField(3)
  e=0.0_dp
  do it=1,n_s
    call field(effField,it,.true.)
    de=dot_product(s_states(:,it),effField)
    e=e-de
  end do
  fullEnergy = e/(2.0_dp*n_s)
end function fullEnergy

!========
!correlation
!========
real(kind=dp) function correlation()
  real(kind=dp)              :: c
  integer                    :: it, ii
  real(kind=dp)              :: state0(3)
  c=0.0_dp
  do it=1,n_s
    state0=0.0_dp ! sum of neighor spins
    do ii=1,nnn !neighbor
      state0 = state0 + s_states(:,n_nns(ii,it))
    end do
    c = c + dot_product(s_states(:,it),state0)
  end do
  correlation = c/(nnn*n_s)
end function correlation

!========
!topological charge
!========
real(kind=dp) function topologicalCharge()
  use constants, only : cmplx_i, cmplx_0, cmplx_1, twopi
  use utility,   only : cross
  complex(kind=dp)  :: exp1
  real(kind=dp)     :: q
  integer           :: it
  real(kind=dp)     :: state0(3),state1(3),state2(3),state3(3)
  q = 0.0_dp
  do it=1,n_s
      state0 = s_states(:,it)
      state1 = s_states(:,n_nns(1,it))
      state2 = s_states(:,n_nns(3,it))

      call cross(state1,state2,state3)

      exp1=cmplx_1+dot_product(state0,state1)+&
         &dot_product(state1,state2)+&
         &dot_product(state2,state0)+&
         &cmplx_i*dot_product(state0,state3)

      q = q + aimag(log(exp1))

      state0 = s_states(:,it)
      state1 = s_states(:,n_nns(2,it))
      state2 = s_states(:,n_nns(4,it))

      call cross(state1,state2,state3)

      exp1=cmplx_1+dot_product(state0,state1)+&
         &dot_product(state1,state2)+&
         &dot_product(state2,state0)+&
         &cmplx_i*dot_product(state0,state3)

      q = q + aimag(log(exp1))
  end do
  topologicalCharge = q/(twopi*n_s)
end function topologicalCharge

real(kind=dp) function magnetisation()
  real(kind=dp)              :: m
  integer                    :: it
  m = 0.0_dp
  do it=1,n_s
    m = m + s_states(3,it)
  end do
  magnetisation = m/n_s
end function magnetisation

subroutine polarization(polar)
  real(kind=dp)             :: polar(3)
  integer                   :: it
  polar = 0.0_dp
  do it=1,n_s
    polar(1)=polar(1)+s_states(2,it)*s_states(3,it)
    polar(2)=polar(2)+s_states(3,it)*s_states(1,it)
    polar(3)=polar(3)+s_states(1,it)*s_states(2,it)
  end do
  polar=polar/n_s
end subroutine polarization

end module sstates
