module lattice
use constants,  only : dp
use parameters, only : n_l
implicit none

integer, parameter     :: nnn = 4 !number of nearest neighhor
integer, parameter     :: nnn2 = 2 !half of nnn

contains
!========
!output spin states
!========
subroutine out_states_l(iounit,s_states)
  integer,      intent(in)    :: iounit
  real(kind=dp),intent(inout) :: s_states(:,:)
  integer                     :: ix,iy,ii,it
  do iy=1,n_l
  do ix=1,n_l
    it = get_index(ix,iy)
    write(iounit,'(1x,2I8,3f20.12)') ix,iy,(s_states(ii,it),ii=1,3)
  end do
  end do
end subroutine out_states_l
!========
!read spin states
!========
subroutine read_states_l(iounit,s_states)
  integer,      intent(in)    :: iounit
  real(kind=dp),intent(inout) :: s_states(:,:)
  integer                     :: ix,iy,ii,it,n_s
  n_s = get_n_s()
  do it=1,n_s
    read(iounit,*) ix, iy, (s_states(ii,it),ii=1,3)
  end do
end subroutine read_states_l
!========
!initialize the neighbor and the direction of DM
!right left up down
!========
subroutine ini_nns_l(n_nns, vec_nD)
  integer                      :: ix,iy,it,tx,ty
  integer,       intent(inout) :: n_nns(:,:)
  real(kind=dp), intent(inout) :: vec_nD(:,:)
  n_nns = 0
  do iy=1,n_l
  do ix=1,n_l
    it = get_index(ix,iy)
! right
    tx = ix+1
    ty = iy
    if (tx > n_l) tx = 1
    n_nns(1,it) = get_index(tx,ty)
! left
    tx = ix-1
    ty = iy
    if (tx < 1) tx = n_l
    n_nns(2,it) = get_index(tx,ty)
! up
    tx = ix
    ty = iy+1
    if (ty > n_l) ty = 1
    n_nns(3,it) = get_index(tx,ty)
! down
    tx = ix
    ty = iy-1
    if (ty < 1) ty = n_l
    n_nns(4,it) = get_index(tx,ty)
  end do
  end do
  vec_nD(:,1)=(/1.0_dp,0.0_dp,0.0_dp/)
  vec_nD(:,2)=(/0.0_dp,1.0_dp,0.0_dp/)
end subroutine ini_nns_l

integer function get_index(ix,iy)
  integer, intent(in) :: ix,iy
  get_index = ix+(iy-1)*n_l
end function get_index

integer function get_n_s()
  get_n_s = n_l*n_l
end function get_n_s

integer function get_n_states()
  get_n_states = n_l*n_l
end function get_n_states

end module lattice
