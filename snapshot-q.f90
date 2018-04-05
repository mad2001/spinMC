program snapshotq

use constants
use utility
use parameters
use sstates

implicit none
real(kind=dp),allocatable  :: q_tchg(:,:)

character(len=50)          :: in_name, out_name

integer           :: narg

complex(kind=dp)  :: exp1
integer           :: ix,iy
real(kind=dp)     :: state0(3),state1(3),state2(3),state3(3)


narg=iargc()
if (narg<2) then
  write(*,*) 'usage: spinmc-snapshot.x in_name out_name'
  stop
endif

call getarg(1, in_name)
call getarg(2, out_name)


open(in_unit,file='INPUT',form='formatted',status='old',action='read')
call ini_param(in_unit) !read parameters
close(in_unit)

call ini_sstates() !initialize spin states
call ini_nns()

open(io_unit,file=trim(in_name),form='formatted',status='old',action='read')
call read_states(io_unit)
close(io_unit)

 
allocate(q_tchg(n_l,n_l))
q_tchg(:,:)=0.0

do iy=1,n_l
do ix=1,n_l
  state0 = s_states(:,ix,iy)
  state1 = s_states(:,n_nns(1,1,ix,iy),n_nns(2,1,ix,iy))
  state2 = s_states(:,n_nns(1,3,ix,iy),n_nns(2,3,ix,iy))

  call cross(state1,state2,state3)
  exp1=cmplx_1+dot_product(state0,state1)+&
     &dot_product(state1,state2)+&
     &dot_product(state2,state0)+&
     &cmplx_i*dot_product(state0,state3)
  q_tchg(ix,iy) = q_tchg(ix,iy) + aimag(log(exp1))

  state0 = s_states(:,ix,iy)
  state1 = s_states(:,n_nns(1,2,ix,iy),n_nns(2,2,ix,iy))
  state2 = s_states(:,n_nns(1,4,ix,iy),n_nns(2,4,ix,iy))
  call cross(state1,state2,state3)
  exp1=cmplx_1+dot_product(state0,state1)+&
     &dot_product(state1,state2)+&
     &dot_product(state2,state0)+&
     &cmplx_i*dot_product(state0,state3)
  q_tchg(ix,iy) = q_tchg(ix,iy) + aimag(log(exp1))
end do
end do

open(io_unit,file=trim(out_name),form='formatted',status='replace')

do iy=1,n_l
do ix=1,n_l
  write(io_unit,'(1x,2I8,E16.8)') ix, iy, q_tchg(ix,iy)
end do
end do

close(io_unit)

deallocate(q_tchg)

end program
