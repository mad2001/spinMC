program spinmc_mpi
include 'mpif.h'

!implicit none
character(len=120)     :: cmd
integer                :: narg, n_proc, ii, ierr
character(len=20)      :: jobx,dirfile, buff,opt
!character, allocatable :: dirlists(:,:)

narg=iargc()
if (narg<2) then
  write(*,'(a)') 'usage: spinmc_mpi.x <n_proc> <serial job x> <dirsuffix> <option>'
  stop
else
  call getarg(1, buff)
  read(buff,*) n_proc
  call getarg(2, jobx)
  call getarg(3, dirfile)
  call getarg(4, opt)
!  call getarg(2, dirfile)
!  open(105,file=trim(dirfile),form='formatted',status='old',action='read')
!  do ii=1,n_proc
!    read(105,'(a)') dirlists(:,ii)
!    write(*,'(I4,a)') ii, dirlists(:,ii)
!  end do
!  close(105)
end if

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world,ii,ierr)
write(cmd,'(a,1x,a,I2.2,1x,a)'),trim(jobx), trim(dirfile), ii+1, trim(opt)
call system(trim(cmd))
call mpi_finalize(ierr)

end program
