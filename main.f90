program spinMC

use constants
use utility
use parameters
use sstates

implicit none
real(kind=dp)     :: temperature, old_temperature
real(kind=dp)     :: energy, d_energy, dsq_energy, sq_energy, cv ! energy
real(kind=dp)     :: magz, d_magz, dsq_magz, sq_magz, susz !susceptibility
real(kind=dp)     :: tcharge, d_tcharge, sq_tcharge, err_tcharge !topological charge
real(kind=dp)     :: orderparam !correlation
real(kind=dp)     :: polar(3),d_polar(3)
integer           :: itt,jtt,ktt,vtt,narg
real(kind=dp)     :: ctime
character(len=20) :: dirname
character(len=50) :: sfilename
character(len=2)  :: opt
!real(kind=dp),dimension(3)  :: s00, s01, s02

narg=iargc()
if (narg<1) then
  write(*,'(a)') 'usage: spinmc.x dirname [option]'
  write(*,'(a)') 'option:'
  write(*,'(a)') '-i : default, both thermal and average, output states before average'
  write(*,'(a)') '-t : only thermal, output states at the end'
  write(*,'(a)') '-a : only average, reading start_states.txt before average'
  stop
endif
call getarg(1, dirname)
open(in_unit,file=trim(dirname)//"/INPUT",form='formatted',status='old',action='read')
call ini_param(in_unit) !read parameters
close(in_unit)

if (narg==2) then
  call getarg(2, opt)
else
  opt = '-i'
end if


open(unit=io_unit,file=trim(dirname)//'/log.txt',form='formatted',status='replace')
call out_param_log(io_unit)
close(io_unit)

call ini_sstates() !initialize spin states

call ini_nns() !initialize neighor maps

open(unit=io_unit,file=trim(dirname)//'/nns.txt',form='formatted',status='replace')
call out_nns(io_unit)
close(io_unit)

if((opt=='-i').or.(opt=='-t')) then
!=== debug ===
!write(*,*) 'begin initial states'
!do jtt = 1, n_l
!  do ktt = 1, n_l
!write(*,'(1x,2I4,3F14.8)') ktt, jtt, (s_states(vtt,ktt,jtt), vtt=1,3)
!  end do
!end do
!write(*,*) 'end initial states'
!=============

temperature = T_h
itt = 0
call init_random_seed() !this is important for gfortran
do while(temperature>=T_u)
!  write(*,*) 'J D B is', r_J, r_D, r_B !debug
  itt = itt+1
!  if(mod(itt,20)==0) call init_random_seed() !this is important for gfortran
  orderparam = 0.0_dp
  call update_beta(temperature)

  do jtt=1, n_throw
    call update_onestep()
  end do !jtt
 
  do jtt=1, n_thermal
    call update_onestep()

    orderparam = orderparam + correlation()
  end do !jtt
  orderparam = orderparam/n_thermal

  call CPU_TIME(ctime)
  open(unit=io_unit,file=trim(dirname)//'/log.txt',form='formatted',status='old',&
      &position='append',action='write')
  write(io_unit,'(1x,I8,f12.6,E16.8,f16.2,a)') itt, temperature, orderparam, ctime, ' THML'
  close(io_unit)
  temperature = temperature*d_thermal
end do

call CPU_TIME(ctime)
open(unit=io_unit,file=trim(dirname)//'/log.txt',form='formatted',status='old',&
    &position='append',action='write')
write(io_unit,'(1x,a,f16.2)') 'time for thermalization:', ctime
close(io_unit)

!========================
  open(unit=states_unit,file=trim(dirname)//'/start_states.txt',&
      &form='formatted',status='replace')
  call out_states(states_unit) 
  close(states_unit)
!========================
  if(opt=='-t') stop !only thermal
else
  if(opt=='-a') then !only average,read start_states.txt

  open(unit=states_unit,file=trim(dirname)//'/start_states.txt',&
      &form='formatted',status='old',action='read')
  call read_states(states_unit)
  close(states_unit)
  else !wrong option arg

    write(*,'(a)') '-i : default, both thermal and average, output states before average'
    write(*,'(a)') '-t : only thermal, output states at the end'
    write(*,'(a)') '-a : only average, reading start_states.txt before average'
    stop
  end if
end if !opt

if (l_energy) then
  open(unit=en_unit,file=trim(dirname)//'/Energy.txt',&
     &form='formatted',status='replace')
  close(en_unit)
end if

if (l_tchg) then
  open(unit=chg_unit,file=trim(dirname)//'/Tcharge.txt',&
     &form='formatted',status='replace')
  close(chg_unit)
end if

if (l_susz) then
  open(unit=susz_unit,file=trim(dirname)//'/Susceptibility_z.txt',&
     &form='formatted',status='replace')
  close(susz_unit)
end if

if (l_pol) then
  open(unit=pol_unit,file=trim(dirname)//'/Polarization.txt',&
     &form='formatted',status='replace')
  close(pol_unit)
end if

temperature = T_u
old_temperature = T_2
itt = 0
call init_random_seed() !this is important for gfortran
do while(temperature>=T_0)
  itt = itt+1
  tcharge = 0.0_dp
  energy = 0.0_dp
  orderparam = 0.0_dp
  magz = 0.0_dp
  polar = 0.0_dp
  sq_energy = 0.0_dp
  sq_tcharge = 0.0_dp
  sq_magz = 0.0_dp
  call update_beta(temperature)

  do jtt=1, n_throw
    call update_onestep()
  end do !jtt

!  jtt = 0
  do jtt=1, n_avrg
    do ktt=1,n_mcs
      call update_onestep()
    end do

    if (l_energy) then
      d_energy = fullEnergy()
      energy = energy + d_energy
      !sq_energy = sq_energy + dsq_energy
      sq_energy = sq_energy + d_energy*d_energy
    end if

    if (l_tchg) then
      d_tcharge = topologicalCharge()
      tcharge = tcharge + d_tcharge
      sq_tcharge = sq_tcharge + d_tcharge*d_tcharge
    end if
    
    if (l_susz) then
      d_magz = magnetisation()
      magz = magz + d_magz
      !sq_magz = sq_magz + dsq_magz    
      sq_magz = sq_magz + d_magz*d_magz
    end if

    if (l_pol) then
      call polarization(d_polar)
      polar = polar + d_polar
    end if

    orderparam = orderparam + correlation()
  end do !avrg

  if (l_energy) then
    energy = energy/n_avrg
    sq_energy = sq_energy/n_avrg
    cv = (sq_energy-energy*energy)*beta*beta
    open(unit=en_unit,file=trim(dirname)//'/Energy.txt',&
       &form='formatted',status='old',position='append',action='write')
    write(en_unit,'(1x,f12.6,2E16.8)') temperature, energy, cv
    close(en_unit)
  end if  

  if (l_tchg) then 
    tcharge = tcharge/n_avrg
    sq_tcharge = sq_tcharge/n_avrg
    err_tcharge = sqrt(sq_tcharge-tcharge*tcharge)
    open(unit=chg_unit,file=trim(dirname)//'/Tcharge.txt',&
       &form='formatted',status='old',position='append',action='write')
    write(chg_unit,'(1x,f12.6,2E16.8)') temperature, tcharge, err_tcharge
    close(chg_unit)
  end if

  if (l_susz) then
    magz = magz/n_avrg
    sq_magz = sq_magz/n_avrg
    susz = (sq_magz-magz*magz)*beta
    open(unit=susz_unit,file=trim(dirname)//'/Susceptibility_z.txt',&
       &form='formatted',status='old',position='append',action='write')
    write(susz_unit,'(1x,f12.6,2E16.8)') temperature, magz, susz
    close(susz_unit)
  end if
  
  if (l_pol) then
    polar = polar/n_avrg
    open(unit=pol_unit,file=trim(dirname)//'/Polarization.txt',&
       &form='formatted',status='old',position='append',action='write')
    write(pol_unit,'(1x,f12.6,3E16.8)') temperature, (polar(vtt),vtt=1,3)
    close(pol_unit)
  end if    

  orderparam = orderparam/n_avrg
  call CPU_TIME(ctime)
  open(unit=io_unit,file=trim(dirname)//'/log.txt',form='formatted',status='old',&
     &position='append',action='write')
  write(io_unit,'(1x,I8,f12.6,E16.8,f16.2,a)') &
     &itt, temperature, orderparam, ctime, '  AVRG'
  close(io_unit)
!======== annealing and snapshot ======== 
  if (temperature > t_2) then
    temperature = temperature*d_avrg2    
  elseif (temperature > t_1) then
    if(old_temperature - temperature > d_snap1) then
      write(sfilename, "(a,'/states_T',f8.6,'.txt')") trim(dirname), temperature
      open(unit=states_unit,file=sfilename,&
         &form='formatted',status='replace')
      call out_states(states_unit)
      close(states_unit)
      old_temperature = temperature
    endif
    temperature = temperature*d_avrg1
  else 
    if(old_temperature - temperature > d_snap0) then
      write(sfilename, "(a,'/states_T',f8.6,'.txt')") trim(dirname), temperature
      open(unit=states_unit,file=sfilename,&
         &form='formatted',status='replace')
      call out_states(states_unit)
      close(states_unit)
      old_temperature = temperature
    endif
    temperature = temperature*d_avrg0
  end if
!========
end do !temperature

open(unit=states_unit,file=trim(dirname)//'/end_states.txt',&
     &form='formatted',status='replace')
call out_states(states_unit)
close(io_unit)

call deallocate_all()

call CPU_TIME(ctime)
open(unit=io_unit,file=trim(dirname)//'/log.txt',form='formatted',status='old',&
    &position='append',action='write')
write(io_unit,'(1x,a,f16.2)') 'time for average:', ctime
close(io_unit)

end program
