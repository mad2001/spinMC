module parameters

use constants,  only : dp

implicit none

!==== control parameters
integer, save       :: n_thermal! = 20000 ! thermal process number for each temperature
integer, save       :: n_avrg ! = 10000 ! update number for each site
integer, save       :: n_throw !number of update will not be used
integer, save       :: n_over ! overrelaxation number per mcupdate
integer, save       :: n_mcs ! = 20 !

!====

!==== interaction parameters
integer, save       :: n_l
real(kind=dp), save :: r_J = 1.00_dp !*Jij interaction, J>0 FM, J<0 AF
real(kind=dp), save :: r_D ! DM interaction
real(kind=dp), save :: r_B ! magnetic field, always >0
real(kind=dp), save :: r_K ! anisotropic field, K>0 zaxis, K<0 xyaxis

!==== temperature parameters
real(kind=dp), save :: d_thermal ! = 0.96
real(kind=dp), save :: d_avrg2 ! = 0.96
real(kind=dp), save :: d_avrg1 ! = 0.92
real(kind=dp), save :: d_avrg0 ! = 0.88

real(kind=dp), save :: T_h ! 10.00
real(kind=dp), save :: T_u !  4.00
real(kind=dp), save :: T_2 !  1.20 
real(kind=dp), save :: T_1 !  0.40
real(kind=dp), save :: T_0 !  0.01


!==== control logical
logical,save        :: l_energy
logical,save        :: l_tchg
logical,save        :: l_susz
logical,save        :: l_pol
!==== control snapshot
real(kind=dp), save :: d_snap1
real(kind=dp), save :: d_snap0

contains

!========
!read parameters from INPUT file
!========
subroutine ini_param(iounit)
  integer :: iounit

  read(iounit,*) n_l
  read(iounit,*) r_J, r_D, r_K, r_B
  read(iounit,*) n_thermal, n_avrg
  read(iounit,*) n_throw, n_over, n_mcs
  read(iounit,*) d_thermal, d_avrg2, d_avrg1, d_avrg0
  read(iounit,*) T_h, T_u, T_2, T_1, T_0
  read(iounit,*) l_energy, l_tchg, l_susz, l_pol
  read(iounit,*) d_snap1, d_snap0
end subroutine ini_param

!========
!output parameters in the logfile
!========
subroutine out_param_log(iounit)
  integer :: iounit

  write(iounit,'(1x,a,i8)') " Size =", n_l
  write(iounit,'(1x,a,f8.4,a,f8.4,a,f8.4,a,f8.4)')&
       & ' J =', r_J, ' DM =', r_D, ' K =', r_K, ' B =', r_B
  write(iounit,'(1x,a,i8,a,i8)')&
       & ' n_thermal =', n_thermal,' n_avrg =', n_avrg
  write(iounit,'(1x,a,i8,a,i4,a,i4)')&
       & ' n_throw =', n_throw, ' n_over =', n_over, ' n_mcs =', n_mcs
  write(iounit,'(1x,a,f8.4,a,f8.4)')&
       & ' d_thermal =',d_thermal,' d_avrg2 =',d_avrg2,' d_avrg1 =',d_avrg1,' d_avrg0 =',d_avrg0
  write(iounit,'(1x,a,f8.4,a,f8.4,a,f8.4,a,f8.4)')&
       & ' T_3 =', T_h, ' T_2 =', T_u, ' T_2 =', T_2,' T_1 =', T_1, ' T_0 =', T_0
  write(iounit,'(1x,a,l2,a,l2,a,l2,a,l2,a,l2)')&
       & ' l_energy =',l_energy,' l_tchg =',l_tchg,' l_susz =',l_susz,' l_pol =',l_pol
  write(iounit,'(1x,a,f8.4,a,f8.4,a,f8.4,a,f8.4)')&
       & ' d_snap1 =', d_snap1, ' d_snap0 =', d_snap0
end subroutine out_param_log

!========
!output parameters as new INPUT file
!========
subroutine out_param_input(iounit)
  integer :: iounit

  write(iounit,'(i8)') n_l
  write(iounit,'(3f8.4)') r_J, r_D, r_K, r_B
  write(iounit,'(3i8)') n_thermal, n_avrg
  write(iounit,'(3i8)') n_throw, n_over, n_mcs
  write(iounit,'(2f8.4)') d_thermal, d_avrg2, d_avrg1, d_avrg0
  write(iounit,'(3f8.4)') T_h, T_u, T_2, T_1, T_0
  write(iounit,'(4l2)') l_energy, l_tchg, l_susz, l_pol
  write(iounit,'(2f8.4)') d_snap1, d_snap0
end subroutine out_param_input

end module
