program main_NR

  USE OMP_LIB
  use Structure
  use Matrix_op
  use gnuplot_fortran
  use bar
  use NR_module
  use local_stage_NR


  implicit none

  real :: start, finish

  ! Solution fields:
  double precision, dimension(:,:), allocatable :: U ! Displacement.
  double precision, dimension(:,:), allocatable :: eps ! Strain.
  double precision, dimension(:,:), allocatable :: sigma ! Stress.

  ! Other fields:
  integer :: nelm
  double precision, dimension(:), allocatable :: plot_vec
  double precision :: x
  double precision, dimension(3, 3) :: m,minv
  double precision, dimension(5) :: v ! vector
  !integer, dimension(10) :: inc ! vector
  integer, dimension(:,:), allocatable :: inc ! vector
  integer :: j,max_iter=15


  ! Latin definitions:
  double precision :: LATIN_error
  integer :: nL = 0 ! Latin iteration

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call cpu_time(start)

  print*, (' ')
  print*, (' ')
  print*, (' ')
  print*, (' ')
  print*, (' ######     ######    ####    ####    ll                    ')
  print*, (' ##    ##  ##    ##   #####  #####    ll                    ')
  print*, (' ##    ##  ##    ##   ##   ##   ##    ll            bb      ')
  print*, (' ########  ##    ##   ##        ##    ll   aaaa     bb      ')
  print*, (' ##    ##  ##    ##   ##        ##    ll  aa  aa    bbbbb   ')
  print*, (' ##    ##  ##    ##   ##        ##    ll  aa  aaa   bb  bb  ')
  print*, (' ##    ##   ######    ##        ##    ll   aaaa aa   bbbb   ')
  print*, (' ')



  print*, "NEWTON RAPHSON SOLVER"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  call Input_assign  ! Constructor of the structure "I"
  call Time_op_assign(I) ! Constructor of the structure "Top"
  call Space_elm_matrices(I)  ! Creation of structure Sop (containing the space elemental matrices)


  ! NR solver:
  call NR_solver(Top,Sop,L,I)




  print*, " "
  print*, "PROBLEM CONVERGED."



  allocate( plot_vec( I%Tngpt )  ) ! vector used to plot.
  nelm = I%Nx ! Element chose to plot.
  plot_vec = reshape( L%d_Ep(1,1,nelm,:) , (/I%Tngpt/) )  ! Plastic deformation
  !plot_vec = reshape( Elast%sigma(1,1,nelm,:) , (/I%Tngpt/) )   ! stress Elastic solution
  call plot(plot_vec)
  call cpu_time(finish)
  print '("RESOLUTION TIME = ",f6.3," seconds.")', finish-start





end program main_NR
