program main

  USE OMP_LIB
  use Structure
  use Matrix_op
  use gnuplot_fortran
  use bar
  use Time_matrices
  !use Elastic_sol
  use Elastic_sol_PGD
  use global_stage
  use local_stage
  !use local_stage_vec
  use LATIN_indicator

  !include 'solver.f03' !use funcion

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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  call Input_assign  ! Constructor of the structure "I"
  call Time_op_assign(I) ! Constructor of the structure "Top"
  call Space_elm_matrices(I)  ! Creation of structure Sop (containing the space elemental matrices)

  !!$OMP PARALLEL

  !call solve_init(I,Sop)  ! Creation of Elast and solve of the initial elastic problem.
  call solve_init_PGD(I,Sop,Top)  ! Creation of Elast and solve of the initial elastic problem.

  ! Initialization of the LATIN iterations:
  LATIN_error = 100

  call assign_elastic(Elast,I) ! Assignation of the elastic solution to the global_stage quantities.

  !call local_stage_vp(I,Sop,G%sigma)
  allocate( plot_vec( I%Tngpt )  ) ! vector used to plot.



  do while (I%LATIN_error_min <= LATIN_error )


      nL = nL + 1

      ! Local stage:
      call local_stage_vp(I,Sop,G%sigma)
      !call local_stage_vp_vec(I,Sop,G%sigma)


      !! TEST PLOTS

        !print*, L%d_Ep
        nelm = 100 ! Element chose to plot.
        plot_vec = reshape( L%d_Ep(1,1,nelm,:) , (/I%Tngpt/) )  ! Plastic deformation
        !plot_vec = reshape( Elast%sigma(1,1,nelm,:) , (/I%Tngpt/) )   ! stress Elastic solution
        !call plot(plot_vec)
        !stop
        !print*, (maxval(plot_vec,1))
        !print*, I%ks

     !! TEST PLOTS

      ! Global stage:
      call global_stage_vp(Elast,Top,Sop,L,I,nL)


      if (I%PGD_eval) then !(.true.) then
        ! Calculate the error of the solver:
        LATIN_error= indicator1(G,L)
        print*, "LATIN ERROR [%]:" , LATIN_error, "ITERATION:", nL
      end if

  end do

  print*, " "
  print*, "PROBLEM CONVERGED."

  nelm = I%Nx ! Element chose to plot.
  plot_vec = reshape( L%d_Ep(1,1,nelm,:) , (/I%Tngpt/) )  ! Plastic deformation
  !plot_vec = reshape( Elast%sigma(1,1,nelm,:) , (/I%Tngpt/) )   ! stress Elastic solution
  call plot(plot_vec)


  call cpu_time(finish)
  print '("RESOLUTION TIME = ",f6.3," seconds.")', finish-start

  !!$OMP END PARALLEL

  !v = linspace(0.D0,10.D0,5)





end program main
