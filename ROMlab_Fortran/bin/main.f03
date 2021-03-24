program main

  USE OMP_LIB
  use Structure
  use Matrix_op
  use gnuplot_fortran
  use bar
  use Time_matrices
  use Elastic_sol
  use global_stage
  use local_stage
  use LATIN_indicator

  !include 'solver.f03' !use funcion

  implicit none


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








  call Input_assign  ! Constructor of the structure "I"
  call Time_op_assign(I) ! Constructor of the structure "Top"

  call Space_elm_matrices(I)  ! Creation of structure Sop (containing the space elemental matrices)

  !!$OMP PARALLEL

  call solve_init(I,Sop)  ! Creation of Elast and solve of the initial elastic problem.


  ! Initialization of the LATIN iterations:
  LATIN_error = 100

  call assign_elastic(Elast,I) ! Assignation of the elastic solution to the global_stage quantities.

  !call local_stage_vp(I,Sop,G%sigma)
  allocate( plot_vec( I%Tngpt )  ) ! vector used to plot.



  do while (I%LATIN_error_min <= LATIN_error )


      nL = nL + 1

      ! Local stage:
      call local_stage_vp(I,Sop,G%sigma)

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
      call global_stage_vp(Top,Sop,L,I,nL)

      ! Calculate the error of the solver:
      LATIN_error= indicator1(G,L)

      print*, "LATIN ERROR:" , LATIN_error, "ITERATION:", nL
  end do

  !!$OMP END PARALLEL

  !v = linspace(0.D0,10.D0,5)





end program main
