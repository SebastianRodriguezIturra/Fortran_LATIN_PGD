module Structure

  use Matrix_op
  implicit none

  type :: Input

     ! Select the solver: "LATIN-PGD" or "Newton-Raphson"
     character(len = 50) :: solver = "LATIN-PGD"

     ! NR solver parameters:
     double precision :: NR_error_min = 1e-4 ! norm of the residual equilibrium.
     double precision :: beta = 1.D0/4.D0
     double precision :: gamma = 1.D0/2.D0
     integer :: ntNR ! NR temporal nodes.


     ! LATIN solver parameters:
     double precision :: LATIN_error_min = 1 ! minimun value of the LATIN error [%].
     double precision :: PGD_ind = 0.1 ! PGD indicator (to stop the PGD iterations).
     double precision :: Upd_PGD_ind = 0.1 ! PGD indicator at the Preliminary step.
     integer :: PGD_iter_max = 4 ! Maximum number of PGD iterations at the global stage.
     double precision :: relax = 0.4
     double precision :: sig_factor = (1e-2)
     integer :: interval_act_modes = 1 ! number of modes to be actualized.

     logical :: PGD_eval = .true.
     double precision :: update_limit = 0.1
     double precision :: update_ind = 100
     integer :: N_up_max = 4 ! maximum number of update iterations at the preliminary step.
     integer :: up_iter = 0

     !character(len = 50) :: title="NEPE PAL QUE LEE"

     double precision :: pi = 4.D0*DATAN(1.D0)
     integer :: dim = 1
     character(len = 50) :: material = "Viscoplasticity" ! Viscoplasticity or Concrete
     integer :: value = 8.5

     ! Material parameters:

     double precision :: rho = 7800  ! Kg/m3
     double precision :: sigma_y = 8e7
     double precision :: np = 5
     double precision :: K = 1220e6
     double precision :: ks



     ! Solver Parameters:

     ! Space definition (1D for the moment):

     double precision :: E = 137600e6 !137600e3
     double precision :: Hooke
     double precision :: d1 = 0.2; ! x dimension [m]
     double precision :: d2 = 0.3; ! y dimension [m]
     double precision :: L = 4 ! z dimension [m]
     double precision :: A   ! [m2]
     integer :: Nx = 250  ! number of fem elements in x direction.
     integer :: n ! total nodes
     integer :: ns ! free nodes
     integer :: ne = 2 ! bar elemental
     double precision :: dx
     integer :: ngps = 2
     integer :: Tngps
     integer :: Nsfs = 2  ! Number of shape functions in space.
     !double precision, dimension(:), allocatable :: gps  ! vector of gauss points on a normalized temporal element.
     !double precision, dimension(:,:), allocatable :: gauss_s ! points and weight matrix.

     ! Temporal definition:

     double precision :: Tend = 10.D0
     double precision :: Ti = 0.D0
     integer :: NT = 1000
     integer :: ntc ! temporal nodes
     double precision :: dt
     double precision, dimension(:), allocatable :: t  ! time vector.
     !double precision, dimension(:), allocatable :: gpt  ! vector of gauss points on a normalized temporal element.
     double precision, dimension(:), allocatable :: vect_gpt ! total vector of gauss points in time.
     !double precision, dimension(:,:), allocatable :: gauss_t ! points and weight matrix.
     integer :: ngpt = 4
     integer :: Tngpt ! Total gauss points in time.
     character(len = 50) :: Tmethod = "Lagrange" !"Linear" ! Lagrange or Hermite
     integer :: Nsft  ! Number of shape functions in time.
     !integer :: Nsft = 2  ! Number of shape functions in time.

     ! Boundary conditions:
     double precision, dimension(:), allocatable :: ud ! imposed displacement vector.
     double precision, dimension(:), allocatable :: udNR ! imposed displacement vector.






  end type Input

  type(Input) :: I



contains



  ! Maybe constructor
  subroutine Input_assign()

    ! Newton Raphson parameters:
    !I%NR_error_min = 1e-4 ! norm of the residual equilibrium.
    !I%beta = 1/6
    !I%gamma = 1/2
    I%ntNR = I%NT + 1 ! Newton Raphson temporal nodes.


    if (I%Tmethod .eq. "Linear") then
      I%Nsft = 2
      I%ntc = I%NT + 1
    else if ( I%Tmethod .eq. "Lagrange" ) then ! Lagrange
      I%Nsft = 3
      I%ntc = I%NT*(2) + 1
    else ! Hermite
      I%Nsft = 4
      I%ntc = 2*(I%NT-1) + 4
    end if

    ! Construction and assign of the left variables:
    I%ks = (1/((I%K)**I%np))

    I%dx = I%L/I%Nx
    I%A = I%d1*I%d2
    I%Hooke = I%E  ! Only true for 1D case.

    I%Tngpt = I%NT*I%ngpt   ! Total number of integration points in time.
    I%Tngps = I%Nx*I%ngps   ! Total number of integration points in space.
    I%dt = (I%Tend - I%Ti)/I%NT




    I%n = I%Nx + 1
    I%ns = I%n - 2

    !allocate ( I%t(I%ndt) ) = linspace(I%Ti,I%Te,I%ndt)
    I%vect_gpt = linspace(I%Ti,I%Tend,I%Tngpt)
    I%t = linspace(I%Ti,I%Tend,I%nt)
    !I%ud = (I%t/I%Tend)*(I%sig_factor)*sin(2*I%pi*I%t)
    I%ud = (I%sig_factor)*sin(2*I%pi*I%vect_gpt)
    I%udNR = (I%sig_factor)*sin(2*I%pi*I%t)
    !allocate ( I%ud(I%ndt) ) = (/(i, i=0,I%Tend, I%ndt)/)


  end subroutine Input_assign









end module Structure
