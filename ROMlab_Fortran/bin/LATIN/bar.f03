module bar

  use Structure
  use Matrix_op
  use Gauss

  implicit none

  type :: Space_op

    double precision, dimension(:,:), allocatable :: gauss_s ! points and weight matrix.
    double precision, dimension(:,:), allocatable :: IntS
    double precision, dimension(:,:), allocatable :: B,Bw
    double precision, dimension(:,:), allocatable :: Be,Bwe

    double precision, dimension(:,:), allocatable :: He
    integer, dimension(:,:), allocatable :: gps
    integer, dimension(:,:), allocatable :: Idvs
    integer, dimension(:), allocatable :: dof_imposed
    integer, dimension(:), allocatable :: dof_free
    double precision, dimension(:,:), allocatable :: Melm,Kelm
    integer :: comp = 1


  end type Space_op

  type(Space_op) :: Sop


  contains



  subroutine Space_elm_matrices(I)

    type(Input), intent(inout) :: I
    integer :: iter,igp,ind
    integer, dimension(1,I%Nx) :: v1,v2
    double precision, dimension(I%ngps,1) :: waux
    double precision, dimension(I%ngps) :: x
    double precision, dimension(I%Nsfs,I%ngps) :: lx,dlx
    double precision, dimension(I%ngps,I%Nsfs) :: laux,dlaux
    double precision, dimension(I%Tngps,I%n) :: Bw

    ! real, dimension(I%ne,1) :: Felm


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !! Vector that integrate any function in the space interval Ix=[0,I.L];

    ! Uni_dimensional case
    allocate(Sop%IntS(I%Tngps,1))
    Sop%gauss_s = gauss_points_weights(I%ngps)
    do igp=1,I%ngps
        waux(igp,1)=(I%dx/2)*Sop%gauss_s(igp,2)
    end do
    !print*, "SHAPE VALUE WAUX:" , shape(waux)

    !print*, "HASTA ACA:", shape(Sop%IntS)
    !print*, I%Nx, I%dx
    !call show(waux)

    Sop%IntS=kron(ones(I%Nx,1),waux)
    !call show(ones(I%Nx,1))
    !print*, "NEPE"
    !call show(Sop%IntS)


    ! Identification matrix in space:
    Sop%Idvs = vertcat( incr(1,I%Nx), incr(2,I%Nx+1) )
    !I.Idvs = vertcat(1:I.Nx,2:I.ns)

    !                       Linear shape functions in space:

    x = ( (I%dx/2)*Sop%gauss_s(:,1)+I%dx/2 )
    !xt = (x)
    lx = vertcat(1-x/I%dx,x/I%dx) ! No derivative.
    dlx = repmat( reshape( (/-1/I%dx, 1/I%dx/), (/2,1/) ), (/1,I%ngps/) ) ! First derivative.



    !                       Linear shape functions in space:


    ! Calculation of the B operator ( eps = B*u):
    laux = transpose(lx)
    dlaux = transpose(dlx)
    allocate(Sop%B(I%Tngps,I%n))  != zeros(I.Ngps,I.ns);
    do iter = 1,I%Nx
        Sop%B( incr(1+(iter-1)*I%ngps,iter*I%ngps) , incr(iter,iter+1) ) = dlaux
    end do

    ! Rigidity matrix calculation:
    Sop%Bw = spread(Sop%IntS(:,1), 2, I%n )*Sop%B
    !I.K = (I%Hooke*I%A)*( matmul(transpose(Sop%B), Bw)

    ! Elemental operators:
    Sop%Be = dlaux
    Sop%Bwe = spread( waux(:,1), 2, I%Nsfs )*Sop%Be



    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! Elemental matrices assignation:

    allocate( Sop%Melm(2,2)  )
    allocate( Sop%Kelm(2,2)  )
    allocate( Sop%He(1,1) )

    Sop%He = I%Hooke   ! Elastic hooke tensor.
    Sop%Kelm = ( (I%Hooke*I%A)/I%dx )*reshape( (/ 1.0, -1.0, -1.0, 1.0/), (/ 2, 2 /) )  ! Elemental stiffness matrix.
    Sop%Melm = ((I%rho*I%A*I%dx)/6)*reshape( (/ 2.0, 1.0, 1.0, 2.0/), (/ 2, 2 /) )  ! Elemental mass matrix.





    ! Identification matrix in space:
    !allocate( Sop%Idvs(2,I%Nx) )
    !v1 = incr2(1,I%Nx) ; v2 = incr2(2,I%n)
    !Sop%Idvs = vertcat(v1, v2)



    ! Dofs fixed and imposed:

    allocate( Sop%dof_imposed(2)  )
    allocate( Sop%dof_free(I%ns)  )

    Sop%dof_imposed = (/1,I%n/)
    Sop%dof_free = incr(2,I%n-1)





    !I%ntc3 = 123 ! Assignation to the "I" structure.
    !print*, I%ntc3



  end subroutine Space_elm_matrices








end module bar
