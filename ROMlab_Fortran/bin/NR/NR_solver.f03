module Elastic_sol


use Structure
use bar
use Matrix_op
use gnuplot_fortran

implicit none


  type :: Init_sol

    double precision, dimension(:,:), allocatable :: U
    double precision, dimension(:,:,:,:), allocatable :: sigma,eps
    double precision, dimension(:,:), allocatable :: Mss,Kss,ChKss


  end type Init_sol

  type(Init_sol) :: Elast



contains

subroutine solve_init(I,Sop)

    type(Input), intent(in) :: I
    type(Space_op), intent(in) :: Sop

    double precision, dimension(I%n,I%n) :: M ! Total mass matrix
    double precision, dimension(I%n,I%n) :: K ! Total stifness matrix
    double precision, dimension(I%ns,I%ns) :: Mss ! Free dofs mass matrix
    double precision, dimension(I%ns,I%ns) :: Kss ! Free dofs stifness matrix
    double precision, dimension(I%ns,I%ns) :: ChKss ! Free dofs stifness matrix

    double precision, dimension(I%ns,I%ns) :: invKss ! inverse of Kss
    double precision, dimension(I%ns,2) :: Kud ! stiffness matrix related to fix dofs
    double precision, dimension(I%ns,I%Tngpt) :: Fud,x
    double precision, dimension(I%n,I%Tngpt) :: U   ! Displacement
    double precision, dimension(2,I%Tngpt) :: Ud   ! imposed Displacement
    double precision, dimension(I%ns,I%Tngpt) :: Us ! free dofs displacement
    double precision, dimension(I%Tngpt) :: vect
    double precision, dimension(I%Tngps,I%Tngpt) :: eps,sigma
    integer, dimension(4) :: siz

    integer :: ielm
    integer :: ns, nrhs, info, lda, ldb   ! dimension of the matrix that must be inverted.
    integer, allocatable, dimension(:) :: ipiv

    ! Assembly:
    do ielm = 1, I%Nx
       M( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm)  ) = M( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm) ) + Sop%Melm
       K( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm)  ) = K( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm) ) + Sop%Kelm
    end do




    ! Imposed displacement:
    U(1,:) = 0 ! Fix condition.
    U(I%n,:) = I%ud ! displacement condition.


    ! Dofs free and imposed of the matrices used:
    Kss=K(Sop%dof_free,Sop%dof_free)
    Kud=K(Sop%dof_free,Sop%dof_imposed)

    ChKss = cholesky(Kss) ! Cholesky decomposition of Kss



    Ud = U(Sop%dof_imposed,:)
    Fud=matmul(Kud,Ud) ! Force due to the imposition of displacements.






    !! Classic solution (code Seba):

    if (.false.) then
      invKss = inv(Kss)
      Us = matmul(invKss,-Fud)

    else
      !! Classic solution (code Blas Lapack):
      ns = size(Kss,1)
      nrhs = size(Fud,2)
      !print*, ns,nrhs
      lda = ns; ldb = ns
      allocate(ipiv(ns))
      x = -Fud
      call dgesv(ns,nrhs,Kss,lda,ipiv,x,ldb,info)
      Us = x
      !print*, "info:"
      !print*, info

    end if

    ! Plots of the solution:
    !vect = Us(5,:)
    !call plot(vect)
    !!call plot(incr(1,10))

    !print*, "Reshape of the multi-dimnsional vector:"
    !print*, shape( reshape(Us, (/ size(Us), 1 /) ) )

    !print*, "norm2 value:"
    !print*, norm2(Us)
    !print*, shape(Us)
    !call show(Us)

    U(Sop%dof_free,:) = Us
    eps = matmul(Sop%B,U)
    sigma = reshape(matmul(Sop%He, reshape( eps, (/Sop%comp, size(eps)/Sop%comp/))), (/size(eps,1),I%Tngpt/) )




    ! Assignation of the final solution:
    siz = (/Sop%comp,I%ngps,I%Nx,I%Tngpt/)
    Elast%U = U
    allocate(Elast%eps(siz(1),siz(2),siz(3),siz(4)))
    allocate(Elast%sigma(siz(1),siz(2),siz(3),siz(4)))


    Elast%eps = reshape(eps, siz )    !matmul(Sop%B,Elast%U) ! Sop%B*Elast%U
    Elast%sigma = reshape(sigma, siz )  !matmul(Sop%He,Elast%eps)  ! Sop%He*Elast%eps

    Elast%Mss = Mss ! mass matrix
    Elast%Kss = Kss ! stiffness matrix
    Elast%ChKss = ChKss ! stiffness matrix


    print*, "ELASTIC SOLUTION DONE, ERROR:", 100*norm2( matmul(Kss,Us) + Fud )/norm2(Fud) , "[%]."



  end subroutine solve_init




end module Elastic_sol
