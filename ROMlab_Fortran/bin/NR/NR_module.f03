module NR_module


use Structure
use bar
use Matrix_op
use gnuplot_fortran

implicit none


  type :: NR_struct

    double precision, dimension(:,:), allocatable :: U
    double precision, dimension(:,:,:,:), allocatable :: sigma,eps,d_Ep
    !double precision, dimension(:,:), allocatable :: Kss ! Mss

  end type NR_struct

  type(NR_struct) :: NR



contains

subroutine NR_solver(I,Sop)

    type(Input), intent(in) :: I
    type(Space_op), intent(in) :: Sop

    !double precision, dimension(I%n,I%n) :: M ! Total mass matrix
    double precision, dimension(I%n,I%n) :: K ! Total stifness matrix
    !double precision, dimension(I%ns,I%ns) :: Mss ! Free dofs mass matrix
    double precision, dimension(I%ns,I%ns) :: Kss ! Free dofs stifness matrix
    double precision, dimension(I%ns,I%ns) :: ChKss ! Free dofs stifness matrix

    double precision, dimension(I%ns,I%ns) :: invKss ! inverse of Kss
    double precision, dimension(I%ns,2) :: Kud ! stiffness matrix related to fix dofs
    double precision, dimension(I%ns,I%Tngpt) :: Fud,x
    double precision, dimension(I%n,I%Tngpt) :: U   ! Displacement
    double precision, dimension(2,I%Tngpt) :: Ud   ! imposed Displacement
    double precision, dimension(I%ns,I%Tngpt) :: Us ! free dofs displacement
    double precision, dimension(I%Tngpt) :: vect
    double precision, dimension(I%Tngps,I%Tngpt) :: eps,sigma,d_Ep
    integer, dimension(4) :: siz
    integer :: itemp


    integer :: ielm
    integer :: ns, nrhs, info, lda, ldb   ! dimension of the matrix that must be inverted.
    integer, allocatable, dimension(:) :: ipiv

    ! Assignation of the final solution to structure NR:
    allocate(NR%U(I%n,I%Tngpt))
    allocate(NR%eps(siz(1),siz(2),siz(3),siz(4)))
    allocate(NR%sigma(siz(1),siz(2),siz(3),siz(4)))

    siz = (/Sop%comp,I%ngps,I%Nx,I%Tngpt/)

    ! Assembly:
    do ielm = 1, I%Nx
       !M( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm)  ) = M( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm) ) + Sop%Melm
       K( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm)  ) = K( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm) ) + Sop%Kelm
    end do


    ! Imposed displacement:
    U(1,:) = 0 ! Fix condition.
    U(I%n,:) = I%ud ! displacement condition.


    ! Dofs free and imposed of the matrices used:
    Kss=K(Sop%dof_free,Sop%dof_free)
    Kud=K(Sop%dof_free,Sop%dof_imposed)
    ChKss = cholesky(Kss) ! Cholesky decomposition of Kss


    ! Imposed displacement and generated force.
    Ud = U(Sop%dof_imposed,:)
    Fud=matmul(Kud,Ud) ! Force due to the imposition of displacements.


    if (.false.) then ! Elastic solution

      !! Resolution (code Blas Lapack):
      ns = size(Kss,1)
      nrhs = size(Fud,2)
      !print*, ns,nrhs
      lda = ns; ldb = ns
      allocate(ipiv(ns))
      x = -Fud
      call dgesv(ns,nrhs,Kss,lda,ipiv,x,ldb,info)
      Us = x


      ! Displacement, strain and stress:
      U(Sop%dof_free,:) = Us
      eps = matmul(Sop%B,U)
      sigma = reshape(matmul(Sop%He, reshape( eps, (/Sop%comp, size(eps)/Sop%comp/))), (/size(eps,1),I%Tngpt/) )

      print*, "SOLUTION DONE, ERROR:", 100*norm2( matmul(Kss,Us) + Fud )/norm2(Fud) , "[%]."

    else ! Elasto-visco-plasticity

      do itemp = 1,I%Tngpt

        ! Initial solution:


        do while (I%LATIN_error_min <= LATIN_error )

          ! Evaluation of the nonlinear constitutive relation:
          call local_stage_vp_NR(I,Sop,sigma)

          ! Linear correction:



        end do ! Iterative process of the NR

      end do ! temporal for

    end if


    NR%U = U
    NR%eps = reshape(eps, siz )    !matmul(Sop%B,Elast%U) ! Sop%B*Elast%U
    NR%sigma = reshape(sigma, siz )  !matmul(Sop%He,Elast%eps)  ! Sop%He*Elast%eps



  end subroutine NR_solver

end module NR_module
