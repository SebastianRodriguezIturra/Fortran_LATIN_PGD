module Elastic_sol_PGD


use Structure
use bar
use Matrix_op
use gnuplot_fortran
use Time_matrices

implicit none


  type :: Init_sol

    double precision, dimension(:,:), allocatable :: U
    double precision, dimension(:,:,:,:), allocatable :: sigma,eps
    double precision, dimension(:,:), allocatable :: Mss,Kss,ChKss


  end type Init_sol

  type(Init_sol) :: Elast



contains

subroutine solve_init_PGD(I,Sop,Top)

    type(Input), intent(in) :: I
    type(Space_op), intent(in) :: Sop
    type(Time_op), intent(inout) :: Top


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


    ! PGD variable definitions:
    double precision, dimension(I%ns,I%Tngpt) :: F
    double precision :: error_PGD_lim, error_PGD, Fud_norm, a, b, snorm, stagn_lim, stagn
    integer ::  PGDm, iter_pgd, max_iter

    double precision, dimension(I%ns,1) :: Si,Fi
    double precision, dimension(I%Tngpt,1) :: Ti, Ta
    double precision, dimension(:,:) , allocatable :: Sl,Tl,S,T


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

      print*, "Classic resolution of elastic problem."

      !! (code Blas Lapack):
      ns = size(Kss,1)
      nrhs = size(Fud,2)
      !print*, ns,nrhs
      lda = ns; ldb = ns
      allocate(ipiv(ns))
      x = -Fud
      call dgesv(ns,nrhs,Kss,lda,ipiv,x,ldb,info)
      Us = x


    else

      ! PGD resolution:

      print*, "PGD resolution of elastic problem."

      ns = size(Kss,1)
      nrhs = 1
      lda = ns; ldb = ns
      allocate(ipiv(ns)); !ipiv = 0


      error_PGD = 100
      error_PGD_lim = 0.01
      stagn_lim = 0.1

      max_iter = 6
      Fud_norm = norm2(Fud)

      ! S, T;  PGD functions related to the displacement solution.
      PGDm = 0 ! PGD mode
      do while ( error_PGD_lim <= error_PGD )

        PGDm = PGDm + 1 ! PGD mode


        ! Update of the external force:
        if (PGDm .eq. 1) then
          F = -Fud
        else
          F = -Fud - matmul( matmul(Kss,S), transpose(T) )
        end if
        ! end update of the external force:

        Si = 0; Ti = reshape(Top%vect_gpt, (/I%Tngpt,1/) ); Ta = Ti ! Initialization.

        iter_pgd = 0
        stagn = 100
        do while ( stagn_lim <= stagn .and. max_iter > iter_pgd )

          iter_pgd = iter_pgd + 1
          !print*, "Iter PGD", iter_pgd

          ! Space problem:
          a = sum( (Ti**2)*Top%IntT )
          Fi = matmul(F,Ti*Top%IntT)
          call dgesv(ns,nrhs,Kss,lda,ipiv,Fi,ldb,info)
          Si = Fi/a


          ! Temporal problem:
          b = sum(matmul( matmul(transpose(Si),Kss), Si ))
          Ti = transpose( matmul( transpose(Si) , F ) )/b


          ! Stagnation of the PGD:
          snorm = norm2(Si)
          Si = Si/snorm; Ti = Ti*snorm
          stagn = abs( norm2(Ti-Ta)/norm2(Ta) )
          Ta = Ti

          !print*, "STAGN" ,stagn


        end do ! end while PGD calculation


        ! Update of the PGD and external excitation:

        if (PGDm .eq. 1) then

          allocate(Sl(ns,PGDm)); allocate(Tl(I%Tngpt,PGDm))
          allocate(S(ns,PGDm)); allocate(T(I%Tngpt,PGDm))

          Sl = Si; Tl = Ti
          S = Si; T = Ti

        else

          deallocate(Sl); deallocate(Tl)
          allocate(Sl(ns,PGDm-1)); allocate(Tl(I%Tngpt,PGDm-1))
          Sl = S; Tl = T

          deallocate(S); deallocate(T)
          allocate(S(ns,PGDm)); allocate(T(I%Tngpt,PGDm))
          S(:,1:PGDm-1) = Sl; T(:,1:PGDm-1) = Tl
          S(:,PGDm) = Si(:,1); T(:,PGDm) = Ti(:,1)

        end if

        ! PGD indicator:
        error_PGD = 100*norm2( matmul( matmul( Kss,S ), transpose(T) )  + Fud)/Fud_norm ! [%]
        print*, "Error PGD", error_PGD

      end do ! end while PGD process

      Us = matmul( S,transpose(T) )

    end if ! end solution choices


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



  end subroutine solve_init_PGD




end module Elastic_sol_PGD
