module NR_module


use Structure
use bar
use Matrix_op
use gnuplot_fortran
use local_stage_NR


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
    double precision, dimension(I%ns,I%ntNR) :: Fud,x
    double precision, dimension(I%ns,I%ntNR-1) :: dfud

    double precision, dimension(I%n,I%ntNR) :: U   ! Displacement
    double precision, dimension(I%n,I%ntNR) :: Ut   ! Displacement
    double precision, dimension(I%n,I%ntNR) :: Utt   ! Displacement

    double precision, dimension(2,I%ntNR) :: Ud   ! imposed Displacement

    double precision, dimension(I%ns,I%ntNR) :: Us ! free dofs displacement
    double precision, dimension(I%ns,I%ntNR) :: Ust ! free dofs displacement
    double precision, dimension(I%ns,I%ntNR) :: Ustt ! free dofs displacement

    double precision, dimension(I%ntNR) :: vect
    double precision, dimension(Sop%comp,I%ngps,I%Nx,I%ntNR) :: eps,sigma,d_Ep,Ep !(I%Tngps,I%ntNR)

    integer, dimension(4) :: siz
    integer, dimension(3) :: sizr
    integer :: iter, iT
    double precision :: error

    double precision, dimension(I%ns,I%ntNR) :: g
    double precision, dimension(I%ns,1) :: res, dUs
    double precision, dimension(I%n,1) :: F
    double precision, dimension(I%ns) :: xs

    integer, dimension(I%ntNR-1) ::  vec1, vec2

    double precision, dimension(Sop%comp,I%ngps,I%Nx) :: sigma_a, eps_a, d_Epa 

    integer :: ielm
    integer :: ns, nrhs, info, lda, ldb   ! dimension of the matrix that must be inverted.
    integer, allocatable, dimension(:) :: ipiv


    ! Assignation of the final solution to structure NR:
    siz = (/Sop%comp,I%ngps,I%Nx,I%ntNR/)
    sizr = (/Sop%comp,I%ngps,I%Nx/)

    allocate(NR%U(I%n,I%ntNR))
    allocate(NR%eps(siz(1),siz(2),siz(3),siz(4)))
    allocate(NR%sigma(siz(1),siz(2),siz(3),siz(4)))
    allocate(NR%d_Ep(siz(1),siz(2),siz(3),siz(4)))

    ! Init:
    eps = 0; sigma = 0; d_Ep = 0; Ep = 0; U = 0;
    Us = 0; Ust = 0; Ustt = 0; 
    g = 0; res = 0; dUs = 0;

    ! Assembly:
    do ielm = 1, I%Nx
       !M( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm)  ) = M( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm) ) + Sop%Melm
       K( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm)  ) = K( Sop%Idvs(:,ielm) , Sop%Idvs(:,ielm) ) + Sop%Kelm
    end do


    ! Imposed displacement:
    U(1,:) = 0 ! Fix condition.
    U(I%n,:) = I%udNR ! displacement condition.

    !print*, I%udNR
    !stop

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
      eps = reshape( matmul(Sop%B,U) , siz)
      sigma = reshape(matmul(Sop%He, reshape( eps, (/Sop%comp, size(eps)/Sop%comp/))), siz )
      !sigma = reshape(matmul(Sop%He, reshape( eps, (/Sop%comp, size(eps)/Sop%comp/))), (/size(eps,1),I%ntNR/) )

      NR%U = U
      NR%eps = eps    
      NR%sigma = sigma    
      !NR%eps = reshape(eps, siz )    
      !NR%sigma = reshape(sigma, siz )  


      print*, "SOLUTION DONE, ERROR:", 100*norm2( matmul(Kss,Us) + Fud )/norm2(Fud) , "[%]."

    else ! Elasto-visco-plasticity

      vec1 = incr(2,I%ntNR)
      vec2 = incr(1,I%ntNR-1)
      dfud = Fud(:,vec1) - Fud(:,vec2);

      !print*, dfud
      !stop

      do iT = 2,I%ntNR

          d_Epa=0; ! derivative in time of the auxiliary plasticity function (Initialization).
          !iT=iT+1
          Ustt(:,iT)=0;
          Us(:,iT) = Us(:,iT-1) + I%dt*Ust(:,iT-1) + (0.5-I%beta)*(I%dt**2)*Ustt(:,iT-1) +(I%dt**2)*I%beta*Ustt(:,iT);
          Ust(:,iT) = Ust(:,iT-1) + (1-I%gamma)*I%dt*Ustt(:,iT-1) + I%dt*I%gamma*Ustt(:,iT);

          ! (Initialization) Boundary conditions imposed:
          res(:,1) = g(:,iT-1);
          !du=transpose(LM)\(LM\(-res - dfud(:,iT-1)));
          xs = -res(:,1) - dfud(:,iT-1)
          call inv_chol(ChKss,xs)
          dUs(:,1) = xs

          !print*, dUs
          !stop  

          Us(:,iT)=Us(:,iT)+dUs(:,1);
          Ust(:,iT)=Ust(:,iT)+dUs(:,1)*(I%gamma/(I%beta*I%dt));
          Ustt(:,iT)=Ustt(:,iT)+dUs(:,1)*(1/(I%beta*I%dt**2));


          error=100; iter=0;

          do while (I%NR_error_min <= error )

            iter=iter+1

            if (iter>1) then
                xs = -res(:,1)
                call inv_chol(ChKss,xs)
                dUs(:,1) = xs

                Us(:,iT) = Us(:,iT) + dUs(:,1);
                Ust(:,iT) = Ust(:,iT) + dUs(:,1)*(I%gamma/(I%beta*I%dt));
                Ustt(:,iT) = Ustt(:,iT) + dUs(:,1)*(1/(I%beta*I%dt**2));
            end if

            U(Sop%dof_free,iT) = Us(:,iT)
            eps_a = reshape( matmul(Sop%B,U(:,iT)), sizr )

           
            ! Stress calculation:
            sigma_a = I%Hooke*(eps_a-(Ep(:,:,:,iT-1)+I%dt*d_Epa))

            
            ! Derivative in time of the plasticity.
            d_Epa = local_stage_vp_NR(I,Sop,sigma_a)

           

            ! Calculation of g:
            F = I%A*( matmul( transpose(Sop%Bw) , reshape(sigma_a, (/size(sigma_a),1/)  ) ) )
            g(:,iT) = F(Sop%dof_free,1)

           

            ! Residual function:
            res(:,1) = g(:,iT)
            !error=sqrt(transpose(res)*res);
            error = norm2(res)

            !print*, error
            !stop

            !print*, error

          end do ! while

          print*, iter, error
          stop

          eps(:,:,:,iT)=eps_a;  ! total deformation.
          Ep(:,:,:,iT) = Ep(:,:,:,iT-1) + I%dt*d_Epa; ! Riemman integration.
          d_Ep(:,:,:,iT)=d_Epa;
          sigma(:,:,:,iT)=sigma_a;

          print*, iT

        end do ! end for iT<=nt

        

        ! Save the results on NR:
        NR%U = U
        NR%eps = eps    !matmul(Sop%B,Elast%U) ! Sop%B*Elast%U
        NR%d_Ep = d_Ep
        NR%sigma = sigma  !matmul(Sop%He,Elast%eps)  ! Sop%He*Elast%eps  
        
        ! Error determination:
        print*, "SOLUTION DONE, ERROR:", norm2( g )

        !NR%eps = reshape(eps, siz)    !matmul(Sop%B,Elast%U) ! Sop%B*Elast%U
        !NR%d_Ep = reshape(d_Ep, siz)
        !NR%sigma = reshape(sigma, siz)  !matmul(Sop%He,Elast%eps)  ! Sop%He*Elast%eps


    end if

  end subroutine NR_solver

end module NR_module
