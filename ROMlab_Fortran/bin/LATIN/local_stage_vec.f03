module local_stage_vec

  USE OMP_LIB
  use Structure
  use Matrix_op
  use bar
  !use global_stage



  implicit none

  type :: local

    double precision, dimension(:,:,:,:), allocatable :: sigma
    double precision, dimension(:,:,:,:), allocatable :: d_Ep

  end type local

  type(local) :: L


  contains






  subroutine local_stage_vp_vec(I,Sop,sigma)

    type(Input), intent(inout) :: I
    !type(global), intent(inout) :: G
    type(Space_op), intent(inout) :: Sop
    double precision, dimension(:,:,:,:) :: sigma
    integer, dimension(4) :: siz
    double precision, dimension(Sop%comp,I%ngps,I%Nx,I%Tngpt) :: sigma_D,d_Ep,tau
    double precision, dimension(1,I%ngps,I%Nx,I%Tngpt) :: J2_sigma,iJ2_sigma,N,nlaw

    integer :: i1,j,k,l1
    double precision, dimension(1,I%ngps,I%Nx,I%Tngpt) :: f_K

    !integer :: ind
    !integer, dimension(1,I%Nx) :: v1,v2


    !!$OMP PARALLEL


    !PRINT *, 'Hello from process: ', OMP_GET_THREAD_NUM()


    siz = (/Sop%comp,I%ngps,I%Nx,I%Tngpt/)
    L%sigma      = sigma
    L%d_Ep = sigma
    !allocate(L%d_Ep( siz(1),siz(2),siz(3),siz(4) ))


    !$omp simd
    do i1=1,Sop%comp
      do j=1,I%ngps
        do k=1,I%Nx
          do l1=1,I%Tngpt
            sigma_D(i1,j,k,l1) = L%sigma(i1,j,k,l1) - L%sigma(i1,j,k,l1)/3
            tau(i1,j,k,l1) = sigma_D(i1,j,k,l1)

          end do
        end do
      end do
    end do
    !$omp end simd


    J2_sigma = reshape( sqrt(1.5*sum((tau)**2,1)), siz )



    !$omp simd
    do j=1,I%ngps
      do k=1,I%Nx
        do l1=1,I%Tngpt
          iJ2_sigma(1,j,k,l1) = 1/J2_sigma(1,j,k,l1)
        end do
      end do
    end do
    !$omp end simd


    where(J2_sigma == 0.D0) iJ2_sigma = 0.D0  !iJ2_sigma(J2_sigma==0) = 0;


    !$omp simd
    do i1=1,Sop%comp
      do j=1,I%ngps
        do k=1,I%Nx
          do l1=1,I%Tngpt
            N(i1,j,k,l1) = iJ2_sigma(1,j,k,l1)*tau(i1,j,k,l1)
          end do
        end do
      end do
    end do
    !$omp end simd


    !$omp simd
    do j=1,I%ngps
      do k=1,I%Nx
        do l1=1,I%Tngpt
          f_K(1,j,k,l1) = J2_sigma(1,j,k,l1)  - I%sigma_y
        end do
      end do
    end do
    !$omp end simd

    where(f_K<0.D0) f_K = 0.D0

    !$omp simd
    do j=1,I%ngps
      do k=1,I%Nx
        do l1=1,I%Tngpt
          nlaw(1,j,k,l1) = I%ks*(f_K(1,j,k,l1)**I%np)
        end do
      end do
    end do
    !$omp end simd



    !$omp simd
    do i1=1,Sop%comp
      do j=1,I%ngps
        do k=1,I%Nx
          do l1=1,I%Tngpt
            L%d_Ep(i1,j,k,l1)  = 3/2*(nlaw(1,j,k,l1)*N(i1,j,k,l1))
          end do
        end do
      end do
    end do
    !$omp end simd
    !stop



    !!$OMP END PARALLEL



  end subroutine local_stage_vp_vec







end module local_stage_vec
