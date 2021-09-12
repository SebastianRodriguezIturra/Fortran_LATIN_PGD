module local_stage

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






  subroutine local_stage_vp(I,Sop,sigma)

    type(Input), intent(inout) :: I
    !type(global), intent(inout) :: G
    type(Space_op), intent(inout) :: Sop
    double precision, dimension(:,:,:,:) :: sigma
    integer, dimension(4) :: siz
    double precision, dimension(Sop%comp,I%ngps,I%Nx,I%Tngpt) :: sigma_D,d_Ep,tau
    double precision, dimension(Sop%comp,I%ngps,I%Nx,I%Tngpt) :: J2_sigma,iJ2_sigma,N,nlaw

    integer :: j
    double precision, dimension(1,I%ngps,I%Nx,I%Tngpt) :: f_K

    !integer :: ind
    !integer, dimension(1,I%Nx) :: v1,v2


    !!$OMP PARALLEL


    !PRINT *, 'Hello from process: ', OMP_GET_THREAD_NUM()






    siz = (/Sop%comp,I%ngps,I%Nx,I%Tngpt/)
    L%sigma      = sigma
    sigma_D      = L%sigma - L%sigma/3

    !!$omp simd
    !do j=1,I%Tngpt
    !     sigma_D(:,:,:,j) = L%sigma(:,:,:,j) - L%sigma(:,:,:,j)/3
    !end do
    !!$omp end simd

    tau = sigma_D
    J2_sigma = reshape( sqrt(1.5*sum((tau)**2,1)), siz )


    !!$omp simd
    !do j=1,I%Tngpt
    !    iJ2_sigma(:,:,:,j) = 1/J2_sigma(:,:,:,j)
    !end do
    !!$omp end simd


    iJ2_sigma = 1/J2_sigma ! only the positive parts.
    where(J2_sigma == 0.D0) iJ2_sigma = 0.D0  !iJ2_sigma(J2_sigma==0) = 0;




    !!$omp simd
    !do j=1,I%Tngpt
    !    N(:,:,:,j) = iJ2_sigma(:,:,:,j)*tau(:,:,:,j)
    !    f_K(:,:,:,j) = J2_sigma(:,:,:,j)  - I%sigma_y
    !end do
    !!$omp end simd



    N = iJ2_sigma*tau
    f_K = J2_sigma  - I%sigma_y

    where(f_K<0.D0) f_K = 0.D0

    !stop

    !!$omp simd
    !do j=1,I%Tngpt
        !nlaw(:,:,:,j) = I%ks*(f_K(:,:,:,j)**I%np)
        !L%d_Ep(:,:,:,j)  = 3/2*(nlaw(:,:,:,j)*N(:,:,:,j))
    !end do
    !!$omp end simd



    !nlaw = I%ks*(f_K**I%np)

    nlaw = Norton_law(J2_sigma,I,Sop)
    L%d_Ep  = 3/2*(nlaw*N)

    !!$OMP END PARALLEL



  end subroutine local_stage_vp




  function Norton_law(J2_sigma,I,Sop) result(nlaw)

    implicit none
    type(Input), intent(inout) :: I
    !type(global), intent(inout) :: G
    type(Space_op), intent(inout) :: Sop
    double precision, dimension(1,I%ngps,I%Nx,I%Tngpt) :: J2_sigma,f_K,nlaw

    f_K = J2_sigma  - I%sigma_y
    where(f_K<0.D0) f_K = 0.D0
    nlaw = I%ks*(f_K**I%np)

  end function Norton_law



end module local_stage
