module local_stage_NR

  USE OMP_LIB
  use Structure
  use Matrix_op
  use bar
  !use global_stage


  implicit none

  !type :: local

  ! double precision, dimension(:,:,:), allocatable :: sigma
  !  double precision, dimension(:,:,:), allocatable :: d_Ep

  !end type local

  !type(local) :: L


  contains


  function local_stage_vp_NR(I,Sop,sigma) result(d_Ep)

    type(Input), intent(in) :: I
    !type(global), intent(inout) :: G
    type(Space_op), intent(in) :: Sop
    double precision, dimension(:,:,:) :: sigma
    integer, dimension(3) :: siz
    double precision, dimension(Sop%comp,I%ngps,I%Nx) :: sigma_D,d_Ep,tau
    double precision, dimension(Sop%comp,I%ngps,I%Nx) :: J2_sigma,iJ2_sigma,N,nlaw

    integer :: j
    double precision, dimension(1,I%ngps,I%Nx) :: f_K


    !!$OMP PARALLEL

    siz = (/1,I%ngps,I%Nx/)
    !L%sigma      = sigma
    sigma_D      = sigma - sigma/3

    tau = sigma_D
    J2_sigma = reshape( sqrt(1.5*sum((tau)**2,1)), siz )

    iJ2_sigma = 1/J2_sigma ! only the positive parts.
    where(J2_sigma == 0.D0) iJ2_sigma = 0.D0  !iJ2_sigma(J2_sigma==0) = 0;

    N = iJ2_sigma*tau
    f_K = J2_sigma  - I%sigma_y
    where(f_K<0.D0) f_K = 0.D0

    nlaw = I%ks*(f_K**I%np)
    !nlaw = Norton_law(I,Sop)
    d_Ep  = 3/2*(nlaw*N)

    !!$OMP END PARALLEL

  end function local_stage_vp_NR


end module local_stage_NR
