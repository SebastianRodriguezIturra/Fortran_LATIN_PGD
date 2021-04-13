
!!!!!!!!!!!!!!!!!!!!!!!!!!! Preliminary step !!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Preliminary_step(H,Q,Sop,Top,I,T,dT)  !result(alfa,alfat) ! [alfa,alfat]

! S,D,E,H,Q,Sop,Top,I


  implicit none

  type(Input), intent(inout) :: I
  type(Space_op), intent(inout) :: Sop
  type(Time_op), intent(inout) :: Top

  double precision, dimension(:,:) :: H  ! constant operator.
  double precision , dimension(:,:), allocatable :: Qr ! reshape of Q
  double precision , dimension(:,:,:,:) :: Q
  double precision, dimension(:,:,:,:) :: S,D  ! S: Plastic deformation tensor; D: Stress tensor.

  ! vector determination:
  double precision, dimension(:,:), allocatable :: aux ! auxiliar variable for vetors calculation.
  double precision, dimension(:,:), allocatable :: Sr
  !double precision, dimension(:,:,:,:), allocatable :: Sg,Dg,auxt
  double precision, dimension(:,:), allocatable :: D1, D0 !

  ! constants determination:
  double precision, dimension(:,:), allocatable :: aux1 ! auxiliar variable for constants calculation.
  double precision, dimension(:,:), allocatable :: A11,A10,A01,A00 ! constants factors.



  !double precision, dimension(I%ntc-1,I%ntc-1) :: Mt
  !double precision, dimension(I%ntc,1) :: Fti  ! assembly vector
  !double precision, dimension(I%ntc-1,1) :: Ft ! assembly vector

  double precision, dimension(I%ngpt,1) :: IntTe
  !double precision, dimension(:,:), allocatable :: IntSg
  double precision, dimension(I%Nsft,I%NT) :: Felm1, Felm2  ! elemental vector
  double precision, dimension(:,:,:), allocatable :: Felm


  double precision, dimension(:,:), allocatable :: Mt
  double precision, dimension(:,:), allocatable :: Ft ! assembly vector
  double precision, dimension(I%ntc,1) :: Fti  ! assembly vector
  !double precision, dimension(:,:), allocatable :: Fti  ! assembly vector



  integer, dimension(4) :: siz
  integer, dimension(3) :: sizr
  integer :: Nelms
  integer :: nbmi, nbme, dnbm


  integer :: i1,j1
  integer :: ns, nrhs, lda, ldb, info
  double precision, dimension(:), allocatable :: ipiv, x
  integer , dimension(I%Nsft*I%NT) :: IDvtc


  double precision, dimension(:,:), allocatable :: alfa_fem
  double precision, dimension(:,:), allocatable :: alfac
  double precision, dimension(:,:), allocatable :: alfacc
  double precision, dimension(:,:), allocatable :: T, dT


  ! multiplication quantities blas:
  double precision, dimension(I%Nsft,I%ngpt) :: a1,a2
  double precision, dimension(I%ngpt, I%NT) :: b1,b2


  ! Assignation of the spatial functions from global_stage:

  ! end Assignation of the spatial functions from global_stage:


  siz = shape(S); sizr = siz(1:3)



  ! Determination  of the modes that need to be actualized
  if (siz(4)==I%interval_act_modes) then
    nbmi = 1
    nbme = I%interval_act_modes

  else if (siz(4)>I%interval_act_modes) then
    nbmi = siz(4)-I%interval_act_modes + 1
    nbme = siz(4)

  else  ! siz(4) < I%interval_act_modes
    nbmi = 1
    nbme = siz(4)
  end if

  dnbm = nbme - nbmi ! number of modes that will be actualized.
  ! End Determination  of the modes that need to be actualized


  S = G%d_Ep(:,:,:,nbmi:nbme)  ! Plastic deformation tensor.
  D = G%sigma(:,:,:,nbmi:nbme) ! Stress tensor.

  ! Constants factors determination:

  allocate(A11(dnbm,dnbm)); allocate(A10(dnbm,dnbm))
  allocate(A01(dnbm,dnbm)); allocate(A00(dnbm,dnbm))

  Nelms = size(S(:,:,:,1)) ! Number of elements in space.

  allocate( aux1( Nelms , 1 ) )
  aux1 = reshape( matmul(H, reshape(D, (/Sop%comp, Nelms/Sop%comp /) ) ), (/Nelms,dnbm/) )

  do i1=1,dnbm
    do j1=1,dnbm
      A11 =  sum( (reshape(S(:,:,:,i1)*S(:,:,:,j1), (/ Nelms,1 /) ))*Sop%IntS )  !matmul( transpose(reshape(S, (/ size(S),1 /) )**2) ,  Sop%IntS )
      A10 = -sum( reshape(S(:,:,:,i1), (/ Nelms,1 /) )*aux1(:,j1)*Sop%IntS )  !matmul( transpose(reshape(S, (/ size(S),1 /) )*reshape(H*D,  (/ size(D),1 /)  )) , Sop%IntS )
      A01 = -sum( aux1(:,i1)*reshape(S(:,:,:,j1), (/ Nelms,1 /) )*Sop%IntS )
      A00 =  sum( (aux1(:,i1)*aux1(:,j1))*Sop%IntS )  !matmul( transpose(reshape(H*D,  (/ size(D),1 /) )**2) , Sop%IntS )
    end do
  end do


  ! Temporal functions of the rhs determination:

  siz = shape(Q);

  allocate(D1(I%Tngpt,dnbm)); allocate(D0(I%Tngpt,dnbm))


  allocate(Qr( size(Q)/I%Tngpt , I%Tngpt ))
  allocate(Sr( size(Q)/I%Tngpt , dnbm ))
  allocate(aux( size(Q)/I%Tngpt , dnbm ))

  aux = reshape( matmul(H, reshape( D,  (/Sop%comp, size(D)/Sop%comp /) ) ) , (/size(Q)/I%Tngpt , dnbm/)  )


  do it = 1,I%Tngpt
    do im = 1,dnbm
      D0(it,im) =  sum( Sop%IntS*aux(:,im)*Qr(:,it),1)
      D1(it,im) = -sum( Sop%IntS*Sr(:,im)*Qr(:,it),1)
    end do
  end do

  !call plot(D0(:,1))
  !call plot(D1(:,1))
  !stop





  ! DETERMINATION OF THE FEM MATRIX AND VECTOR:

  allocate(Mt( (I%ntc-1)*dnbm , (I%ntc-1)*dnbm ))
  Mt = kron(Top%I11Tr,A11) + kron(Top%I01Tr,A01) + kron(transpose(Top%I01Tr),A10) + kron(Top%I00Tr,A00)
  !Mt = (Top%I11Tr)*A11 + (Top%I01Tr)*A10 + (Top%I10Tr)*A10 + (Top%I00Tr)*A00


  ! calculate the elementary finite vectors of rhs
  if (.true.) then

    a1 = Top%vwalfat
    b1 = reshape(D1,  (/ I%ngpt, I%NT*dnbm /) )
    a2 = Top%vwalfa
    b2 = reshape(D0,  (/ I%ngpt, I%NT*dnbm /) )


    !CALL DGEMM (transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
    call dgemm('n','n',size(Felm1,1),size(Felm1,2),size(a1,2),1.D0,a1,size(a1,1),b1,size(b1,1),0.D0,Felm1,size(Felm1,1))
    call dgemm('n','n',size(Felm2,1),size(Felm2,2),size(a2,2),1.D0,a2,size(a2,1),b2,size(b2,1),0.D0,Felm2,size(Felm2,1))

  else

    Felm1 = matmul( Top%vwalfat , reshape(D1,  (/ I%ngpt, I%NT*dnbm /) )  )   ! N shape functions x ngpt
    Felm2 = matmul( Top%vwalfa  , reshape(D0,  (/ I%ngpt, I%NT*dnbm /) )  ) ! N shape functions x ngpt

  end if

  allocate(Felm( I%Nsft,I%NT, dnbm ))
  Felm =  reshape( Felm1 + Felm2,  (/ I%Nsft,I%NT,dnbm /) )


  allocate(Fti( I%ntc , dnbm ))
  allocate(Ft( (I%ntc-1)*dnbm , 1 )) ! rhs vector


  Felm = Felm1 + Felm2  ! Elemental vectors.
  Fti = 0; Ft = 0
  do j=1,I%NT
    do im=1,dnbm
      Fti(Top%IDvtc(:,j),im) = Fti(Top%IDvtc(:,j),im) + Felm(:,j,im)
    end do
  end do
  Ft =  reshape( Fti(Top%dof_free,:), (/(I%ntc-1)*dnbm,1/) ) ! Null initial condition.


  ! FEM RESOLUTION:

  allocate( alfa_fem( I%ntc, dnbm ) )


  !!! Fast inverse matrix:
  !ns, nrhs, lda, ldb, ipiv, x
  ns = size(Mt,1)
  nrhs = size(Ft,2)
  lda = ns; ldb = ns
  allocate(ipiv((I%ntc-1)*dnbm)); allocate(x((I%ntc-1)*dnbm))
  x = Ft(:,1)
  call dgesv(ns,nrhs,Mt,lda,ipiv,x,ldb,info)
  alfa_fem(Top%dof_free,:) = reshape(x, (/ I%ntc-1,dnbm /) )
  !!! end Fast inverse matrix:








  IDvtc = reshape(Top%IDvtc, (/ size(Top%IDvtc) /)  )

  allocate( alfac( I%Nsft*I%NT,dnbm ) )
  alfac = reshape(alfa_fem(IDvtc,:), (/I%Nsft*I%NT,dnbm/) )

  allocate( alfacc( I%Nsft,I%NT*dnbm ) )
  alfacc=reshape( alfac,  (/I%Nsft,I%NT*dnbm/) )


  allocate(T(I%Tngpt,dnbm))
  allocate(dT(I%Tngpt,dnbm))

  T = reshape( matmul( transpose(Top%valfa) , alfacc ) , (/I%Tngpt,dnbm/) )  ! alfa in each gauss point
  dT =  reshape( matmul(  transpose(Top%valfat), alfacc ), (/I%Tngpt,dnbm/) )  ! alfa in each gauss point

  !call plot(T(:,1))
  !stop

  !call plot(dT(:,1))
  !stop


end subroutine Preliminary_step




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
