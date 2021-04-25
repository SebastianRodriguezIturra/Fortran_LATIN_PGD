
module global_stage

  use Structure
  use Matrix_op
  use local_stage
  !use local_stage_vec
  use Elastic_sol
  use bar
  use Time_matrices

  implicit none

  type :: global

    ! Global stage elastic initial solution:

    double precision, dimension(:,:), allocatable :: Utot_o ! Displacement.
    double precision, dimension(:,:,:,:), allocatable :: eps_o ! Strain.
    double precision, dimension(:,:,:,:), allocatable :: sigma_o ! Stress.
    double precision, dimension(:,:,:,:), allocatable :: d_Ep_o ! Plastic deformation.


    ! Global stage general solution:
    double precision, dimension(:,:), allocatable :: Utot ! Displacement.
    double precision, dimension(:,:,:,:), allocatable :: eps ! Strain.
    double precision, dimension(:,:,:,:), allocatable :: sigma ! Stress.
    double precision, dimension(:,:,:,:), allocatable :: d_Ep ! Plastic deformation.

    ! PGD modes:
    double precision, dimension(:,:), allocatable :: U, T_s, dT_s ! Space and Temporal modes
    double precision, dimension(:,:,:,:), allocatable :: S_eps, S_Sigma, S_d_Ep ! Spatial modes.

    ! Constants factors:
    logical :: saved_const = .false.
    double precision, dimension(:,:), allocatable :: A11,A10,A01,A00

  end type global

  type(global) :: G


  contains


subroutine assign_elastic(Elast,I)
  implicit none
  type(Init_sol), intent(inout) :: Elast
  type(Input), intent(inout) :: I

  ! Initial elastic solution:
  allocate(G%Utot_o( I%n, I%Tngpt ))
  allocate(G%eps_o( Sop%comp,I%ngps,I%Nx, I%Tngpt ))
  allocate(G%sigma_o( Sop%comp,I%ngps,I%Nx, I%Tngpt ))
  allocate(G%d_Ep_o( Sop%comp,I%ngps,I%Nx, I%Tngpt ))

  ! Initialize the general global stage solution:
  allocate(G%Utot( I%n, I%Tngpt ))
  allocate(G%eps( Sop%comp,I%ngps,I%Nx, I%Tngpt ))
  allocate(G%sigma( Sop%comp,I%ngps,I%Nx, I%Tngpt ))
  allocate(G%d_Ep( Sop%comp,I%ngps,I%Nx, I%Tngpt ))


  G%Utot_o = Elast%U
  G%eps_o = Elast%eps
  G%sigma_o = Elast%sigma
  G%d_Ep_o = 0


  G%Utot = Elast%U
  G%eps = Elast%eps
  G%sigma = Elast%sigma
  G%d_Ep = 0

  !print*, "BLAAAAA" ,shape(G%d_Ep)

  allocate(G%A11(1,1))
  allocate(G%A10(1,1))
  allocate(G%A01(1,1))
  allocate(G%A00(1,1))

  G%A11 = 0
  G%A10 = 0
  G%A01 = 0
  G%A00 = 0

end subroutine assign_elastic





subroutine global_stage_vp(Elast,Top,Sop,L,I,iter)  ! [global_stage]=Global_stage_plasticity(global_stage,local,I,iter)

  type(Init_sol), intent(inout) :: Elast
  type(Input), intent(inout) :: I
  type(local), intent(inout) :: L
  type(Space_op), intent(inout) :: Sop
  type(Time_op), intent(inout) :: Top

  integer :: iter,iter_pgd
  double precision, dimension(Sop%comp,Sop%comp) :: H
  double precision, dimension(Sop%comp,I%ngps,I%Nx,I%Tngpt) :: Q
  double precision, dimension(I%n,1) :: U
  double precision, dimension(I%ns,1) :: Uss
  double precision, dimension(Sop%comp,I%ngps,I%Nx,1) :: S,D,E
  double precision, dimension(I%Tngpt,1) :: Ta,T,dT
  double precision :: snorm,stagn
  double precision, dimension(:,:), allocatable :: Tn,Td ! numerator and denominator.
  double precision, dimension(:,:), allocatable :: T_up,dT_up ! Temporal functions of preliminary step.
  integer :: nbmi,nbme,dnbm

  ! Choice of the constant operator:

  H = inv(I%Hooke)  ! Constant operator.
  Q = G%d_Ep - L%d_Ep ! Difference between the global and local quantities.


  if (I%PGD_eval) then !(.true.) then

    ! PGD method to calculate the corrections functions:
    iter_pgd = 0; stagn = 100
    T = reshape(Top%vect_gpt, (/I%Tngpt,1/) ); dT = ones(I%Tngpt,1)
    Ta=T  ! Initial assignation to the temporal function.

    do while (iter_pgd<I%PGD_iter_max .and. stagn > I%PGD_ind)

      iter_pgd = iter_pgd + 1

      ! Space problem:
      call Space_problem(T,dT,H,Q,Sop,Top,I,E,D,S)        ![S,D,E,U] = Space_problem(T,dT,H,Q)

      ! Time problem:
      call Time_problem(S,D,E,H,Q,Sop,Top,I,T,dT)  ![T,dT] = Time_problem(S,D,E,H,Q)

      snorm = norm2(S)
      S = S/snorm; D = D/snorm; E = E/snorm; U = U/snorm
      T = T*snorm; dT = dT*snorm

      stagn = abs( norm2(T-Ta)/norm2(Ta) )
      !print*, "STAGNATION PGD", stagn
      Ta = T

    end do ! end while
    print*, "PGD mode determined."

    T = I%relax*T ; dT = I%relax*dT
    call Add_modes(S,D,E,U,T,dT,iter) ! Add the PGD modes (verificar)

    I%PGD_eval = .false.
    G%saved_const = .false.
  else ! Preliminary step:

    I%up_iter = I%up_iter + 1
    call Preliminary_step(H,Q,Sop,Top,I,nbmi,nbme,dnbm,T_up,dT_up)
    G%saved_const = .true. ! re-use the constants

    !call plot(dT_up(:,1))
    !stop

    ! Update indicator calculation:
    if (size(G%T_s,2)==1) then
      I%update_ind = (sum(0.5*abs(T_up)/( 1e-12 + abs(G%T_s) + abs(T_up)))/I%Tngpt)
      !call plot(dT_up(:,1))
      !stop
    else
      !I%update_ind = (maxval(reshape(sum(0.5*abs(T_up)/( 1e-12 + abs(G%T_s) + abs(T_up)),1)/I%Tngpt,(/size(G%T_s,2)/)),1))
      I%update_ind = maxval(sum(0.5*abs(T_up)/( 1e-12 + abs(G%T_s) + abs(T_up)),1)/I%Tngpt,1)
      !call plot(dT_up(:,1))
      !stop
    end if
    print*, I%update_ind

    if (I%update_ind < I%update_limit .or. I%up_iter > I%N_up_max) then
      I%PGD_eval = .true.
      I%up_iter = 0
      print*, "Preliminary step done."
    end if


    ! Add final relaxed corrections:
    G%T_s = G%T_s + I%relax*T_up
    G%dT_s = G%dT_s + I%relax*dT_up

  end if

  call assign_global(S,D,E,U,T,dT)  ! Assignation of the global quantities.


end subroutine global_stage_vp



subroutine Space_problem(T,dT,H,Q,Sop,Top,I,E,D,S) ! [S,D,E,U] output

  type(Input), intent(inout) :: I
  type(Space_op), intent(inout) :: Sop
  type(Time_op), intent(inout) :: Top

  double precision, dimension(I%n,I%n) :: M
  double precision, dimension(I%ns,I%ns) :: Mss
  double precision, dimension(I%n,1) :: F
  double precision, dimension(I%ns,1) :: Fss
  double precision :: term1, term2

  double precision, dimension(:,:,:,:) :: Q
  double precision, dimension(I%Tngps,1) :: delta
  double precision, dimension(Sop%comp,Sop%comp) :: H,W,invW

  double precision, dimension(I%Tngpt,1) :: T, dT ! Time modes.
  double precision, dimension(I%n,1) :: U ! space function of the dispacement correction.
  double precision, dimension(Sop%comp,I%ngps,I%Nx,1) :: E,D,S ! output spatial modes.
  double precision, dimension(:,:), allocatable :: E_aux, D_aux, S_aux
  double precision, dimension(:,:), allocatable :: aux1, aux2

  integer, dimension(4) :: siz
  integer, dimension(4) :: sizr
  integer :: ns, nrhs, lda, ldb, info
  double precision, dimension(:), allocatable :: ipiv, x
  double precision, dimension(:,:), allocatable :: plot_vec
  double precision, dimension(:,:), allocatable :: Integrant


  siz = shape(Q)
  term1 = sum( (T**2)*Top%IntT )
  term2 = sum( (T*dT)*Top%IntT )

  !print*, "NORM T", norm2(T)
  !print*, "NORM dT", norm2(dT)
  !print*, "NORM T INTEGRAL", norm2(T*Top%IntT)

  allocate(Integrant( size(Q)/I%Tngpt , I%Tngpt ))
  Integrant = spread( reshape( transpose(T*Top%IntT), (/I%Tngpt/) ) ,1,size(Q)/I%Tngpt )

  !print*, "SHAPE INTEGRANT:", shape(spread( reshape( transpose(T*Top%IntT), (/I%Tngpt/) ) ,1,size(Q)/I%Tngpt ))
  !print*, "NORMA PTM", norm2(sum( reshape(Q, (/size(Q)/I%Tngpt,I%Tngpt/) )*Integrant,2))

  invW = H*term1 + inv(Sop%He)*term2
  delta = reshape(sum( reshape(Q, (/size(Q)/I%Tngpt,I%Tngpt/) )*Integrant,2), (/size(delta),1/) )

  !call plot(delta(:,1))
  !stop


  ! Spatial Correction in displacement.
  U = 0 ! Initialization


  W = inv(invW)
  F = I%A*W*( matmul( transpose(Sop%Bw) , delta ) )
  Fss = F(Sop%dof_free,:)


  if (.false.) then

    M = I%A*W*( matmul( transpose(Sop%B) , Sop%Bw ) )
    Mss = M(Sop%dof_free,Sop%dof_free)


    !!!!!!!!!! Fast inverse:
    !ns, nrhs, lda, ldb, ipiv, x
    ns = size(Mss,1)
    nrhs = size(Fss,2)
    lda = ns; ldb = ns
    allocate(ipiv(ns)); allocate(x(ns))
    x = -Fss(:,1)
    call dgesv(ns,nrhs,Mss,lda,ipiv,x,ldb,info)
    U(Sop%dof_free,1) = x

  else

    allocate(x(ns))
    x = -Fss(:,1)
    call inv_chol(Elast%ChKss,x)
    U(Sop%dof_free,1) = (1/(term1 + term2))*x

  end if





  allocate( E_aux( size(delta),1 ) )
  allocate( D_aux( size(delta),1 ) )
  allocate( S_aux( size(delta),1 ) )


  ! Spatial Correction in deformation:
  E_aux = matmul(Sop%B,U)


  ! Spatial Correction in stress:
  allocate( aux1(Sop%comp,size(delta)/Sop%comp)  )
  aux1 = matmul( W, reshape( (E_aux + delta), (/Sop%comp, size(delta)/Sop%comp /) ))
  D_aux = reshape( aux1 , (/size(delta),1/) ) ! matmul( W, (E + delta) )  ! inv(invW)*(E + delta)


  ! Spatial Correction in plastic deformation:
  allocate( aux2(Sop%comp,size(delta)/Sop%comp)  )
  aux2 = matmul(H, reshape(D_aux, (/Sop%comp, size(delta)/Sop%comp/)))
  S_aux = ( reshape( aux2 , (/size(delta),1/) )*term1 - delta )/term2     ! ( (H*D)*term1 - delta )/term2


  ! Reshape of all the spatial fields:
  sizr = siz; sizr(4) = 1

  E = reshape(E_aux,sizr)
  D = reshape(D_aux,sizr)
  S = reshape(S_aux,sizr)


end subroutine Space_problem








subroutine Time_problem(S,D,E,H,Q,Sop,Top,I,T,dT)  !result(alfa,alfat) ! [alfa,alfat]

! S,D,E,H,Q,Sop,Top,I

  implicit none

  type(Input), intent(inout) :: I
  type(Space_op), intent(inout) :: Sop
  type(Time_op), intent(inout) :: Top

  double precision , dimension(:,:,:,:) :: Q
  double precision, dimension(:,:,:,:) :: S,D,E
  double precision, dimension(:,:,:,:), allocatable :: Sg,Dg,Eg,aux
  double precision, dimension(:,:), allocatable :: aux1
  double precision, dimension(:,:) :: H
  integer, dimension(4) :: siz
  integer, dimension(3) :: sizr
  double precision :: A11,A10,A00
  double precision, dimension(I%ngpt,1) :: IntTe
  double precision, dimension(I%ntc-1,I%ntc-1) :: Mt
  double precision, dimension(I%ntc,1) :: Fti
  double precision, dimension(I%ntc-1,1) :: Ft

  integer :: j
  integer :: ns, nrhs, lda, ldb, info
  double precision, dimension(:), allocatable :: ipiv, x
  double precision, dimension(I%ntc,1) :: alfa_fem
  integer , dimension(I%Nsft*I%NT) :: IDvtc
  double precision, dimension(I%Nsft*I%NT,1) :: alfac
  double precision, dimension(I%Nsft,I%NT) :: alfac2, alfacc
  double precision, dimension(I%ngpt,I%NT) :: alfae, alfate
  double precision, dimension(I%Tngpt,1) :: T, dT
  double precision, dimension(I%Tngpt,1) :: D1, D0

  !double precision, dimension(:,:), allocatable :: IntS
  !double precision, dimension(:,:,:,:), allocatable :: t1
  double precision, dimension(:,:), allocatable :: IntSg
  double precision, dimension(:,:), allocatable :: to,t1
  double precision, dimension(I%Nsft,I%NT) :: Felm, Felm1, Felm2
  double precision, dimension(:,:), allocatable :: oper

  ! multiplication quantities blas:
  double precision, dimension(I%Nsft,I%ngpt) :: a1,a2
  double precision, dimension(I%ngpt, I%NT) :: b1,b2



  siz = shape(Q); sizr = siz(1:3)

  ! Constants:
  allocate( aux1( size(D) , 1 ) )
  aux1 = reshape( matmul(H, reshape(D, (/Sop%comp, size(D)/Sop%comp /) ) ), (/size(D),1/) )
  A11 =  sum( (reshape(S, (/ size(S),1 /) )**2)*Sop%IntS )  !matmul( transpose(reshape(S, (/ size(S),1 /) )**2) ,  Sop%IntS )
  A10 = -sum( reshape(S, (/ size(S),1 /) )*aux1*Sop%IntS )  !matmul( transpose(reshape(S, (/ size(S),1 /) )*reshape(H*D,  (/ size(D),1 /)  )) , Sop%IntS )
  A00 =  sum( (aux1**2)*Sop%IntS )  !matmul( transpose(reshape(H*D,  (/ size(D),1 /) )**2) , Sop%IntS )


  allocate(Sg( siz(1),siz(2),siz(3),siz(4) ))
  allocate(Dg( siz(1),siz(2),siz(3),siz(4) ))
  allocate(Eg( siz(1),siz(2),siz(3),siz(4) ))

  Sg = spread( reshape(S, (/siz(1),siz(2),siz(3)/) ) ,4,I%Tngpt)  ! dim = 4, count = I%Tngpt
  Dg = spread( reshape(D, (/siz(1),siz(2),siz(3)/) ) ,4,I%Tngpt)
  Eg = spread( reshape(E, (/siz(1),siz(2),siz(3)/) ) ,4,I%Tngpt)

  !print*, "NEEEPEEE   2", shape(spread( reshape(S, (/siz(1),siz(2),siz(3)/) ) ,4,I%Tngpt))

  !print*, "HASTA AQUI" ,shape(Sop%IntS*reshape( (Sg*Q), (/size(Q)/I%Tngpt,I%Tngpt/)  ))

  allocate(IntSg( size(Sop%IntS,1), I%Tngpt  ))  ! operator used to integrate in space.

  !IntSg = spread(  reshape(Sop%IntS, (/size(Sop%IntS,1)/) ) , 2, I%Tngpt )
  IntSg = spread(  Sop%IntS(:,1) , 2, I%Tngpt )


  allocate(t1( 1 , I%Tngpt ))
  !t1 =  matmul( transpose(Sop%IntS), reshape( (Sg*Q), (/size(Q)/I%Tngpt,I%Tngpt/)  ) )
  D1(:,1) = -reshape(sum( IntSg*reshape( (Sg*Q), (/size(Q)/I%Tngpt,I%Tngpt/)),1), (/I%Tngpt/))
  !D1(:,1) = - reshape( t1, (/I%Tngpt/) )    !sum( t1 ,1 )

  !call plot(D1(:,1))
  !stop

  allocate( to( 1 , I%Tngpt) )
  allocate(aux( siz(1),siz(2),siz(3),siz(4) ))
  aux = reshape( matmul(H, reshape( Dg,  (/Sop%comp, size(Dg)/Sop%comp /) ) ) , (/siz(1),siz(2),siz(3),siz(4)/)  )
  !to = matmul( transpose(Sop%IntS) , reshape( aux*Q , (/size(Q)/I%Tngpt , I%Tngpt/) ) )
  D0(:,1) = reshape(sum( IntSg*reshape( aux*Q , (/size(Q)/I%Tngpt , I%Tngpt/)),1), (/I%Tngpt/))
  !D0(:,1) =  reshape( to, (/I%Tngpt/) ) ! sum(to,1)

  !call plot( D1(:,1) + D0(:,1)  )
  !stop

  Mt = (Top%I11Tr)*A11 + (Top%I01Tr)*A10 + (Top%I10Tr)*A10 + (Top%I00Tr)*A00

  if (.true.) then

    a1 = Top%vwalfat
    b1 = reshape(D1,  (/ I%ngpt, I%NT /) )
    a2 = Top%vwalfa
    b2 = reshape(D0,  (/ I%ngpt, I%NT /) )


    !CALL DGEMM (transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
    call dgemm('n','n',size(Felm1,1),size(Felm1,2),size(a1,2),1.D0,a1,size(a1,1),b1,size(b1,1),0.D0,Felm1,size(Felm1,1))
    call dgemm('n','n',size(Felm2,1),size(Felm2,2),size(a2,2),1.D0,a2,size(a2,1),b2,size(b2,1),0.D0,Felm2,size(Felm2,1))



  else

    Felm1 = matmul( Top%vwalfat , reshape(D1,  (/ I%ngpt, I%NT /) )  )   ! N shape functions x ngpt
    Felm2 = matmul( Top%vwalfa  , reshape(D0,  (/ I%ngpt, I%NT /) )  ) ! N shape functions x ngpt

    !stop

  end if

  Felm = Felm1 + Felm2  ! Elemental vectors.
  Fti = 0; Ft = 0
  do j=1,I%NT
      Fti(Top%IDvtc(:,j),1) = Fti(Top%IDvtc(:,j),1) + Felm(:,j)
  end do
  Ft = Fti(Top%dof_free,:)  ! Null initial condition.


  ! Obtention of the time functions over all the gauss points:

  !!!!!!!!!! Fast inverse matrix:
  !ns, nrhs, lda, ldb, ipiv, x
  ns = size(Mt,1)
  nrhs = size(Ft,2)
  lda = ns; ldb = ns
  allocate(ipiv(I%ntc-1)); allocate(x(I%ntc-1))
  x = Ft(:,1)
  call dgesv(ns,nrhs,Mt,lda,ipiv,x,ldb,info)
  alfa_fem(Top%dof_free,1) = x


  IDvtc = reshape(Top%IDvtc, (/ size(Top%IDvtc) /)  )
  alfac = reshape(alfa_fem(IDvtc,1), (/size(alfac),1/) )


  alfacc=reshape( alfac,  (/I%Nsft,I%NT/) )


  allocate( oper( I%ngpt,I%Nsft ) )
  oper = transpose(Top%valfa)
  T = reshape( matmul( oper , alfacc ) , (/I%Tngpt,1/) )  ! alfa in each gauss point

  !call plot(T(:,1))
  !stop

  !alfate =  matmul(  transpose(Top%Tmat_vectauxalfat), alfacc )
  dT =  reshape( matmul(  transpose(Top%valfat), alfacc ), (/I%Tngpt,1/) )  ! alfa in each gauss point

  !call plot(dT(:,1))
  !stop


end subroutine Time_problem





!!!!!!!!!!!!!!!!!!!!!!!!!!! Preliminary step !!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Preliminary_step(H,Q,Sop,Top,I,nbmi,nbme,dnbm,T,dT)  !result(alfa,alfat) ! [alfa,alfat]

! S,D,E,H,Q,Sop,Top,I


  implicit none

  type(Input), intent(inout) :: I
  type(Space_op), intent(inout) :: Sop
  type(Time_op), intent(inout) :: Top

  double precision, dimension(:,:) :: H  ! constant operator.
  double precision , dimension(:,:), allocatable :: Qr ! reshape of Q
  double precision , dimension(:,:,:,:) :: Q
  double precision, dimension(:,:,:,:), allocatable :: S,D  ! S: Plastic deformation tensor; D: Stress tensor.

  ! vector determination:
  double precision, dimension(:,:), allocatable :: aux ! auxiliar variable for vetors calculation.
  double precision, dimension(:,:), allocatable :: Sr
  double precision, dimension(:,:), allocatable :: D1, D0 !

  ! constants determination:
  double precision, dimension(:,:), allocatable :: aux1 ! auxiliar variable for constants calculation.
  double precision, dimension(:,:), allocatable :: A11,A10,A01,A00 ! constants factors.


  double precision, dimension(I%ngpt,1) :: IntTe
  !double precision, dimension(I%Nsft,I%NT) :: Felm1, Felm2  ! elemental vector
  double precision, dimension(:,:), allocatable :: Felm1, Felm2  ! elemental vector
  double precision, dimension(:,:,:), allocatable :: Felm


  double precision, dimension(:,:), allocatable :: Mt
  double precision, dimension(:,:), allocatable :: Ft,Fti ! assembly vector
  !double precision, dimension(I%ntc,1) :: Fti  ! assembly vector



  integer, dimension(4) :: siz
  integer, dimension(3) :: sizr
  integer :: Nelms
  integer :: nbmi, nbme, dnbm, nbm


  integer :: i1,j1,it,im
  integer :: ns, nrhs, lda, ldb, info
  double precision, dimension(:), allocatable :: ipiv, x
  integer , dimension(I%Nsft*I%NT) :: IDvtc


  double precision, dimension(:,:), allocatable :: alfa_fem
  double precision, dimension(:,:), allocatable :: alfac
  double precision, dimension(:,:), allocatable :: alfacc
  double precision, dimension(:,:), allocatable :: T, dT


  ! multiplication quantities blas:
  double precision, dimension(I%Nsft,I%ngpt) :: a1,a2
  double precision, dimension(:,:), allocatable :: b1,b2

  double precision, dimension(:), allocatable :: IntS





  ! Determination  of the modes that need to be actualized
  siz = shape(G%S_d_Ep); sizr = siz(1:3)
  nbm = siz(4)
  !print*, siz

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

  dnbm = nbme - nbmi + 1 ! number of modes that will be actualized.

  !print*, "dnbm:", dnbm
  !print*, "nbm:", siz(4)

  ! End Determination  of the modes that need to be actualized

  ! Assignation of the spatial functions from global_stage:

  allocate( S( siz(1),siz(2),siz(3),dnbm ) )
  allocate( D( siz(1),siz(2),siz(3),dnbm ) )




  S = G%S_d_Ep(:,:,:,nbmi:nbme)  ! Plastic deformation tensor.
  D = G%S_Sigma(:,:,:,nbmi:nbme) ! Stress tensor.

  ! end Assignation of the spatial functions from global_stage:



  ! Constants factors determination:

  allocate(A11(dnbm,dnbm)); allocate(A10(dnbm,dnbm))
  allocate(A01(dnbm,dnbm)); allocate(A00(dnbm,dnbm))

  Nelms = size(S(:,:,:,1)) ! Number of elements in space.

  !print*, Nelms
  !stop


  if (G%saved_const .eqv. .false.) then

    allocate( aux1( Nelms , dnbm ) )
    aux1 = reshape( matmul(H, reshape(D, (/Sop%comp, size(D)/Sop%comp /) ) ), (/Nelms,dnbm/) )


    do i1=1,dnbm
      do j1=1,dnbm
        A11(i1,j1) =  sum( (reshape(S(:,:,:,i1)*S(:,:,:,j1), (/ Nelms,1 /) ))*Sop%IntS)
        A10(i1,j1) = -sum( reshape(S(:,:,:,i1), (/ Nelms,1 /) )*reshape(aux1(:,j1),(/Nelms,1/))*Sop%IntS)
        A01(i1,j1) = -sum( reshape(aux1(:,i1),(/Nelms,1/))*reshape(S(:,:,:,j1), (/ Nelms,1 /) )*Sop%IntS)
        A00(i1,j1) =  sum( ( reshape(aux1(:,i1),(/Nelms,1/))*reshape(aux1(:,j1),(/Nelms,1/)) )*Sop%IntS)
      end do
    end do

    deallocate(G%A11); deallocate(G%A10);
    deallocate(G%A01); deallocate(G%A00);

    allocate(G%A11(dnbm,dnbm)); allocate(G%A10(dnbm,dnbm));
    allocate(G%A01(dnbm,dnbm)); allocate(G%A00(dnbm,dnbm));
    G%A11 = A11; G%A10 = A10; G%A01 = A01; G%A00 = A00

  else ! re-use the constants:
    !print*, "Reutilization of the constants"
    A11 = G%A11; A10 = G%A10; A01 = G%A01; A00 = G%A00
  end if

  !call show(A11)
  !call show(A10)
  !call show(A01)
  !call show(A00)

  !print*, "_______________________"

  ! Temporal functions of the rhs determination:

  siz = shape(Q);

  allocate(D1(I%Tngpt,dnbm)); allocate(D0(I%Tngpt,dnbm))
  D1 = 0.D0; D0 = 0.D0

  allocate(Qr( size(Q)/I%Tngpt , I%Tngpt ))
  allocate(Sr( size(Q)/I%Tngpt , dnbm ))
  allocate(aux( size(Q)/I%Tngpt , dnbm ))

  Sr = reshape( S, (/Nelms,dnbm/) )
  Qr = reshape( Q, (/size(Q)/I%Tngpt , I%Tngpt/) )
  aux = reshape( matmul(H, reshape( D,  (/Sop%comp, size(D)/Sop%comp /) ) ) , (/size(Q)/I%Tngpt , dnbm/)  )

  allocate( IntS( size(Sop%IntS) ) )
  IntS = reshape( Sop%IntS, (/size(Sop%IntS)/) )

  do it = 1,I%Tngpt
    do im = 1,dnbm
      D0(it,im) =  sum( IntS*aux(:,im)*Qr(:,it))
      D1(it,im) = -sum( IntS*Sr(:,im)*Qr(:,it))
    end do
  end do


  if (.false.) then ! (dnbm .eq. 1) then ! (.false.) then
    print*, "TEST PRELIMINAR:"
    call plot(D0(:,1)+D1(:,1))
    stop
  end if


  !print*, it
  !call show(D0)

  !call plot(D1(:,1))
  !call plot(D1(:,1))
  !stop





  ! DETERMINATION OF THE FEM MATRIX AND VECTOR:

  allocate(Mt( (I%ntc-1)*dnbm , (I%ntc-1)*dnbm ))
  Mt = kron(Top%I11Tr,A11) + kron(Top%I01Tr,A01) + kron(transpose(Top%I01Tr),A10) + kron(Top%I00Tr,A00)
  !Mt = (Top%I11Tr)*A11 + (Top%I01Tr)*A10 + (Top%I10Tr)*A10 + (Top%I00Tr)*A00

  allocate( Felm1( I%Nsft,I%NT*dnbm )  )
  allocate( Felm2( I%Nsft,I%NT*dnbm )  )

   Felm1 = 0.D0; Felm2 = 0.D0

  ! calculate the elementary finite vectors of rhs
  if (.false.) then

    a1 = Top%vwalfat
    a2 = Top%vwalfa

    allocate( b1( I%ngpt, I%NT*dnbm ) )
    allocate( b2( I%ngpt, I%NT*dnbm ) )

    b1 = reshape(D1,  (/ I%ngpt, I%NT*dnbm /) )
    b2 = reshape(D0,  (/ I%ngpt, I%NT*dnbm /) )


    !CALL DGEMM (transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
    call dgemm('n','n',size(Felm1,1),size(Felm1,2),size(a1,2),1.D0,a1,size(a1,1),b1,size(b1,1),0.D0,Felm1,size(Felm1,1))
    call dgemm('n','n',size(Felm2,1),size(Felm2,2),size(a2,2),1.D0,a2,size(a2,1),b2,size(b2,1),0.D0,Felm2,size(Felm2,1))

  else

    Felm1 = matmul( Top%vwalfat , reshape(D1,  (/ I%ngpt, I%NT*dnbm /) )  )   ! N shape functions x ngpt
    Felm2 = matmul( Top%vwalfa  , reshape(D0,  (/ I%ngpt, I%NT*dnbm /) )  ) ! N shape functions x ngpt

  end if


  allocate(Felm( I%Nsft,I%NT, dnbm ))
  !Felm =  reshape( Felm1 + Felm2,  (/ I%Nsft,I%NT,dnbm /) )


  allocate(Fti( I%ntc , dnbm ))
  allocate(Ft( (I%ntc-1)*dnbm , 1 )) ! rhs vector

  !print*, shape(Felm1)
  !print*, shape(Felm2)
  !print*, (/I%Nsft,I%NT, dnbm/)
  !print*, "HASTA ACA:"
  Felm = 0.D0
  Felm = reshape( Felm1 + Felm2, (/I%Nsft,I%NT,dnbm/) ) ! Elemental vectors.

  Fti = 0.D0; Ft = 0.D0
  do j1=1,I%NT
    do im=1,dnbm
      Fti(Top%IDvtc(:,j1),im) = Fti(Top%IDvtc(:,j1),im) + Felm(:,j1,im) !reshape(Felm(:,j1,im),(/I%Nsft/))
    end do
  end do
  Ft =  reshape( Fti(Top%dof_free,:), (/(I%ntc-1)*dnbm,1/) ) ! Null initial condition.


  if (.false.) then !(nbm .eq. 2) then !(.false.) then

    call show(Felm(:,I%NT-1,2))
    call show(Felm(:,I%NT,2))

    print*, "TEST PRELIMINAR Fti:"
    call plot(Fti(:,1))
    stop
  end if


  ! FEM RESOLUTION:
  allocate( alfa_fem( I%ntc, dnbm ) )
  alfa_fem = 0


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


  allocate(T(I%Tngpt,nbm))
  allocate(dT(I%Tngpt,nbm))

  T = 0; dT=0

  T(:,nbmi:nbme) = reshape( matmul( transpose(Top%valfa) , alfacc ) , (/I%Tngpt,dnbm/) )  ! alfa in each gauss point
  dT(:,nbmi:nbme) = reshape( matmul(  transpose(Top%valfat), alfacc ), (/I%Tngpt,dnbm/) )  ! alfa in each gauss point


  if (.false.) then ! (dnbm .eq. 2) then
    print*, "TEST PRELIMINAR Ft:"
    call plot(dT(:,2))
    stop
  end if

  !stop
  !call plot(T(:,1))
  !stop

  !call plot(dT(:,1))
  !stop



end subroutine Preliminary_step




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










! Reviewed:
subroutine Add_modes(S,D,E,U,T,dT,iter) ![global_stage]=Add_modes(global_stage,S,D,E,U,T,dT,iter)

  implicit none

  integer :: nmodes, iter
  integer, dimension(4) :: siz


  double precision, dimension(:,:,:,:) :: S, D, E
  double precision, dimension(:,:) :: U
  double precision, dimension(:,:) :: T, dT
  ! If we have already some modes, variables used to cat the arrays on G%
  double precision, dimension(:,:,:,:), allocatable :: G_Ep, G_Sigma, G_E
  double precision, dimension(:,:), allocatable :: G_U, G_T_s, G_dT_s

  !integer :: ns



  !ns = size(U,1)

  !print*, "Modes Added"

  if (iter==1) then

      siz = shape(S)

      ! PGD modes:
      !double precision, dimension(:,:), allocatable :: U, T_s, dT_s ! Space and Temporal modes
      !double precision, dimension(:,:,:,:), allocatable :: S_eps, S_Sigma, S_d_Ep ! Spatial modes.


      allocate( G%S_d_Ep( siz(1), siz(2), siz(3), 1 ) )
      allocate( G%S_Sigma( siz(1), siz(2), siz(3), 1 ) )
      allocate( G%S_eps( siz(1), siz(2), siz(3), 1 ) )

      allocate( G%U( size(U,1),1) )
      allocate( G%T_s( size(T,1),1) )
      allocate( G%dT_s( size(T,1),1) )



      G%S_d_Ep = S
      G%S_Sigma = D
      G%S_eps = E
      G%U = U
      G%T_s  = T
      G%dT_s = dT

  else

      siz = shape(G%S_d_Ep)

      ! Save last arrays:
      allocate( G_Ep( siz(1), siz(2), siz(3), siz(4) ) )
      allocate( G_Sigma( siz(1), siz(2), siz(3), siz(4) ) )
      allocate( G_E( siz(1), siz(2), siz(3), siz(4) ) )
      allocate( G_U( size(U,1), siz(4) ) )

      allocate( G_T_s( size(T,1), siz(4) ) )
      allocate( G_dT_s( size(dT,1), siz(4) ) )

      G_Ep = G%S_d_Ep
      G_Sigma = G%S_Sigma
      G_E = G%S_eps
      G_U = G%U

      G_T_s = G%T_s
      G_dT_s = G%dT_s

      ! end Save last arrays.



      ! New sized arrays:
      nmodes = siz(4) + 1   ! New quantity of modes.
      siz(4) = nmodes

      deallocate(G%S_d_Ep); deallocate(G%S_Sigma);
      deallocate(G%S_eps); deallocate(G%U);
      deallocate(G%T_s); deallocate(G%dT_s);


      allocate( G%S_d_Ep( siz(1), siz(2), siz(3), siz(4) ) )
      allocate( G%S_Sigma( siz(1), siz(2), siz(3), siz(4) ) )
      allocate( G%S_eps( siz(1), siz(2), siz(3), siz(4) ) )
      allocate( G%U( size(U,1), siz(4) ) )

      allocate( G%T_s( size(T,1), siz(4) ) )
      allocate( G%dT_s( size(dT,1), siz(4) ) )

      ! Adjoint of already computed modes:

      G%S_d_Ep(:,:,:,1:nmodes-1) = G_Ep !cat(4,G%Ep,S)
      G%S_Sigma(:,:,:,1:nmodes-1) = G_Sigma !cat(4,G%Sigma,D)
      G%S_eps(:,:,:,1:nmodes-1) = G_E  ! cat(4,G%E,E)
      G%U(:,1:(nmodes-1)) = G_U !horzcat(G%U,U)

      G%T_s(:,1:(nmodes-1))  = G_T_s !horzcat(G%T_s,T)
      G%dT_s(:,1:(nmodes-1)) = G_dT_s !horzcat(G%dT_s,dT)

      ! Concatenation of the new mode:

      G%S_d_Ep(:,:,:,nmodes) = S(:,:,:,1) !cat(4,G%Ep,S)
      G%S_Sigma(:,:,:,nmodes) = D(:,:,:,1) !cat(4,G%Sigma,D)
      G%S_eps(:,:,:,nmodes) = E(:,:,:,1)  ! cat(4,G%E,E)
      G%U(:,nmodes) = U(:,1) !horzcat(G%U,U)

      G%T_s(:,nmodes)  = T(:,1) !horzcat(G%T_s,T)
      G%dT_s(:,nmodes) = dT(:,1) !horzcat(G%dT_s,dT)

  end if


end subroutine Add_modes ! end Add_modes







subroutine assign_global(S,D,E,U,T,dT)

  implicit none

  double precision, dimension(:,:) :: U,T,dT
  double precision, dimension(:,:,:,:) :: S,D,E

  !call plot(T(:,1))

  if (.false.) then

    G%Utot = G%Utot + full_field(U,T)
    G%eps =  G%eps + full_field(E,T)
    G%sigma = G%sigma + full_field(D,T)
    G%d_Ep = G%d_Ep + full_field(S,dT)

  else

    ! PGD modes:
    !double precision, dimension(:,:), allocatable :: U, T_s, dT_s ! Space and Temporal modes
    !double precision, dimension(:,:,:,:), allocatable :: S_eps, S_Sigma, S_d_Ep ! Spatial modes.


    G%Utot = G%Utot_o + full_field(G%U,G%T_s)
    G%eps =  G%eps_o + full_field(G%S_eps,G%T_s)
    G%sigma = G%sigma_o + full_field(G%S_Sigma,G%T_s)
    G%d_Ep = full_field(G%S_d_Ep,G%dT_s)

  end if


end subroutine assign_global





end module global_stage
