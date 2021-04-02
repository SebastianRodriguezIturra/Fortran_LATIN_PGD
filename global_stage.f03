
module global_stage

  use Structure
  use Matrix_op
  use local_stage
  use Elastic_sol
  use bar
  use Time_matrices

  implicit none

  type :: global

    double precision, dimension(:,:), allocatable :: Utot ! Displacement.
    double precision, dimension(:,:,:,:), allocatable :: eps ! Strain.
    double precision, dimension(:,:,:,:), allocatable :: sigma ! Stress.
    double precision, dimension(:,:,:,:), allocatable :: d_Ep ! Plastic deformation.

    ! PGD modes:
    double precision, dimension(:,:), allocatable :: U, T_s, dT_s ! Space and Temporal modes
    double precision, dimension(:,:,:,:), allocatable :: S_eps, S_Sigma, S_d_Ep ! Spatial modes.


  end type global

  type(global) :: G


  contains


subroutine assign_elastic(Elast,I)
  implicit none
  type(Init_sol), intent(inout) :: Elast
  type(Input), intent(inout) :: I

  allocate(G%Utot( I%n, I%Tngpt ))
  allocate(G%eps( Sop%comp,I%ngps,I%Nx, I%Tngpt ))
  allocate(G%sigma( Sop%comp,I%ngps,I%Nx, I%Tngpt ))
  allocate(G%d_Ep( Sop%comp,I%ngps,I%Nx, I%Tngpt ))



  G%Utot = Elast%U
  G%eps = Elast%eps
  G%sigma = Elast%sigma
  G%d_Ep = 0

  !print*, "BLAAAAA" ,shape(G%d_Ep)

end subroutine assign_elastic





subroutine global_stage_vp(Top,Sop,L,I,iter)  ! [global_stage]=Global_stage_plasticity(global_stage,local,I,iter)

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


  ! Choice of the constant operator:

  H = inv(I%Hooke)  ! Constant operator.
  Q = G%d_Ep - L%d_Ep ! Difference between the global and local quantities.

  ! PGD method to calculate the corrections functions:

  iter_pgd = 0; stagn = 100
  T = reshape(Top%vect_gpt, (/I%Tngpt,1/) ); dT = ones(I%Tngpt,1)

  Ta=T  ! Initial assignation to the temporal function.

  do while (iter_pgd<I%PGD_iter_max .and. stagn > I%PGD_ind)

    iter_pgd = iter_pgd + 1

    !print*, "PGD ITERS" , iter_pgd

    ! Space problem:
    call Space_problem(T,dT,H,Q,Sop,Top,I,E,D,S)        ![S,D,E,U] = Space_problem(T,dT,H,Q)

    !print*, "NORM S", norm2(S)
    ! Time problem:
    call Time_problem(S,D,E,H,Q,Sop,Top,I,T,dT)  ![T,dT] = Time_problem(S,D,E,H,Q)


    snorm = norm2(S)
    !print*, "NORM S", norm2(S)

    !stop

    S = S/snorm; D = D/snorm; E = E/snorm; U = U/snorm
    T = T*snorm; dT = dT*snorm



    stagn = abs( norm2(T-Ta)/norm2(Ta) )
    !print*, "STAGNATION PGD", stagn
    Ta = T

  end do ! end while

  T = I%relax*T ; dT = I%relax*dT


  ![global_stage] = Add_modes(global_stage,S,D,E,U,I.relax*T,I.relax*dT,iter)
  !call Add_modes(S,D,E,U,T,dT,iter) ! Add the PGD modes (verificar)
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

  W = inv(invW)
  !print*, W
  M = I%A*W*( matmul( transpose(Sop%B) , Sop%Bw ) )
  F = I%A*W*( matmul( transpose(Sop%Bw) , delta ) )

  !call plot(F(:,1))
  !stop

  Mss = M(Sop%dof_free,Sop%dof_free)
  Fss = F(Sop%dof_free,:)

  !call show(M)
  !stop

  !call show(F)
  !stop
  !call plot(F(:,1))
  !stop
  ! Spatial Correction in displacement.
  U = 0 ! Initialization


  !!!!!!!!!! Fast inverse:
  !ns, nrhs, lda, ldb, ipiv, x
  ns = size(Mss,1)
  nrhs = size(Fss,2)
  lda = ns; ldb = ns
  allocate(ipiv(ns)); allocate(x(ns))
  x = -Fss(:,1)
  !print*, 1e10*sum(-Fss,1), 1e10*sum(x)
  !print*, max(-Fss), max(x)
  call dgesv(ns,nrhs,Mss,lda,ipiv,x,ldb,info)
  U(Sop%dof_free,1) = x

  !call show(x)
  !call plot(x)
  !stop

  !print*, shape(U)
  !!!!!!!!!!

  allocate( E_aux( size(delta),1 ) )
  allocate( D_aux( size(delta),1 ) )
  allocate( S_aux( size(delta),1 ) )

  ! Spatial Correction in deformation:

  !print*, "NORM W", norm2(W)
  !print*, "NORM DELTA", norm2(delta)
  !print*, "NORM Q", norm2(Q)

  E_aux = matmul(Sop%B,U)

  !print*, "NORM E", norm2(E_aux)
  !stop

  ! Spatial Correction in stress:
  allocate( aux1(Sop%comp,size(delta)/Sop%comp)  )
  aux1 = matmul( W, reshape( (E_aux + delta), (/Sop%comp, size(delta)/Sop%comp /) ))
  D_aux = reshape( aux1 , (/size(delta),1/) ) ! matmul( W, (E + delta) )  ! inv(invW)*(E + delta)

  !print*, "NORM D", norm2(D_aux)

  ! Spatial Correction in plastic deformation:
  allocate( aux2(Sop%comp,size(delta)/Sop%comp)  )
  aux2 = matmul(H, reshape(D_aux, (/Sop%comp, size(delta)/Sop%comp/)))
  S_aux = ( reshape( aux2 , (/size(delta),1/) )*term1 - delta )/term2     ! ( (H*D)*term1 - delta )/term2

  !print*, "NORM TEST:", norm2(U)
  !print*, "NORM TEST:", norm2(E_aux)
  !print*, "NORM TEST:", (norm2(D_aux))
  !print*, "NORM S:", (norm2(S_aux))

  !stop

  ! Reshape of all the spatial fields:
  sizr = siz; sizr(4) = 1

  E = reshape(E_aux,sizr)
  D = reshape(D_aux,sizr)
  S = reshape(S_aux,sizr)

  !print*, "NORM TEST:", (norm2(S))

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

  !print*, norm2(S), norm2(D), norm2(E)
  !print*, A11, A10, A00
  !stop
  !print*, "NEEEPEEE   1"

  !print*, norm2(Q(:,:,:,200))
  !print*, shape(reshape( Q, (/size(Q)/I%Tngpt,I%Tngpt/)  ))



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


  !D0 =  reshape(sum(reshape(IntS.*(H*D).*Q,[],I.Ngpt),1),[],1);
  allocate( to( 1 , I%Tngpt) )
  allocate(aux( siz(1),siz(2),siz(3),siz(4) ))
  aux = reshape( matmul(H, reshape( Dg,  (/Sop%comp, size(Dg)/Sop%comp /) ) ) , (/siz(1),siz(2),siz(3),siz(4)/)  )
  !to = matmul( transpose(Sop%IntS) , reshape( aux*Q , (/size(Q)/I%Tngpt , I%Tngpt/) ) )
  D0(:,1) = reshape(sum( IntSg*reshape( aux*Q , (/size(Q)/I%Tngpt , I%Tngpt/)),1), (/I%Tngpt/))
  !D0(:,1) =  reshape( to, (/I%Tngpt/) ) ! sum(to,1)

  !call plot( D1(:,1) + D0(:,1)  )
  !stop

  !print*, shape(D0)
  !call plot(D0(:,1))
  !stop
  ! Time problem resolution:

  !I00Tr = Top%I00Tr
  !I01Tr = Top%I01Tr
  !I11Tr = Top%I11Tr

  Mt = (Top%I11Tr)*A11 + (Top%I01Tr)*A10 + (Top%I10Tr)*A10 + (Top%I00Tr)*A00


  ! RHS:
  !IntTe = reshape( Top%IntT(1:I%ngpt,1) ,  (/ I%ngpt,1 /) )
  !Felm1 = matmul( (IntTe*(Top%Tmat_vectauxalfat)) , reshape(D1,  (/ I%ngpt, I%NT /) )  )   ! N shape functions x ngpt
  !Felm2 = matmul( (IntTe*(Top%Tmat_vectauxalfa))  , reshape(D0,  (/ I%ngpt, I%NT /) )  ) ! N shape functions x ngpt




  if (.true.) then

    a1 = Top%vwalfat
    b1 = reshape(D1,  (/ I%ngpt, I%NT /) )
    a2 = Top%vwalfa
    b2 = reshape(D0,  (/ I%ngpt, I%NT /) )


    !CALL DGEMM (transa, transb, l, n, m, alpha, a, lda, b, ldb, beta, c, ldc)
    call dgemm('n','n',size(Felm1,1),size(Felm1,2),size(a1,2),1.D0,a1,size(a1,1),b1,size(b1,1),0.D0,Felm1,size(Felm1,1))
    call dgemm('n','n',size(Felm2,1),size(Felm2,2),size(a2,2),1.D0,a2,size(a2,1),b2,size(b2,1),0.D0,Felm2,size(Felm2,1))



  else


    !print*, "ACA"

    Felm1 = matmul( Top%vwalfat , reshape(D1,  (/ I%ngpt, I%NT /) )  )   ! N shape functions x ngpt
    Felm2 = matmul( Top%vwalfa  , reshape(D0,  (/ I%ngpt, I%NT /) )  ) ! N shape functions x ngpt

    !stop

  end if

  !call show(Top%vwalfat)
  !call show(Top%vwalfa)
  !call show(Felm2)

  !print*, norm2(Felm1), norm2(Felm2)
  !stop


  Felm = Felm1 + Felm2  ! Elemental vectors.

  !print*, norm2(Felm1), norm2(Felm2), norm2(Felm)
  !print*, maxval(Felm1), maxval(Felm2), maxval(Felm)
  !stop


  !call show(Top%IDvtc)
  !stop

  !Ft = zeros(I%ntc,1)  ! Total vector.
  Fti = 0; Ft = 0
  do j=1,I%NT
      Fti(Top%IDvtc(:,j),1) = Fti(Top%IDvtc(:,j),1) + Felm(:,j)
  end do
  Ft = Fti(Top%dof_free,:)  ! Null initial condition.

  !print*, Fti
  !call plot(Fti(:,1))
  !stop
  !
  !print*, "QUE PASA?"
  !print*, Top%IDvtc(:,1), Top%IDvtc(:,2), Top%IDvtc(:,I%NT)
  !stop

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


  !call plot(x)
  !call plot(alfa_fem(:,1))
  !stop

  ! CORRECTO HASTA ACA!
  !!!!!!!!!!



  !alfa_fem=(Mt\Ft)
  !alfa_fem = vertcat(0,alfa_fem)
  IDvtc = reshape(Top%IDvtc, (/ size(Top%IDvtc) /)  )

  !print*, IDvtc
  !stop


  alfac = reshape(alfa_fem(IDvtc,1), (/size(alfac),1/) )


  !call plot(reshape(alfac, (/size(alfac)/) ))
  !call show(alfac)
  !call plot(alfac(:,1))
  !stop

  alfacc=reshape( alfac,  (/I%Nsft,I%NT/) )

  !print*, "STRANGE VALUE max"
  !print*, maxval(alfacc)
  !call plot(reshape(alfacc, (/size(alfac)/) ))


  !stop
  !print*, maxval( reshape(alfac2,  (/ size(alfac2) /)) )
  !print*, matmul( transpose(Top%Tmat_vectauxalfa), reshape(alfac2(:,I%NT), (/I%Nsft,1/)))
  !stop
  !call plot( reshape(alfac2,  (/ size(alfac2) /)  ) )
  !stop
  ! Correcto

  allocate( oper( I%ngpt,I%Nsft ) )
  oper = transpose(Top%valfa)
  !call show(oper)
  !print*, "NEPE"
  !print*, reshape(alfac2, (/size(alfac2)/) )
  !call plot(reshape(alfacc, (/size(alfacc)/) ))
  !stop
  !alfae = matmul( oper , alfacc )

  !print*, "STRANGE VALUE max"
  !print*, maxval( matmul( oper , alfacc ) )
  !stop

  !alfai = matmul( transpose(Top%Tmat_vectauxalfa) , reshape( alfac,  (/I%Nsft,I%NT/) ) )


  !print*, "SIZE:", shape( matmul( transpose(Top%Tmat_vectauxalfa) , alfac2 ) )
  !stop

  !call plot( reshape( matmul( transpose(Top%Tmat_vectauxalfa) , alfac2 ), (/I%Tngpt/)  )  )
  !stop

  T = reshape( matmul( oper , alfacc ) , (/I%Tngpt,1/) )  ! alfa in each gauss point

  !call plot(T(:,1))
  !stop

  !alfate =  matmul(  transpose(Top%Tmat_vectauxalfat), alfacc )
  dT =  reshape( matmul(  transpose(Top%valfat), alfacc ), (/I%Tngpt,1/) )  ! alfa in each gauss point

  !call plot(dT(:,1))
  !stop


end subroutine Time_problem





!!!!!!!!!!!!!!!!!!!!!!!!!!! Preliminary step !!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Preliminary_step(S,D,E,H,Q,Sop,Top,I,T,dT)  !result(alfa,alfat) ! [alfa,alfat]

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


  allocate(IntSg( size(Sop%IntS,1), I%Tngpt  ))  ! operator used to integrate in space.

  !IntSg = spread(  reshape(Sop%IntS, (/size(Sop%IntS,1)/) ) , 2, I%Tngpt )
  IntSg = spread(  Sop%IntS(:,1) , 2, I%Tngpt )


  allocate(t1( 1 , I%Tngpt ))
  !t1 =  matmul( transpose(Sop%IntS), reshape( (Sg*Q), (/size(Q)/I%Tngpt,I%Tngpt/)  ) )
  D1(:,1) = -reshape(sum( IntSg*reshape( (Sg*Q), (/size(Q)/I%Tngpt,I%Tngpt/)),1), (/I%Tngpt/))
  !D1(:,1) = - reshape( t1, (/I%Tngpt/) )    !sum( t1 ,1 )

  !call plot(D1(:,1))
  !stop


  !D0 =  reshape(sum(reshape(IntS.*(H*D).*Q,[],I.Ngpt),1),[],1);
  allocate( to( 1 , I%Tngpt) )
  allocate(aux( siz(1),siz(2),siz(3),siz(4) ))
  aux = reshape( matmul(H, reshape( Dg,  (/Sop%comp, size(Dg)/Sop%comp /) ) ) , (/siz(1),siz(2),siz(3),siz(4)/)  )
  !to = matmul( transpose(Sop%IntS) , reshape( aux*Q , (/size(Q)/I%Tngpt , I%Tngpt/) ) )
  D0(:,1) = reshape(sum( IntSg*reshape( aux*Q , (/size(Q)/I%Tngpt , I%Tngpt/)),1), (/I%Tngpt/))

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

  end if



  Felm = Felm1 + Felm2  ! Elemental vectors.


  Fti = 0; Ft = 0
  do j=1,I%NT
      Fti(Top%IDvtc(:,j),1) = Fti(Top%IDvtc(:,j),1) + Felm(:,j)
  end do
  Ft = Fti(Top%dof_free,:)  ! Null initial condition.


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


  !call plot(reshape(alfac, (/size(alfac)/) ))
  !call show(alfac)
  !call plot(alfac(:,1))
  !stop

  alfacc=reshape( alfac,  (/I%Nsft,I%NT/) )


  allocate( oper( I%ngpt,I%Nsft ) )
  oper = transpose(Top%valfa)



  T = reshape( matmul( oper , alfacc ) , (/I%Tngpt,1/) )  ! alfa in each gauss point
  dT =  reshape( matmul(  transpose(Top%valfat), alfacc ), (/I%Tngpt,1/) )  ! alfa in each gauss point

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


  if (iter==1) then



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


      ! New sized arrays:
      nmodes = siz(4) + 1   ! New quantity of modes.
      siz(4) = nmodes

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

  G%Utot = G%Utot + full_field(U,T)
  G%eps =  G%eps + full_field(E,T)
  G%sigma = G%sigma + full_field(D,T)
  G%d_Ep = G%d_Ep + full_field(S,dT)


end subroutine assign_global





end module global_stage
