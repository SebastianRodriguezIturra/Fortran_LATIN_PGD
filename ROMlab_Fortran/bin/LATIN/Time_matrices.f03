module Time_matrices


    use Structure
    use Matrix_op
    use Gauss

    implicit none


    type :: Time_op

      integer, dimension(:), allocatable :: dof_imposed
      integer, dimension(:), allocatable :: dof_free
      double precision, dimension(:,:), allocatable :: gauss_t ! points and weight matrix.
      double precision, dimension(:,:), allocatable :: IntT,vect_gpt ! Displacement.
      double precision, dimension(:,:), allocatable :: Q00k,Q01k,Q10k,Q11k ! Strain.
      double precision, dimension(:,:), allocatable :: I00Tr,I01Tr,I10Tr,I11Tr ! Stress.
      integer, dimension(:,:), allocatable :: IDvtc ! Plastic deformation.
      double precision, dimension(:,:), allocatable :: l,lt,lw,lwt
      double precision, dimension(:,:), allocatable :: Tmat_l, Tmat_lt, Tmat_timefunc,Tmat_timefunct
      double precision, dimension(:,:), allocatable :: valfa,valfat
      double precision, dimension(:,:), allocatable :: vwalfa,vwalfat

    end type Time_op

    type(Time_op) :: Top

    contains





    subroutine Time_op_assign(I)

      type(Input), intent(inout) :: I
      double precision, dimension(I%ngpt,1) :: waux
      double precision, dimension(I%ngpt) :: t
      double precision, dimension(I%ngpt,I%Nsft) :: IntTe
      double precision, dimension(I%Nsft,I%Nsft) :: Q00k, Q01k, Q10k, Q11k
      double precision, dimension(I%ntc-1,I%ntc-1) :: I00Tr,I01Tr,I10Tr,I11Tr
      double precision, dimension(I%ntc,I%ntc) :: I00T, I01T, I10T, I11T
      integer :: igp,iNt,ind,k,igpt
      double precision :: a,b
      integer :: elmt


      ! Vector that integrate any function in time over I=[0,Tend]
      allocate(Top%IntT(I%Tngpt,1))
      Top%gauss_t = gauss_points_weights(I%ngpt)
      do igp = 1,I%ngpt
        waux(igp,1)=(I%dt/2)*Top%gauss_t(igp,2)
      end do
      Top%IntT=kron(ones(I%NT,1),waux)


      allocate(Top%vect_gpt(1,I%Tngpt))
      a=-I%dt; b=0; ind=0
      do iNt = 1,I%NT
        a=a+I%dt;  b=b+I%dt
        do igpt = 1,I%ngpt
            ind=ind+1
            Top%vect_gpt(1,ind)=((b-a)/2)*Top%gauss_t(igpt,1) + ((b+a)/2)
        end do ! igpt
      end do! iNt


      ! Gauss points for a normalized time element:

      t = ((I%dt/2)*Top%gauss_t(:,1)+I%dt/2)

      allocate( Top%IDvtc( I%Nsft, I%NT ) )

      if ( I%Tmethod .eq. "Linear" ) then

        ! Identification matrix in time (for linear shape functions):
        Top%IDvtc = vertcat( incr(1,I%NT), incr(2,I%NT+1) )


        !                       Linear shape functions in time:

        allocate(Top%l(I%Nsft,I%ngpt))
        allocate(Top%lt(I%Nsft,I%ngpt))
        Top%l = vertcat(1-t/I%dt,t/I%dt) ! Polinomias that multiply each value of the node of each time element
        Top%lt = repmat( reshape( (/-1/I%dt, 1/I%dt/), (/2,1/) ), (/1,I%ngpt/) )

        ! Functions used to calculate the classic PGD mode:
        Top%Tmat_l = repmat(Top%l, (/1,I%NT/) )
        Top%Tmat_lt = repmat( Top%lt, (/1,I%NT/) )

        Top%valfa = Top%l
        Top%valfat = Top%lt


      else if ( I%Tmethod .eq. "Lagrange" )  then  ! Lagrange Polinomias


        do elmt=1,I%NT
          Top%IDvtc(:,elmt)= incr(2*(elmt-1)+1,2*elmt+1)
        end do


        !call show(Top%IDvtc(:,1:3))
        !stop

        allocate(Top%l(I%Nsft,I%ngpt))
        allocate(Top%lt(I%Nsft,I%ngpt))



        Top%l(1,:) = (2*t/I%dt-1)*(t/I%dt-1)
        Top%l(2,:) = 4*(t/I%dt)*(1-t/I%dt)
        Top%l(3,:) = (t/I%dt)*(2*t/I%dt-1)


        Top%lt(1,:) = ((2/I%dt)*(t/I%dt-1)+(2*t/I%dt -1)*(1/I%dt))
        Top%lt(2,:) = ((4/I%dt)*(1-t/I%dt)+(4*t/I%dt)*(-1/I%dt))
        Top%lt(3,:) = ((1/I%dt)*(2*t/I%dt -1) + (t/I%dt)*(2/I%dt))



        ! Functions used to calculate the classic PGD mode:
        Top%Tmat_l = repmat(Top%l, (/1,I%NT/) )
        Top%Tmat_lt = repmat( Top%lt, (/1,I%NT/) )

        Top%valfa = Top%l
        Top%valfat = Top%lt


      else ! Hermite


        !term1=1:2:1+2*(I.NT-1);
        !term2=2:2:2+2*(I.NT-1);
        !I.IDvtc=vertcat(term1,term2,term1+2,term2+2);



        !l=[ ( 1-(3/(dt^2))*t.^2 + (2/(dt^3))*t.^3 ); (t-(2/dt)*t.^2 + (1/(dt^2))*t.^3 ) ; ((3/(dt^2))*t.^2 + -(2/(dt^3))*t.^3 ) ; (-(1/dt)*t.^2 + (1/(dt^2))*t.^3)]; % gives the value at one point of all the shape functions.
        !lt=[(-2*(3/(dt^2))*t + 3*(2/(dt^3))*t.^2 ); (1-2*(2/dt)*t+3*(1/(dt^2))*t.^2) ; (2*(3/(dt^2))*t+-3*(2/(dt^3))*t.^2) ; (-2*(1/dt)*t+3*(1/(dt^2))*t.^2)];
        !ltt=[(-2*(3/(dt^2)) + 6*(2/(dt^3))*t) ; (-2*(2/dt) + 6*(1/(dt^2))*t) ; (2*(3/(dt^2)) + -6*(2/(dt^3))*t) ; (-2*(1/dt) + 6*(1/(dt^2))*t)];


        !l0=[( 1-(3/(dt^2))*ti.^2 + (2/(dt^3))*ti.^3) , (ti -(2/dt)*ti.^2 + (1/(dt^2))*ti.^3) , ((3/(dt^2))*ti.^2 + -(2/(dt^3))*ti.^3) , (-(1/dt)*ti.^2 + (1/(dt^2))*ti.^3)];
        !ldt=[( 1 -(3/(dt^2))*tf.^2 + (2/(dt^3))*tf.^3) , (tf -(2/dt)*tf.^2 + (1/(dt^2))*tf.^3) , ((3/(dt^2))*tf.^2 + -(2/(dt^3))*tf.^3) , (-(1/dt)*tf.^2 + (1/(dt^2))*tf.^3)];
        !lt0=[(-2*(3/(dt^2))*ti + 3*(2/(dt^3))*ti.^2) , (1 -2*(2/dt)*ti + 3*(1/(dt^2))*ti.^2) , (2*(3/(dt^2))*ti + -3*(2/(dt^3))*ti.^2) , (-2*(1/dt)*ti + 3*(1/(dt^2))*ti.^2)];
        !ltdt=[(-2*(3/(dt^2))*tf + 3*(2/(dt^3))*tf.^2) , (1 -2*(2/dt)*tf + 3*(1/(dt^2))*tf.^2) , (2*(3/(dt^2))*tf + -3*(2/(dt^3))*tf.^2) , (-2*(1/dt)*tf + 3*(1/(dt^2))*tf.^2)];

      end if

      !call show(Top%l)
      !print*, "espacio"
      !call show(Top%lt)
      !stop


      ! Calculation of the elementals matrices:
      call Elm_calc_matrices(Top%l,Top%lt,waux,Q00k,Q01k,Q10k,Q11k)


      ! %%%%%%%%%%%%%%%%  Save the differents matrix %%%%%%%%%%%%%%%%%

      Top%Q00k = Q00k
      Top%Q01k = Q01k
      Top%Q10k = Q10k
      Top%Q11k = Q11k



      !call show(Top%Q00k)
      !call show(Top%Q01k)
      !print*, "separacion:"
      !call show(Top%Q10k)
      !print*, "separacion:"
      !call show(Top%Q11k)

      !print*, "STOP AQUI"
      !stop


      ! Global matrix in time (Quasi_static calculation):


      ! % Definition of the global integral matrices of the LATIN - PGD integration steps.

      ! %% Ensamble de matrices en la matriz global.
      do k = 1,I%NT
          I00T(Top%IDvtc(:,k),Top%IDvtc(:,k)) = I00T(Top%IDvtc(:,k),Top%IDvtc(:,k)) + Q00k
          I01T(Top%IDvtc(:,k),Top%IDvtc(:,k)) = I01T(Top%IDvtc(:,k),Top%IDvtc(:,k)) + Q01k
          I10T(Top%IDvtc(:,k),Top%IDvtc(:,k)) = I10T(Top%IDvtc(:,k),Top%IDvtc(:,k)) + Q10k
          I11T(Top%IDvtc(:,k),Top%IDvtc(:,k)) = I11T(Top%IDvtc(:,k),Top%IDvtc(:,k)) + Q11k
      end do

      ! Free dofs:
      allocate( Top%dof_imposed(1)  )
      allocate( Top%dof_free(I%ntc-1)  )

      Top%dof_imposed = (/1/)
      Top%dof_free = incr(2,I%ntc)

      I00Tr = I00T(Top%dof_free,Top%dof_free)
      I01Tr = I01T(Top%dof_free,Top%dof_free)
      I10Tr = I10T(Top%dof_free,Top%dof_free)
      I11Tr = I11T(Top%dof_free,Top%dof_free)

      allocate(Top%I00Tr(I%ntc-1,I%ntc-1))
      allocate(Top%I01Tr(I%ntc-1,I%ntc-1))
      allocate(Top%I10Tr(I%ntc-1,I%ntc-1))
      allocate(Top%I11Tr(I%ntc-1,I%ntc-1))

      Top%I00Tr = I00Tr
      Top%I01Tr = I01Tr
      Top%I10Tr = I10Tr
      Top%I11Tr = I11Tr



      ! Auxiliar Vector (Passing from the nodal values of the time function to the gauss points)

      Top%Tmat_timefunc = transpose(Top%l)
      Top%Tmat_timefunct = transpose(Top%lt)



      IntTe = spread( Top%IntT(1:I%ngpt,1), 2, I%Nsft )  !reshape( Top%IntT( incr(1,I%ngpt) , 1 ), (/ I%ngpt , 1 /) )
      !print*, "INTe:"
      !call show( IntTe )
      Top%vwalfa =  transpose( IntTe*(Top%Tmat_timefunc)  )
      Top%vwalfat = transpose( IntTe*(Top%Tmat_timefunct) )


    end subroutine Time_op_assign




    subroutine Elm_calc_matrices(l,lt,waux,Q00k,Q01k,Q10k,Q11k) !result(Q00k,Q01k,Q11k)

        implicit none
        double precision, dimension(:,:) :: l, lt, waux
        double precision, dimension(:,:), allocatable :: wauxg
        double precision, dimension( size(l,1), size(l,2) )  :: lm,ltm
        double precision, dimension( size(l,2), size(l,1) ) ::  lmw, ltmw
        double precision, dimension( size(l,1), size(l,1) ) :: Q00k,Q01k,Q10k,Q11k

        allocate(wauxg( I%ngpt, I%Nsft ))
        wauxg = spread( reshape(waux, (/I%ngpt/) ) , 2, I%Nsft )

        !call show(wauxg)
        !print*, "espacio:"
        !call show(waux)
        !print*, "l:"
        !call show(l)

        !stop

        lm=l; ltm=lt
        lmw = ( wauxg*transpose(l) );  ltmw = ( wauxg*transpose(lt) )


        Q00k = matmul(lm,lmw)
        Q01k = matmul(ltm,lmw)
        Q10k = matmul(lm,ltmw)
        Q11k = matmul(ltm,ltmw)



    end subroutine Elm_calc_matrices



end module Time_matrices
