module Matrix_op



  implicit none



  interface full_field
    module procedure full_field_vector, full_field_tensor
  end interface

  interface show
    module procedure show1i, show2i, show1d, show2d
  end interface

  interface kron
    module procedure kron2ri,kron2ir,kron2i,kron2r
  end interface

  interface ones
    module procedure ones1,ones2
  end interface

  interface inv
    module procedure invs,invm
  end interface


  interface horzcat
    module procedure h_1,h_2
  end interface

  interface vertcat
    module procedure v_1real,v_1int,v_2real,v_2int
  end interface

  interface mmat
    module procedure mmat1, mmat2
  end interface

  !interface incr
  !  module procedure incr1,incr2
  !end interface



  contains





!!!!! full_field !!!!!

function full_field_vector(m1,m2) result(m)

  ! Function that calculate from the PGD to the full solution:
  !integer, dimension(:), allocatable :: siz, new_siz
  integer :: n,Nmode,Tngpt
  double precision, dimension(:,:) :: m1,m2
  double precision, dimension(:,:), allocatable :: m  ! Displacement result.

  n = size(m1,1)
  Nmode = size(m1,2)
  Tngpt = size(m2,1)

  allocate( m(n,Tngpt) )
  m = matmul( m1, transpose(m2)  )

end function full_field_vector




function full_field_tensor(m1,m2) result(m)

  ! Function that calculate from the PGD to the full solution:
  integer, dimension(:), allocatable :: siz, new_siz
  integer :: Nmode,Tngpt
  double precision, dimension(:,:,:,:) :: m1
  double precision, dimension(:,:) :: m2
  double precision, dimension(:,:), allocatable :: mi
  double precision, dimension(:,:,:,:), allocatable :: m  ! Displacement result.

  Nmode = size(m2,2)
  Tngpt = size(m2,1)

  allocate( siz(4) )
  allocate( new_siz(4) )

  siz = shape(m1)
  new_siz = siz; new_siz(4) = Tngpt

  allocate(  m(new_siz(1), new_siz(2), new_siz(3), new_siz(4))  )
  allocate( mi(size(m1)/Nmode,Tngpt) )

  mi = matmul( reshape(m1,  (/size(m1)/Nmode, Nmode/) ), transpose(m2) )
  m  = reshape( mi , (/new_siz(1), new_siz(2), new_siz(3), new_siz(4)/) )

end function full_field_tensor




!!!!! repmat 2D: !!!!!


function repmat(m,vecdim) result(r)

  implicit none
  double precision, dimension(:,:) :: m
  double precision, dimension(:,:), allocatable :: r
  integer, dimension(2) :: vecdim, sizem , sizef
  integer :: i,jr,j,ir

  sizem = shape(m)

  do i = 1,2 !(int i=0;i<siz.length;i++){
      sizef(i) = sizem(i)*vecdim(i)   ! sizf[i] = siz[i]*repdim[i];
  end do ! determine the final size desired after repmat. !} // determine the final size desired after repmat.


  allocate(r(sizef(1),sizef(2)))

  do jr = 1, vecdim(2) !for (int jr = 1; jr <= repdim[1]; jr++) {
    do j = 1, sizem(2) ! for (int j = 0; j < siz[1]; j++) {
      do ir = 1, vecdim(1) !for (int ir = 1; ir <= repdim[0]; ir++) {
        do i = 1,sizem(1)   ! for (int i = 0; i < siz[0]; i++) {
          r( (ir-1)*sizem(1)+i,(jr-1)*sizem(2)+j ) = m(i,j) ! result[ (ir-1)*siz[0]+i ][ (jr-1)*siz[1]+j ] = m[i][j];
        end do !} // 1st dim
      end do !} // 1st rep
    end do  !} // 2nd dim
  end do !} // 2nd rep


end function repmat


!public double[][] repmat(double m[][],int[] repdim){
!
!    int[] siz = size(m);
!    int[] sizf = new int[siz.length];
!    for (int i=0;i<siz.length;i++){
!        sizf[i] = siz[i]*repdim[i];
!    } // determine the final size desired after repmat.
!
!    double result[][] = new double[sizf[0]][sizf[1]];
!
!
!    for (int jr = 1; jr <= repdim[1]; jr++) {
!        for (int j = 0; j < siz[1]; j++) {
!            for (int ir = 1; ir <= repdim[0]; ir++) {
!                for (int i = 0; i < siz[0]; i++) {
!                    result[ (ir-1)*siz[0]+i ][ (jr-1)*siz[1]+j ] = m[i][j];
!                } // 1st dim
!            } // 1st rep
!        } // 2nd dim
!    } // 2nd rep

!    return result;
!} // end repmat 2D





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function kron2ri(m1,m2) result(m)
  implicit none

  double precision, dimension(:,:) :: m1 ! input
  integer, dimension(:,:) :: m2 ! input


  integer :: d1
  integer :: d2
  integer :: row1, column1, row2, column2
  integer :: i,j,k,l
  integer, allocatable, dimension(:) :: ind1, ind2 ! index that change over the iterations.
  double precision, allocatable, dimension(:,:) :: m ! output


  row1 = size(m1,1);  column1 = size(m1,2)
  row2 = size(m2,1);  column2 = size(m2,2)

  d1 = size(m1,1)*size(m2,1)
  d2 = size(m1,2)*size(m2,2)

  allocate( m(d1,d2) )

  ind1 = incr(1,size(m2,1)); ind2 = incr(1,size(m2,2))

  !call show(ind1)  ! shows the vector

  do i = 1,row1
    do j = 1,column1
      do k = 1,row2
        do l = 1,column2
          m((i-1)*row2+k,(j-1)*column2+l) =  m1(i,j)*m2(k,l)
        end do
      end do
    end do
  end do


end function kron2ri





  function kron2ir(m1,m2) result(m)
    implicit none

    integer, dimension(:,:) :: m1 ! input
    double precision, dimension(:,:) :: m2 ! input


    integer :: d1
    integer :: d2
    integer :: row1, column1, row2, column2
    integer :: i,j,k,l
    integer, allocatable, dimension(:) :: ind1, ind2 ! index that change over the iterations.
    double precision, allocatable, dimension(:,:) :: m ! output


    row1 = size(m1,1);  column1 = size(m1,2)
    row2 = size(m2,1);  column2 = size(m2,2)

    d1 = size(m1,1)*size(m2,1)
    d2 = size(m1,2)*size(m2,2)

    allocate( m(d1,d2) )

    ind1 = incr(1,size(m2,1)); ind2 = incr(1,size(m2,2))

    !call show(ind1)  ! shows the vector

    do i = 1,row1
      do j = 1,column1
        do k = 1,row2
          do l = 1,column2
            m((i-1)*row2+k,(j-1)*column2+l) =  m1(i,j)*m2(k,l)
          end do
        end do
      end do
    end do


  end function kron2ir





  function kron2i(m1,m2) result(m)
    implicit none

    integer, dimension(:,:) :: m1,m2 ! input
    integer :: d1
    integer :: d2
    integer :: row1, column1, row2, column2
    integer :: i,j,k,l
    integer, allocatable, dimension(:) :: ind1, ind2 ! index that change over the iterations.
    double precision, allocatable, dimension(:,:) :: m ! output


    row1 = size(m1,1);  column1 = size(m1,2)
    row2 = size(m2,1);  column2 = size(m2,2)

    d1 = size(m1,1)*size(m2,1)
    d2 = size(m1,2)*size(m2,2)

    allocate( m(d1,d2) )

    ind1 = incr(1,size(m2,1)); ind2 = incr(1,size(m2,2))

    !call show(ind1)  ! shows the vector

    do i = 1,row1
      do j = 1,column1
        do k = 1,row2
          do l = 1,column2
            m((i-1)*row2+k,(j-1)*column2+l) =  m1(i,j)*m2(k,l)
          end do
        end do
      end do
    end do


  end function kron2i





function kron2r(m1,m2) result(m)
  implicit none

  double precision, dimension(:,:) :: m1,m2 ! input
  integer :: d1
  integer :: d2
  integer :: row1, column1, row2, column2
  integer :: i,j,k,l
  integer, allocatable, dimension(:) :: ind1, ind2 ! index that change over the iterations.
  double precision, allocatable, dimension(:,:) :: m ! output


  row1 = size(m1,1);  column1 = size(m1,2)
  row2 = size(m2,1);  column2 = size(m2,2)

  d1 = size(m1,1)*size(m2,1)
  d2 = size(m1,2)*size(m2,2)

  allocate( m(d1,d2) )

  ind1 = incr(1,size(m2,1)); ind2 = incr(1,size(m2,2))

  !call show(ind1)  ! shows the vector

  do i = 1,row1
    do j = 1,column1
      do k = 1,row2
        do l = 1,column2
          m((i-1)*row2+k,(j-1)*column2+l) =  m1(i,j)*m2(k,l)
        end do
      end do
    end do
  end do


end function kron2r


!!!!!!!!!!!!!!!!!!


function ones1(s1) result(v)
  implicit none
  integer :: s1
  integer, dimension(s1) :: v ! output
  v = 1
end function ones1

function ones2(s1,s2) result(v)
  implicit none
  integer :: s1,s2
  integer, dimension(s1,s2) :: v ! output
  v = 1
end function ones2




! mmat:

  function mmat1(m1,m2) result(m)
    implicit none

    double precision, dimension(:,:,:) :: m1,m2 ! input
    integer :: d3
    integer :: i
    double precision, dimension(size(m1,1),size(m2,2),size(m1,3)) :: m ! output
    !double precision :: aux,auxp

    d3 = size(m1,3) ! third dimension of the matrices (assume same d3 for both)

    do i = 1,d3
        m(:,:,i) = matmul(m1(:,:,i),m2(:,:,i))
    end do

  end function mmat1




  function mmat2(m1,m2) result(m)
    implicit none

    double precision, dimension(:,:) :: m1 ! input
    double precision, dimension(:,:,:) :: m2 ! input
    integer :: d3
    integer :: i
    double precision, dimension(size(m1,1),size(m2,2),size(m2,3)) :: m ! output
    !double precision :: aux,auxp

    d3 = size(m2,3) ! third dimension of the matrices (assume same d3 for both)

    do i = 1,d3
        m(:,:,i) = matmul(m1,m2(:,:,i))
    end do

  end function mmat2


! sort:

recursive subroutine sort(array)
    double precision, intent(inout)::array(:)
    double precision :: temp,pivot
    integer :: i,j,last,left,right

    last=size(array)

    if (last.lt.50) then ! use insertion sort on small arrays
       do i=2,last
          temp=array(i)
          do j=i-1,1,-1
             if (array(j).le.temp) exit
             array(j+1)=array(j)
          enddo
          array(j+1)=temp
       enddo
       return
    endif
    ! find median of three pivot
    ! and place sentinels at first and last elements
    temp=array(last/2)
    array(last/2)=array(2)
    if (temp.gt.array(last)) then
       array(2)=array(last)
       array(last)=temp
    else
       array(2)=temp
    endif
    if (array(1).gt.array(last)) then
       temp=array(1)
       array(1)=array(last)
       array(last)=temp
    endif
    if (array(1).gt.array(2)) then
       temp=array(1)
       array(1)=array(2)
       array(2)=temp
    endif
    pivot=array(2)

    left=3
    right=last-1
    do
       do while(array(left).lt.pivot)
          left=left+1
       enddo
       do while(array(right).gt.pivot)
          right=right-1
       enddo
       if (left.ge.right) exit
       temp=array(left)
       array(left)=array(right)
       array(right)=temp
       left=left+1
       right=right-1
    enddo
    if (left.eq.right) left=left+1
    call sort(array(1:left-1))
    call sort(array(left:))

  end subroutine sort












! inv matrix:

  function  invs(ms) result(minv)
    implicit none
    double precision ::  ms ! input
    double precision ::  minv ! output

    minv = 1/ms
  end function invs


  function  invm(ms) result(minv)
    implicit none
    ! double == double precision (arreglar esto, aumentar precision!)
    !real, intent(in), dimension(:,:) :: m ! input
    double precision, dimension(:,:) :: ms ! input
    integer :: d1,d2,i,i1,j
    double precision, dimension(size(ms,1),size(ms,2)) :: minv,m ! output
    double precision :: aux,auxp


    d1 = size(ms,1);d2 = size(ms,2)
    m = ms

    ! Creation of the identity matrix:
    minv = 0
    do i=1,d1
      minv(i,i) = 1
    end do


    ! Calculation of the Inverse matrix:

    do i=1,d1
      aux = m(i,i)
      do j=1,d2
        m(i,j) = m(i,j)/aux
        minv(i,j) = minv(i,j)/aux
      end do ! j

      do i1 = 1,d1
        if (i1 /= i) then ! i1 not equal to i
          auxp = -m(i1,i);

          do j=1,d2
            m(i1,j) = m(i1,j) + m(i,j)*auxp;
            minv(i1,j) = minv(i1,j) + minv(i,j)*auxp;
          end do

        end if
      end do ! i1

    end do ! i


  end function invm



  function  linspace(vi,vf,N) result(m)
    implicit none
    double precision :: vi,vf,dx
    integer :: N,i
    double precision, dimension(N) :: m ! output

    dx = (vf-vi)/(N-1)
    m = (/(i, i=0,N-1, 1)/)*dx

  end function linspace


  ! integer array formed: (equivalent matlab a:b)

  function  incr(a,b) result(m)
    implicit none
    integer :: a,b,i,N
    integer, allocatable, dimension(:) :: m ! output

    N = b-a+1
    allocate( m(N) )  != (/(i, i=a, b, 1)/)
    m = (/(i, i=a, b, 1)/)

  end function incr


  function  incr2(a,b) result(m)
    implicit none
    integer :: a,b,i,N
    integer, allocatable, dimension(:,:) :: m ! output

    N = b-a+1
    allocate( m(1,N) )  != (/(i, i=a, b, 1)/)
    m(1,:) = (/(i, i=a, b, 1)/)
    !DEALLOCATE(m)
  end function incr2






! horzcat:

function h_1(m1,m2) result(m)
  implicit none
  double precision, dimension(:) :: m1,m2 ! input
  double precision, dimension(size(m1)+size(m2)) :: m ! output

  m = [m1,m2]
  !m(1:size(m1,1)) = m1
  !m(size(m1,1)+1:(size(m1,1)+size(m2,1))) = m2

end function h_1


function h_2(m1,m2) result(m)
  implicit none
  double precision, dimension(:,:) :: m1,m2 ! input
  double precision, dimension(size(m1,1),size(m1,2)+size(m2,2)) :: m ! output


  m(:,1:size(m1,2)) = m1
  m(:,size(m1,2)+1:size(m1,2)+size(m2,2)) = m2

end function h_2




! vertcat:

function v_1real(m1,m2) result(m)
  implicit none
  double precision, dimension(:) :: m1,m2 ! input
  double precision, dimension(2,size(m1)) :: m ! output


  m(1,:) = m1
  m(2,:) = m2

end function v_1real


function v_1int(m1,m2) result(m)
  implicit none
  integer, dimension(:) :: m1,m2 ! input
  integer, dimension(2,size(m1)) :: m ! output


  m(1,:) = m1
  m(2,:) = m2

end function v_1int




function v_2real(m1,m2) result(m)
  implicit none
  double precision, dimension(:,:) :: m1,m2 ! input
  double precision, dimension(size(m1,1)+size(m2,1),size(m1,2)) :: m ! output


  m(1:size(m1,1),:) = m1
  m(size(m1,1)+1:size(m1,1)+size(m2,1),:) = m2

end function v_2real


function v_2int(m1,m2) result(m)
  implicit none
  integer, dimension(:,:) :: m1,m2 ! input
  integer, dimension(size(m1,1)+size(m2,1),size(m1,2)) :: m ! output


  m(1:size(m1,1),:) = m1
  m(size(m1,1)+1:size(m1,1)+size(m2,1),:) = m2

end function v_2int



!!!!!!!!


subroutine show1i(m)
  implicit none
  integer, dimension(:) :: m ! input
  integer :: row

  do row = 1,size(m)
    print*, m(row)
  end do

end subroutine show1i


subroutine show2i(m)
  implicit none
  integer, dimension(:,:) :: m ! input
  integer :: row

  do row = 1,size(m,1)
    print*, m(row,:)
  end do

end subroutine show2i





subroutine show1d(m)
  implicit none
  double precision, dimension(:) :: m ! input
  integer :: row

  do row = 1,size(m)
    print*, m(row)
  end do

end subroutine show1d


subroutine show2d(m)
  implicit none
  double precision, dimension(:,:) :: m ! input
  integer :: row

  do row = 1,size(m,1)
    print*, m(row,:)
  end do

end subroutine show2d







end module Matrix_op
