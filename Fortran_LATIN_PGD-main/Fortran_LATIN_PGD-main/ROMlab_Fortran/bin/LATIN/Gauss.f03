module Gauss


  contains

  function gauss_points_weights(ngps) result(vect_gauss)

   !! Gauss points and weights
   ! Give the exact integral for 2n-1 degree of polinomials.
   ! n: Number of gauss's points.
   implicit none
   integer :: ngps
   double precision :: term1, term2, term3, term4
   !integer, dimension(:) :: sizeg
   double precision, dimension(ngps,2) :: vect_gauss

   !sizeg = (/ ngps, ngps /)

   if (ngps==1) then

      vect_gauss(1,:) = (/0.D0 , 2.D0 /) ![[0],[2]];

   else if (ngps==2) then

      vect_gauss(:,1) =  (/ -sqrt(1/3.D0) , sqrt(1/3.D0) /)                     ! [[-sqrt(1/3);sqrt(1/3)],[1;1]];
      vect_gauss(:,2) =  (/ 1.D0, 1.D0 /)

   else if (ngps==3) then

      vect_gauss(:,1) =  (/ -sqrt(3/5.D0), 0.D0, sqrt(3/5.D0) /)
      vect_gauss(:,2) =  (/ 5/9.D0, 8/9.D0, 5/9.D0 /)

      !vect_gauss=[[-sqrt(3/5);0;sqrt(3/5)],[5/9;8/9;5/9]];
   else !(ngps==4)

      term1 = -sqrt((3+2*sqrt(6/5.D0))/7)
      term2 = -sqrt((3-2*sqrt(6/5.D0))/7)
      term3 = sqrt((3-2*sqrt(6/5.D0))/7)
      term4 = sqrt((3+2*sqrt(6/5.D0))/7)

      vect_gauss(:,1) =  (/ term1 , term2 , term3 , term4 /)
      vect_gauss(:,2) =  (/ (18-sqrt(30.D0))/36 , (18+sqrt(30.D0))/36 , (18+sqrt(30.D0))/36 , (18-sqrt(30.D0))/36 /)

!(18-sqrt(30))/36;   (18-sqrt(30))/36 ;      (18+sqrt(30))/36 ;       (18+sqrt(30))/36
!      [2,4,3,1]
!        gaussn4 =[[sqrt((3+2*sqrt(6/5))/7);-sqrt((3+2*sqrt(6/5))/7);sqrt((3-2*sqrt(6/5))/7);-sqrt((3-2*sqrt(6/5))/7)],[(18-sqrt(30))/36;(18-sqrt(30))/36;(18+sqrt(30))/36;(18+sqrt(30))/36]];
        ![ord,I1]=sort(gaussn4(:,1));
        !I.gauss_s=gaussn4(I1,:);

   end if


 end function gauss_points_weights


end module Gauss
