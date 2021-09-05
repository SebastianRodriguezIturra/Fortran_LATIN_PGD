module gnuplot_fortran

    use Matrix_op

    implicit none

    interface plot
      module procedure plot1i, plot2i, plot1r, plot2r
    end interface




    contains



      subroutine  plot1i(y) ! only one argument
        integer, intent(in), dimension(:) :: y
        integer, allocatable , dimension(:) :: x ! vector of x
        integer :: size_y, i
        size_y = size(y)

        ! creation of x:
        allocate(x(size_y))
        x = incr(1,size_y)

        open( unit=1, file = 'plot.dat' )
        do i = 1,size_y
           write(1,*) x(i), ' ', y(i)
        end do


        call system('gnuplot -p plot.plt')

      end subroutine plot1i





      subroutine  plot1r(y) ! only one argument
        double precision, intent(in), dimension(:) :: y
        integer, allocatable , dimension(:) :: x ! vector of x
        integer :: size_y, i
        size_y = size(y)

        ! creation of x:
        allocate(x(size_y))
        x = incr(1,size_y)

        open( unit=1, file = 'plot.dat' )
        do i = 1,size_y
           write(1,*) x(i), ' ', y(i)
        end do


        call system('gnuplot -p plot.plt')

      end subroutine plot1r




      subroutine  plot2i(x,y) ! two arguments
        integer, intent(in), dimension(:) :: x,y
        integer :: size_x, size_y, i
        size_x = size(x)
        size_y = size(y)
        if ( size_x /= size_y ) then
          print*, "Array size missmatch"
        else
          open( unit=1, file = 'plot.dat' )
          do i = 1,size(x)
             write(1,*) x(i), ' ', y(i)
          end do
        end if

        call system('gnuplot -p plot.plt')

      end subroutine plot2i





    subroutine  plot2r(x,y) ! two arguments
      double precision, intent(in), dimension(:) :: x,y
      integer :: size_x, size_y, i
      size_x = size(x)
      size_y = size(y)
      if ( size_x /= size_y ) then
        print*, "Array size missmatch"
      else
        open( unit=1, file = 'plot.dat' )
        do i = 1,size(x)
           write(1,*) x(i), ' ', y(i)
        end do
      end if

      call system('gnuplot -p plot.plt')

    end subroutine plot2r









end module gnuplot_fortran
