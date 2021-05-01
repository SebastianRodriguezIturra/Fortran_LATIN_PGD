PROGRAM Parallel_Hello_World
USE OMP_LIB
!OMP_NUM_THREADS=4;
!export OMP_NUM_THREADS
!$OMP PARALLEL

    PRINT *, 'Hello from process:', OMP_GET_THREAD_NUM()

!$OMP END PARALLEL

!!gfortran -o hello parallel_hello_world.f03 -fopenmp (to compile)


END
