
VPATH = bin:bin/LATIN:bin/Input:bin/gnu_plot_use

FC = gfortran
FFLAGS=-O3 -march=native
SRC2 = Matrix_op.f03 Structure.f03 gnuplot_fortran.f03 Gauss.f03
#SRC2 = Input_ROMlab.f03 I_structure.f03 gnuplot_fortran.f03 Matrix_op.f03
SRC1 = bar.f03 Time_matrices.f03 Elastic_sol.f03 Elastic_sol_PGD.f03
SRC = local_stage.f03 local_stage_vec.f03 global_stage.f03  LATIN_indicator.f03 main.f03
OBJ = ${SRC2:.f03=.o} ${SRC1:.f03=.o} ${SRC:.f03=.o}
BLIBS1 = -L/usr/lib/x86_64-linux-gnu/blas -lblas
BLIBS2 = -L/usr/lib/x86_64-linux-gnu/lapack -llapack




%.o: %.f03
	$(FC) $(FFLAGS) -o $@ -c $< -fopenmp


maths: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(BLIBS1) $(BLIBS2) -fopenmp
	rm *.mod *.o
	./maths


clean: maths
	rm *.mod *.o maths
	echo Clean done
