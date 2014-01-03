all: NumSim

NumSim: computation.o gridfunction.o IO.o NumSimFluidPara.o solver.o stencil.o communication.o derivatives.o
	mpic++ computation.o gridfunction.o IO.o NumSimFluidPara.o solver.o stencil.o communication.o derivatives.o -o NumSim

computation.o: computation.cpp
	mpicxx -c computation.cpp

gridfunction.o: gridfunction.cpp
	mpicxx -c gridfunction.cpp

IO.o: IO.cpp
	mpicxx -c IO.cpp

NumSimFluidPara.o: NumSimFluidPara.cpp
	mpicxx -c NumSimFluidPara.cpp

solver.o: solver.cpp
	mpicxx -c solver.cpp

stencil.o: stencil.cpp
	mpicxx -c stencil.cpp

communication.o: communication.cpp
	mpicxx -c communication.cpp

derivatives.o: derivatives.cpp
	mpicxx -c derivatives.cpp

clean:
	rm -rf *o hello
