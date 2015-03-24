# ****************************************************
#Variable to control make file operation
ICC=mpiCC
IXX=g++ 
mpi=icpc
INC1=-I$(MKLROOT)/include
#INC1=-I/global/software/intel/composerxe-2011-sp1.10.319/mkl/include  -I/global/software/include
#LIB1=-L/global/software/intel/composerxe-2011-sp1.10.319/mkl/lib/intel64

LIB2= -mkl=sequential

runfile: main.o matrix.o cmatrix.o kron.o problem.o utilities.o DE.o EvalUT.o printi.o
	$(ICC) $(LIB2) printi.o EvalUT.o DE.o utilities.o problem.o kron.o cmatrix.o matrix.o main.o -o runfile -O3

main.o: main.cpp matrix.h
	$(ICC) -c main.cpp -O3

matrix.o:matrix.cpp matrix.h
	$(ICC) -c matrix.cpp  -O3

cmatrix.o:cmatrix.cpp
	$(CXX) -c cmatrix.cpp -O3

kron.o:kron.cpp
	$(CXX) -c kron.cpp -O3

problem.o:problem.cpp
	$(CXX) -c problem.cpp

utilities.o:utilities.cpp
	$(CXX) -c utilities.cpp

DE.o:DE.cpp
	$(CXX) -c DE.cpp

EvalUT.o:EvalUT.cpp
	$(CXX) -c EvalUT.cpp

printi.o:printi.cpp
	$(CXX) -c printi.cpp


clean :
	rm ./*.o ./*~ p1
