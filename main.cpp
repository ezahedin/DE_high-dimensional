#include "cmatrix.h"
#include "kron.h"
#include <math.h>
#include <iostream>
#include "utilities.h"
#include <time.h>
#include "DE.h"
#include "problem.h"
#include "EvalUT.h"
#include "printi.h"
#include <mpi.h>


int main(int argc, char* argv[])
{
    //==============================================================Initializing MPI Gloabl Variable=================================================================

    int totalnode,mynode;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&totalnode);
    MPI_Comm_rank(MPI_COMM_WORLD,&mynode);
    MPI_Status status;

    //==============================================================Initializing Variables for DE====================================================================
    unsigned int Generation=200000;
    float Crossover_CR=0.5;
    float Mutation_F=0.9;
    double NPDim_p=5.78;
    double ctrpar=27;
    //defining system paremeters
    unsigned int D=4;   //Dimension of each qubits
    unsigned int Computational_D=20;
    unsigned int Hamitonian_M=3;  //numner of control Hamiltonain
    float g=0.03;  //coupling factor
    double eta=-0.2; //anharmonocity
    double Time=27;
    double dt=2*M_PI*Time/(ctrpar);
    int mesh=10;
    //========================================================Declaring qubits' individual Hamiltonains (Real)===================================================

    matrix  etamat(D,D),X(D,D),eye(D,D),eye2(2,2),HXtot(pow(D,D),pow(D,D)),H1(D,D),EtaH(pow(D,D),pow(D,D));
    matrix  Hw1(pow(D,D),pow(D,D)),Hw2(pow(D,D),pow(D,D)),Hw3(pow(D,D),pow(D,D)),CNot(8,8),CNot_vec(1,8),Utarget(8,8),HYtot(pow(D,D),pow(D,D)),H(pow(D,D),pow(D,D));
    
    etamat(0,0)=0;etamat(0,1)=0;etamat(0,2)=0;etamat(0,3)=0;
    etamat(1,0)=0;etamat(1,1)=0;etamat(1,2)=0;etamat(1,3)=0;
    etamat(2,0)=0;etamat(2,1)=0;etamat(2,2)=eta;etamat(2,3)=0;
    etamat(3,0)=0;etamat(3,1)=0;etamat(3,2)=0;etamat(3,3)=3*eta;
    
    
    X(0,0)=0;X(0,1)=1;X(0,2)=0;X(0,3)=0;
    X(1,0)=1;X(1,1)=0;X(1,2)=sqrt(2);X(1,3)=0;
    X(2,0)=0;X(2,1)=sqrt(2);X(2,2)=0;X(2,3)=sqrt(3);
    X(3,0)=0;X(3,1)=0;X(3,2)=sqrt(3);X(3,3)=0;
    
    
    H1(0,0)=0;H1(0,1)=0;H1(0,2)=0;H1(0,3)=0;
    H1(1,0)=0;H1(1,1)=1;H1(1,2)=0;H1(1,3)=0;
    H1(2,0)=0;H1(2,1)=0;H1(2,2)=2;H1(2,3)=0;
    H1(3,0)=0;H1(3,1)=0;H1(3,2)=0;H1(3,3)=3;

    
    eye(0,0)=1;eye(0,1)=0;eye(0,2)=0;eye(0,3)=0;
    eye(1,0)=0;eye(1,1)=1;eye(1,2)=0;eye(1,3)=0;
    eye(2,0)=0;eye(2,1)=0;eye(2,2)=1;eye(2,3)=0;
    eye(3,0)=0;eye(3,1)=0;eye(3,2)=0;eye(3,3)=1;

    
    eye2(0,0)=1;eye2(0,1)=0;eye2(1,0)=0;eye2(1,1)=1;
    
    //================================================================Target Gate (CZ gate)====================================================================
    CNot_vec.array1d[0]=1;CNot_vec.array1d[1]=1;CNot_vec.array1d[2]=1;CNot_vec.array1d[3]=1;
    CNot_vec.array1d[4]=1;CNot_vec.array1d[5]=1;CNot_vec.array1d[6]=1;CNot_vec.array1d[7]=-1;
    CNot=diagmat(CNot_vec);
    //========================================================Declaring qubits' indovidual Hamiltonains (Complex)===================================================
    cmatrix Y(D,D);
    Y(0,0)=Complex(0);Y(0,1)=Complex(0,-1);Y(0,2)=Complex(0);Y(0,3)=Complex(0);
    Y(1,0)=Complex(0,1);Y(1,1)=Complex(0);Y(1,2)=Complex(0,-sqrt(2));Y(1,3)=Complex(0);
    Y(2,0)=Complex(0);Y(2,1)=Complex(0,sqrt(2));Y(2,2)=Complex(0);Y(2,3)=Complex(0,-sqrt(3));
    Y(3,0)=Complex(0);Y(3,1)=Complex(0);Y(3,2)=Complex(0,sqrt(3));Y(3,3)=Complex(0);
    
    //=================================================Define the total Hamiltonian along X,Y,Z of the system==========================================
    HYtot=g/2*(kron(kron(Y,Y),eye)+kron(eye,kron(Y,Y))).cmatrix_to_matrix();
    HXtot=g/2*(kron(kron(X,X),eye)+kron(eye,kron(X,X)));
    EtaH=kron(etamat,kron(eye,eye))+kron(kron(eye,eye),etamat)+kron(kron(eye,etamat),eye);
    H=HXtot+HYtot+EtaH;
    Hw1=kron(H1,kron(eye,eye));
    Hw2=kron(eye,kron(H1,eye));
    Hw3=kron(eye,kron(eye,H1));
    Utarget=CNot;
    mVector Htot(4,matrix(20,20));Htot[0]=projection(H);Htot[1]=projection(Hw1);Htot[2]=projection(Hw2);Htot[3]=projection(Hw3);
    mVector U(1,matrix(8,8)); U[0]=Utarget;
    //=====================================Define the problem to be sent to DE function, Pre-initialization of DE optimization===========================
    problem Problem(Hamitonian_M,Crossover_CR,Mutation_F,Computational_D,NPDim_p,ctrpar,Generation,Htot,U,dt,mesh);
    int dist1=1;
    int dist2=2;
    unsigned int NPDIM=(Problem.NPDim)*(Problem.K);
    unsigned int each_popsize=Problem.K*Problem.M;
	int seed[4];
if (fmod((double)mynode,2)==0) {
    	seed[0]=mynode+451;seed[1]=mynode+3425;seed[2]=mynode+2501;seed[3]=mynode+3217;
} else {
	seed[0]=mynode+1078;seed[1]=mynode+3960;seed[2]=mynode+1601;seed[3]=mynode+1202;
}
    matrix _pop(NPDIM,each_popsize),Fidelity(1,NPDIM),single_pop(1,each_popsize),Fidelity_itr(1,Generation);
    rand_matrix(dist2,seed,single_pop);
    
    double err=EvalUT(single_pop,Problem);

    if (mynode==0) {
        Fidelity.array1d[0]=err;
        for (int j=0; j<each_popsize; j++)
            _pop(mynode,j)=single_pop.array1d[j];
    }
    
    //=====================================================================DE rates definition===========================================================
    double tau1=.1;double tau2=.1;double Fl=.1;double Fu=.9;
    matrix rates(1,2);
    rand_matrix(dist1,seed,rates);
    rates.array1d[0]=rates.array1d[0]*Fu+Fl;

    //======================================================================MPI part of the program==========================================================
    //1-1 transferring data from salves to master-Initialization
    if (mynode!=0) {
        MPI_Send(single_pop.array1d,each_popsize,MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);
        MPI_Send(&err,1,MPI_DOUBLE, 0, 1,MPI_COMM_WORLD);
    } else {
        for (int j=1; j<totalnode; j++) {
            //1-2 recieving the population data
            MPI_Recv(single_pop.array1d,each_popsize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&status);
            for (int k=0; k<each_popsize; k++)
                _pop(status.MPI_SOURCE,k)=single_pop.array1d[k];
            //4-3 recieving the fidelity data, objective function output
            MPI_Recv(&err,1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&status);
            Fidelity.array1d[status.MPI_SOURCE]=err;
        }
    }

    int Gen=0;
    //======================================================================Iterative part of DE==========================================================

    while (Gen<Problem.G) {
        //2-1 transferring "allpop" and Fidelity to slaves
        if (mynode==0) {
            for(int j=1; j<totalnode; j++) {
                MPI_Send(_pop.array1d,NPDIM*each_popsize,MPI_DOUBLE, j, 0,MPI_COMM_WORLD);
                MPI_Send(Fidelity.array1d,NPDIM,MPI_DOUBLE, j, 0,MPI_COMM_WORLD);
            }
            //2-2 recieving "allpop" and Fidelity from slaves
            
        } else {
            MPI_Recv(_pop.array1d,NPDIM*each_popsize, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&status);
            MPI_Recv(Fidelity.array1d,NPDIM, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&status);
        }

        err=DE(Problem,_pop,single_pop,mynode,seed,Fidelity,rates,Gen);

        //3-1 running DE on each node and returning the results
        //updating master node information
        //saving data from node zero on global variable on node zero
        if (mynode==0) {
            Fidelity.array1d[0]=err;
            
            for(int j=0; j<each_popsize; j++)
                _pop(0,j)=single_pop.array1d[j];
        }
        //4 Updating allpop contents by single_pop from slaves
        if (mynode!=0) {
            MPI_Send(&err,1,MPI_DOUBLE, 0, 1,MPI_COMM_WORLD);
            MPI_Send(single_pop.array1d,each_popsize,MPI_DOUBLE, 0, 0,MPI_COMM_WORLD);

        } else {
            for (int j=1; j<totalnode; j++) {
                //6-1 recieving the fidelity data, objective function output
                MPI_Recv(&err,1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&status);
                         
                         Fidelity.array1d[status.MPI_SOURCE]=err;
                //6-2 recieving the population data
                MPI_Recv(single_pop.array1d,each_popsize, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&status);
                         
                         for (int k=0; k<each_popsize; k++)
                         _pop(status.MPI_SOURCE,k)=single_pop.array1d[k];
            }
            
            //printing out the results---from node zero
            
            Fidelity_itr.array1d[Gen]=(Fidelity.min()).array1d[0];
            
            //prints(single_pop);
        }
         Gen++;
    }
    
    
    //================================================================================TestArea=========================================================================
    //prints(X.rowi(2));
    //matrix test(2,3);
    //time_t time1=clock();
    //std::cout<<(double)(clock()-time1)/CLOCKS_PER_SEC<<'\n';

    
    if (mynode==0) {
        printi(Fidelity_itr,"Fidelity_itr.dat");
        printi(_pop.rowi((int)(Fidelity.min()).array1d[1]),"Optimal_solution.dat");
        printi(_pop,"Population.dat");
    }

    //matrix A(1,3),B(3,3);
    //eig(A,B,X);
    
       /*
    for (int i=0; i<eye3.rows; i++) {
        for (int j=0; j<eye3.cols; j++) {
            std::cout<<i+1<<"  "<<j+1<<"   "<<eye3(i,j)<<'\n';;
        }
    }



       */
    MPI_Finalize();


    
    
        return 0;
}


