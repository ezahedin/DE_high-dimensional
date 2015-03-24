#ifndef PROBLEM_H
#define PROBLEM_H
#include <complex>
#include "cmatrix.h"
#include "matrix.h"
#include <math.h>
#include <vector>

typedef std::vector<matrix> mVector;
typedef std::vector<cmatrix> cmVector;


class problem {
public:
    //M defines the number of control Hamitonian
    // Crossover rate
    // F mutation rate
    //dimension: is the Hilbert-space dimension
    //NPDim is the number that if multiply with K gives us the number of population
    //G Number of Generation
    //K number of control parameter
    unsigned int M;                                                           //h1
    float CR;                                                                 //h2
    float F;                                                                  //h3
    unsigned int dimension;                                                   //h4
    unsigned int NPDim;                                                       //h5
    int K;                                                                    //h6
    unsigned int G;                                                           //h7
    double dt;
    //Having Hamiltonian defined Here. Their number may vary as the number of control Hamiltonian changes
    mVector H;   //h8
    
    //matrix* H1,*H2,*H3,*U;
    //Target Matrix
    mVector U;                                                            //h9
    int mesh;
public:
    //Constructor: Here we construct all the element of problem [h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11]=[M,CR,F,dimesnion,NPDim,K,G,H,U,dt,mesh]
problem(unsigned int h1,float h2,float h3,unsigned int h4,double h5,double h6,unsigned int h7, mVector h8,mVector h9,double h10,int h11);
    //void kronker(complex** h9,complex** h10);
};

#endif //PROBLEM_H


