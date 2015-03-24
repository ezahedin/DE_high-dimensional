#include "problem.h"
#include <iostream>
using namespace std;

problem::problem(unsigned int h1,float h2,float h3,unsigned int h4,double h5,double h6,unsigned int h7, mVector h8,mVector h9,double h10,int h11) {

    M=h1;
    CR=h2;
    F=h3;
    dimension=h4;
    NPDim=h5;
    K=h6;
    G=h7;
    H.push_back(h8[0]);
    H.push_back(h8[1]);
    H.push_back(h8[2]);
    H.push_back(h8[3]);
    U.push_back(h9[0]);
    dt=h10;
    mesh=h11;
    
    //U=h9;
}


//[h1,h2,h3,h4,h5,h6,h7,h8,h9]=[M,CR,F,dimesnion,NPDim,K,G,H,U]
//Remember: I have to define the destructore for this class in the future.
//Although it is not needed at this time since We only generate One problem for one run and it
//willbe destryed at the end of the run.