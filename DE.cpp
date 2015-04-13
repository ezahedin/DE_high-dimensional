#include "DE.h"
double DE(problem& Problem,matrix& _pop,matrix& single_pop,int mynode,int* seed,matrix& Fidelity,matrix& rates,int Gen) {
    //===================================================Initializing DE parameters============================================================
    //Population size NPDIM
    unsigned int NPDIM=(Problem.NPDim)*(Problem.K);
    unsigned int each_popsize=Problem.K*Problem.M;
    //mutatted population
    matrix _mpop(1,each_popsize);
    _mpop=single_pop;
    //Crossed over elements
    matrix _cpop(1,each_popsize);
    matrix index(1,each_popsize);
    //=============================================Random number for Crossover and Mutation evolution======================================
    int dist1=1;
    int dist2=2;
    matrix temprand(1,4);
    rand_matrix(dist1,seed,temprand);
    //===================================================Step One-- Randomly Initialize the population (_pop)==================================
    int indexid[NPDIM],idrand[each_popsize];
    for (int i=0; i<NPDIM; i++) {
        indexid[i]=i;
    }

    for (int i=0; i<each_popsize; i++) {
        idrand[i]=i;
    }
    //===================================Define the Self-adaptive constants and and self adaptive array containing the crossover rate and mutation rate
    double tau1=.1;double tau2=.1;double Fl=.1;double Fu=.9;
    //=============decide whether we wanna do the single element or all elements optimization===================
        random_shuffle1(idrand,each_popsize,(mynode+1)*seed[0]);
         int S=1;
        if (fmod(double(Gen),7)==0) {
            S=each_popsize;
        }
    //================================mutation part using mutation rate, F==========================================================================

            random_shuffle1(indexid,NPDIM,mynode*seed[2]);
            if ((indexid[0]==mynode)|(indexid[1]==mynode)|(indexid[2]==mynode)) {
                int tmp1=indexid[0];int tmp2=indexid[1];int tmp3=indexid[2];
                indexid[0]=indexid[3];indexid[1]=indexid[4];indexid[2]=indexid[5];
                indexid[3]=tmp1;indexid[4]=tmp2;indexid[5]=tmp3;
            }

    
            for (int j=0; j<S; ++j)
                _mpop.array1d[idrand[idrand[j]]]=_pop(indexid[0],idrand[idrand[j]])+rates.array1d[0]*(_pop(indexid[1],idrand[idrand[j]])-_pop(indexid[2],idrand[idrand[j]]));
    //==========================if this is a constrained problem run this part, Lets write the function=============================================
        constrain(_mpop,-2.5,2.5);
    //=============================================Crossover part of the DE using Cross over rate, G ===============================================
    _cpop=_mpop;
    
            rand_matrix(dist1,seed,index);
            for (int j=0; j<S; ++j) {
                if (index.array1d[j]<rates.array1d[1]) {
                    _cpop.array1d[idrand[idrand[j]]]=single_pop.array1d[idrand[idrand[j]]];
                }

            }
     
    //============================================================Selection part of the DE============================================================
    double err=Fidelity.array1d[mynode];
    double tmp=EvalUT(_cpop,Problem);
    
            if (tmp<Fidelity.array1d[mynode]) {
                  err=tmp;
                single_pop=_cpop;}
    //============================================================Updating the DE constants rates===========================================================
    rates.array1d[0]=(temprand.array1d[0]<tau1) ? (Fl+temprand.array1d[1]*Fu) : rates.array1d[0];
    rates.array1d[1]=(temprand.array1d[2]<tau2) ? (temprand.array1d[3]) : rates.array1d[1];
    
    return err;
}


