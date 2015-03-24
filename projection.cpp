#include "projection.h"

matrix projection(matrix& t1)
{
    int G[20]={0 1 2 3 4 5 6 8 9 12 16 17 18 20 21 24 32 33 36 48};
    matrix tmp(20,20);
    
    for (int k=0; k<20; k++)
        for (int j=0; j<20; j++)
            tmp(k,j)=t1(G[k],G[k]);

        return tmp;
    }
}


