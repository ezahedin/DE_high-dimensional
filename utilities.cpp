#include "utilities.h"
//=======================================================================
void rand_matrix(int a,int* seed,matrix& t1)
{
    
    int rows=t1.rows;
    int cols=t1.cols;
    int size=rows*cols;
    dlarnv (&a, seed, &size,t1.array1d);
    
}

//=======================================================================
void random_shuffle1(int* mat,int size,int mynode)
{
    srand (unsigned(time(0)*mynode));
    std::random_shuffle( &mat[0], &mat[size]);
   /* for (int i=0; i<size; i++)
        std::cout<<mat[i]<<'\n';*/
}

//=======================================================================
void constrain(matrix& pop,double min,double max)
{
    int rows=pop.rows;
    int cols=pop.cols;
    int size=cols*rows;
    for (int i=0; i<size; i++) {
        if(pop.array1d[i]>max){
            pop.array1d[i]=max;
        } else if (pop.array1d[i]<min) {
            pop.array1d[i]=min;
        }
    }
}
//=======================================================================
matrix projection(matrix& t1)
{
    int G[20]={0 ,1, 2, 3, 4, 5, 6, 8, 9, 12, 16, 17, 18, 20, 21, 24, 32, 33, 36, 48};
    matrix tmp(20,20);
    
    for (int k=0; k<20; k++)
        for (int j=0; j<20; j++)
            tmp(k,j)=t1(G[k],G[j]);
    
    return tmp;
}

//=======================================================================
matrix erfun(problem& Problem,matrix& f)
{
    int K=Problem.K;
    int M=Problem.M;
    int mesh=Problem.mesh;
    int dim=(K-1)*mesh;
    double Pdt=Problem.dt;
    double width=5;
    matrix pulse1(1,K),pulse2(1,K),pulse3(1,K);
    matrix sp1(1,dim),sp2(1,dim),sp3(1,dim);
    matrix out(1,3*dim);
    
    for (int main=0; main<K; main++) {
        pulse1.array1d[main]=f.array1d[main*M];
        pulse2.array1d[main]=f.array1d[main*M+1];
        pulse3.array1d[main]=f.array1d[main*M+2];
    }
    
    
    double dt=(double)1/mesh*2*M_PI;
    int num=0;
    double X,Y1p1,Y1p2,Y1p3,Y2p1,Y2p2,Y2p3;
    for (int i=0; i<K-1; i++) {
        X=((double)(i+0.5)*Pdt)+(dt/2);
        Y1p1=0.5*(pulse1.array1d[i]+pulse1.array1d[i+1]);Y2p1=0.5*(pulse1.array1d[i+1]-pulse1.array1d[i]);
        Y1p2=0.5*(pulse2.array1d[i]+pulse2.array1d[i+1]);Y2p2=0.5*(pulse2.array1d[i+1]-pulse2.array1d[i]);
        Y1p3=0.5*(pulse3.array1d[i]+pulse3.array1d[i+1]);Y2p3=0.5*(pulse3.array1d[i+1]-pulse3.array1d[i]);
        
        for (int j=0; j<mesh; j++) {
            double err=erf((width*(X-0.5*(Pdt*2*(i+1))))/(Pdt));
            sp1.array1d[num]=Y1p1+Y2p1*err;
            sp2.array1d[num]=Y1p2+Y2p2*err;
            sp3.array1d[num]=Y1p3+Y2p3*err;
            num+=1;
            X+=dt;
        }
    }
    
    num=0;
    for (int k=0; k<dim; k++) {
        out.array1d[num]=sp1.array1d[k];
        out.array1d[num+1]=sp2.array1d[k];
        out.array1d[num+2]=sp3.array1d[k];
        num+=3;
    }
    
//printi(out,"out.dat");

    return out;
}//End of function

