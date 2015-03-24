#include "EvalUT.h"
#include "iostream"

double EvalUT(matrix m,problem& Problem)
{
    int N=Problem.dimension;
    int M=Problem.M;
    int K=Problem.K;
    int mesh=Problem.mesh;
    matrix m1=erfun(Problem,m);
    cmatrix Utot(N,N);Utot.eye();
    matrix V(N,N),Er(1,N);
    cmatrix Ec(1,N);
    double INTEGDT=(double)1/mesh*2*M_PI;
    for (int main=0; main<(K-1)*mesh; main++) {

        matrix X=Problem.H[0];
        

        for (int i=0; i<M; i++) {
            //std::cout<<m.array1d[M*main+i]<<'\n';
            X=X+m1.array1d[M*main+i]*Problem.H[i+1];
        }
        
        eig(Er,V,X);

        for (int j=0; j<N; j++) {
            Ec.array1d[j]=exp(Complex(INTEGDT,0)*Complex(0,-1)*Er.array1d[j]);   //check to see if the Er number are always real?? It should be
        }
        
        Utot=(V*diagmat(Ec)*transpose(V))*Utot;
    }
    
    double val = abs(sum_mat(innerproduct(project_null(Utot),(Problem.U[0]))));
    double err= 1-val/8;
    return err;
}




cmatrix project_null(cmatrix& t1)
{
    int N=8;
    int G[]={0,1,4,5,10,11,13,14};
    cmatrix UTp(N,N);

    
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            UTp(i,j)=t1(G[i],G[j]);
        }
    }

    
    Complex  g(0,-1);
    double teta1=arg(UTp(4,4));
    double teta2=arg(UTp(2,2));
    double teta3=arg(UTp(1,1));
    double teta4=arg(UTp(0,0));
    cmatrix V_tmp(N,N);V_tmp.eye();
    V_tmp(0,0)=Complex(1,0);
    V_tmp(1,1)=exp(g*teta3);
    V_tmp(2,2)=exp(g*teta2);
    V_tmp(3,3)=exp(g*(teta2+teta3));
    V_tmp(4,4)=exp(g*teta1);
    V_tmp(5,5)=exp(g*(teta1+teta3));
    V_tmp(6,6)=exp(g*(teta1+teta2));
    V_tmp(7,7)=exp(g*(teta1+teta2+teta3));
    V_tmp=exp(g*teta4)*V_tmp;
    UTp=UTp*V_tmp;
    return UTp;
}


