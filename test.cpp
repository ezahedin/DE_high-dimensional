/* frexp example
#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

int main ()
{
    int i=0;
    for (int i=0; i<10; ++i) {
        cout<<i<<'\n';
    }
    cout<<i<<'\n';

    return 0;
}

*/


#include <iostream>
#include <Eigen/Dense>
#include "cmatrix.h"
#include "matrix.h"
#include <complex>
using namespace std;
using namespace Eigen;

void testfun1(double** G)
{
    cout<<G[0][0]<<'\n';
}


int main()
{
    cmatrix T(3,3);
    T.carray2d[0][0]=complex<double>(10,2);
    T.carray2d[0][1]=complex<double>(2,3);
    T.carray2d[0][2]=complex<double>(1,2);
    T.carray2d[1][0]=complex<double>(2,2);
    T.carray2d[1][1]=complex<double>(2,2);
    T.carray2d[1][2]=complex<double>(20,3);
    T.carray2d[2][0]=complex<double>(2,2);
    T.carray2d[2][1]=complex<double>(8,2);
    T.carray2d[2][2]=complex<double>(2,2);
    matrix T1(3,3);
        
    int m=3;
    cmatrix *mat[m];
    mat[2]=new cmatrix(3,3);
    *mat[2]=T;
    *mat[2]=*mat[2]*(*mat[2]);
   // cout<<mat[2]->carray2d[0][0]<<'\n';
    
    MatrixXcd* P;
    P=new MatrixXcd(3,3);
    
    MatrixXcd A(3,3),B(3,3);
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            A(i,j)=T.carray2d[i][j];
        }
    }
    //A<< complex<double>(1,-1), complex<double>(1,2), complex<double>(2,1),
    //complex<double>(13,10), complex<double>(11,5), complex<double>(1,9),
    //complex<double>(1.1,15), complex<double>(6,100), complex<double>(8,9);
    Vector3cd b;
    //A << 1,2,3,  4,5,6,  7,8,10;
    b << complex<double>(1,-12), complex<double>(1.1,1), complex<double>(1.1,3);
    //B=A;
    //double C[3][3]={ { -1, 9, 12}, { 12, 4, 5 }, { 1, 2, 1} };
    Vector3cd x=A.colPivHouseholderQr().solve(A.col(0));
    
    T.carray2d[0][0]=x(0);
   // cout << "T forst element is" <<T.carray2d[0][0]<< endl;
   // cout << "The solution is:\n" <<x<< endl;
    
    cmatrix* Apowers=new cmatrix[10];
    Apowers[0]=T;
    Apowers[1]=T;
    //Apowers[0].~cmatrix();
    
    cout<<Apowers[0].carray2d[1][2]<<'\n';
    cout<<Apowers[1].carray2d[1][2]<<'\n';

    int *p=new int[10];
    cout<<p[0]<<'\n';
    
    
    delete[] p;
    delete[] Apowers;
    //Apowers=NULL;
 
    cout<<p[0]<<'\n';
    cout<<T.carray2d[1][2]<<'\n';

    cout<<Apowers[0].carray2d[0][1]<<'\n';
    cout<<Apowers[1].carray2d[0][1]<<'\n';

    
    //for (int i=0; i<10; i++) {
    //    Apowers[i].~cmatrix();
    //}
    
        //T1.array2d[2][2]=10;
    //testfun1(T1.array2d);
    
    
    return 0;
}

