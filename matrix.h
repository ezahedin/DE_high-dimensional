#ifndef MATRIX
#define MATRIX
#include <algorithm>                                        // std::copy
#include <complex>

#define MKL_Complex16 std::complex<double>

#include "mkl.h"

/*Function for Matric<double>*Matrix<double> multiplication "Blas routine"*/
//extern "C" void dgemm(char* TransA,char* TransB,int* m,int* n,int* k,double* alpha,double* a,int* lda,double* b,int* ldb,double* beta,double* c,int* ldc);




class matrix
{
public:
    unsigned int rows,cols;
    double* array1d;
    matrix(unsigned int t1,unsigned int t2);                //Constructor
    matrix(const matrix &t);                                //Matrix copy constructor
    matrix& operator=(const matrix& t);                     //Matrix copy assignment
//    matrix reshape(unsigned int t1,unsigned int t2);        //Reshape Matrix
    ~matrix();                                              //deallocate memroy
    double& operator()(unsigned int t1,unsigned int t2);    //Forming Mat(i,j) instead of Mat.array[i][j]
    matrix rowi(unsigned int i);                            //return the ith row of matrix
    void eye();                                             //Making Identity Matrices
    matrix min();
    matrix operator+(double a);
};
matrix operator+(double a,matrix& t1);
matrix operator+(const matrix& m1,const matrix& m2);        //Matrix+  overload
matrix operator*(const matrix& t1,const matrix& t2);        //Matrix* overload
matrix operator*(const double c,const matrix& m);           //Matrix* overload
matrix operator*(const matrix& m,const double c);           //Matrix* overload
void eig(matrix& eigevalue,matrix& eigenvector,matrix& t1); //Calculate the eigenvalue and eigenvector of a symetric matrix
matrix transpose(const matrix& t1);
matrix conjugate(const matrix& t1);
matrix diagmat(matrix& t1);
#endif






