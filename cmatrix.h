#ifndef COMPLEX_MATRIX
#define COMPLEX_MATRIX
#include <complex>
#include "matrix.h"
#include <iostream>

typedef std::complex<double> Complex;

class cmatrix
{
public:
    unsigned int rows,cols;
    std::complex<double>* array1d;
    
    cmatrix(unsigned int t1,unsigned int t2);               //Complex matrix constructor
    ~cmatrix();                                             //array destructor
    cmatrix(const cmatrix &t);                              //Complex copy constructor
    cmatrix& operator=(const cmatrix& t);                   //Complex copy assignment
//    double norm1();
    std::complex<double>& operator()(unsigned int t1,unsigned int t2);
    matrix cmatrix_to_matrix();
    void eye();
};
cmatrix operator+(const cmatrix& m1,const cmatrix& m2);     //Matrix+  overload
cmatrix operator+(const matrix& m1,const cmatrix& m2);      //Matrix+  overload
cmatrix operator+(const cmatrix& m1,const matrix& m2);      //Matrix+  overload
cmatrix operator*(Complex beta,const cmatrix t2);

cmatrix operator*(const cmatrix& t1,const cmatrix& t2);        //Matrix* overload
cmatrix operator*(const cmatrix& t1,const matrix& t2);      //Matrix* overload
cmatrix operator*(const matrix& t1,const cmatrix& t2);      //Matrix* overload


cmatrix operator*(const double c,const cmatrix& m);           //Matrix* overload
cmatrix operator*(const cmatrix& m,const double c);           //Matrix* overload
//double max(double* p,int i);
cmatrix diagmat(cmatrix& t1);                                                 //Making matrix diagonal
cmatrix transpose(const cmatrix& t1);                         //Transpose of a Matrix
std::complex<double> sum_mat(cmatrix t1);
cmatrix innerproduct(cmatrix t1,matrix t2);
#endif //#_COMPLEX_MATRIX

