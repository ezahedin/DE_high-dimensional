#include "cmatrix.h"

/*
 version 0.00-- Date=Sep 25 2014
 Written by Ehsan Zahedinejad. PhD at U of Calgary. Department of Physics...
 This is the std::complex<double> matrix class implementation. A short explanation is given for each function
 in the matrix.h file. Note that I am treating 1d aray as 2d array in my program.
 This is just to increase the cashe waste and memory miss and to spped up the algorithms
 performance. The other reason is that you can not pass 2d array in MPI library between nodes se this could be a good reason to have 1d and deal as 2d arrays.
 */
//=============================================================================//
cmatrix::cmatrix(unsigned int t1,unsigned int t2) {
    rows=t1;
    cols=t2;
    array1d=new std::complex<double>[rows*cols];
}

//=============================================================================//

cmatrix::~cmatrix()
{
    delete [] array1d;
}

//=============================================================================//
cmatrix::cmatrix(const cmatrix &t)
{
    rows=t.rows;
    cols=t.cols;
    int size=rows*cols;
    int one=1;
    array1d= new std::complex<double>[size];
    zcopy(&size,t.array1d,&one,array1d,&one);
}

//=============================================================================//
cmatrix& cmatrix::operator=(const cmatrix& t) {
    if (this!=&t) {
        this->~cmatrix();
        rows=t.rows;
        cols=t.cols;
        int size=rows*cols;
        int one=1;
        array1d= new std::complex<double>[rows*cols];
        zcopy(&size,t.array1d,&one,array1d,&one);

}
    return *this;
}
//=============================================================================//
cmatrix operator+(const cmatrix& m1,const cmatrix& m2) {
    
    int r=m1.rows;
    int c=m1.cols;
    cmatrix sum(r,c);
    for (int i=0; i<r*c; i++)
        sum.array1d[i]=m1.array1d[i]+m2.array1d[i];
    
    return sum;
    
}

//=============================================================================//
cmatrix operator+(const matrix& m1,const cmatrix& m2)
{
    int r=m1.rows;
    int c=m1.cols;
    cmatrix sum(r,c);
    for (int i=0; i<r*c; i++)
        sum.array1d[i]=std::complex<double>(m1.array1d[i])+m2.array1d[i];
    
    return sum;
}


//=============================================================================//
cmatrix operator+(const cmatrix& m1,const matrix& m2) {
    
    return(m2+m1);
    
}

//=============================================================================//
cmatrix operator*(const cmatrix& t1,const cmatrix& t2) {
    cmatrix product(t1.rows,t2.cols);
    
    int m=t1.rows,n=m,k=m,lda=m,ka=m,ldc=m,ldb=m;
    Complex alpha=Complex(1.0),beta=Complex(0); //C<--alpha*A*B_beta*C
    zgemm("NSA","NSA",&m,&n,&k,&alpha,t1.array1d,&lda,t2.array1d,&ldb,&beta,product.array1d,&ldc);
    
    return product;
}

//=============================================================================//
cmatrix operator*(const matrix& t1,const cmatrix& t2) {
    
    int rows=t1.rows;
    int cols=t1.cols;
    int lda=rows;
    int ldc=rows;
    int ldb=rows;
    cmatrix product(rows,cols);
    double* rwork=new double[2*rows*cols];

    zlarcm(&rows,&cols,t1.array1d,&lda,t2.array1d,&ldb,product.array1d,&ldc,rwork);
    delete [] rwork;
    return product;
}

//=============================================================================//
cmatrix operator*(const cmatrix& t1,const matrix& t2) {
    
    int rows=t1.rows;
    int cols=t1.cols;
    int lda=rows;
    int ldc=rows;
    int ldb=rows;
    cmatrix product(rows,cols);
    double* rwork=new double[2*rows*cols];
    
    zlacrm(&rows,&cols,t1.array1d,&lda,t2.array1d,&ldb,product.array1d,&ldc,rwork);
    delete [] rwork;
    
    return product;
    
}

//=============================================================================//
cmatrix operator*(const double c,const cmatrix& m) {
    cmatrix product(m.rows,m.cols);
    
    for (int i=0; i<m.rows*m.cols; i++) {
        product.array1d[i]=std::complex<double>(c)*m.array1d[i];
    }
    
    return product;
    
}

//=============================================================================//
cmatrix operator*(const cmatrix& m,const double c)
{
    return(c*m);
}


//=============================================================================//
std::complex<double>& cmatrix::operator()(unsigned int t1,unsigned int t2)
{
    unsigned int r=this->rows;
    return(this->array1d[r*t2+t1]);
}

//=============================================================================//
matrix cmatrix::cmatrix_to_matrix()
{
    int r=this->rows;
    int c=this->cols;

    matrix tmp(r,c);
    
    for (int i=0; i<r*c; i++)
        tmp.array1d[i]=real(this->array1d[i]);
    
    return tmp;
}

//=============================================================================//
void cmatrix::eye()
{
    int r=this->rows;
    for (int i=0 ; i<r*r; i++) {
        this->array1d[i]=Complex(0);
    }
    
    for (int j=0; j<r; j++) {
        this->array1d[cols*j+j]=Complex(1.0);
    }
}

//=============================================================================//
cmatrix diagmat(cmatrix& t1) {
    
    int cols=t1.cols;
    cmatrix tmp(cols,cols);
    for (int i=0 ; i<cols*cols; i++) {
        tmp.array1d[i]=Complex(0);
    }
    
    for (int j=0; j<cols; j++) {
        tmp.array1d[cols*j+j]=t1.array1d[j];
    }
    return tmp;
}

//=============================================================================//
cmatrix transpose(const cmatrix& t1) {
    
    //first making the identity matrix
    int rows=t1.rows;
    int cols=t1.cols;
    cmatrix eyec(rows,cols);
    for (int i=0 ; i<rows*rows; i++) {
        eyec.array1d[i]=Complex(0);
    }
    
    for (int j=0; j<rows; j++) {
        eyec.array1d[cols*j+j]=Complex(1.0);
    }

    
    cmatrix product(cols,rows);
    
    int n=rows,k=rows,lda=rows,ka=rows,ldc=rows,ldb=rows;
    Complex alpha=Complex(1.0),beta=Complex(0); //C<--alpha*A*B_beta*C
    zgemm("TSA","NSA",&rows,&n,&k,&alpha,t1.array1d,&lda,eyec.array1d,&ldb,&beta,product.array1d,&ldc);
    
    return product;
}

//=============================================================================//
std::complex<double>  sum_mat(cmatrix t1) {
    
    std::complex<double> tmp(0,0);
    int r=t1.rows;
    int c=t1.cols;
    for (int i=0; i<r*c; i++) {
        tmp+=t1.array1d[i];

    }
    return tmp;
}

//=============================================================================//
cmatrix operator*(Complex beta,cmatrix t2) {
    cmatrix product=t2;
    
    int m=t2.rows,n=m,k=m,lda=m,ka=m,ldc=m,ldb=m;
    Complex alpha=Complex(0); //C<--alpha*A*B_beta*C
    zgemm("NSA","NSA",&m,&n,&k,&alpha,product.array1d,&lda,product.array1d,&ldb,&beta,product.array1d,&ldc);
    
    return product;
}

//=============================================================================//
cmatrix innerproduct(cmatrix t1,matrix t2) {
    
    int r=t1.rows;
    int c=t1.cols;
    cmatrix tmp(r,c);
    
    
    for (int i=0; i<r; i++) {
        for (int j=0; j<c; j++) {
            tmp(i,j)=t1(i,j)*Complex(t2(i,j));
        }
    }
    
    return tmp;
    
}














