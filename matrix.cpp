#include "matrix.h"
#include <iostream>
/*
 version 0.00-- Date=Sep 25 2014
 Written by Ehsan Zahedinejad. PhD at U of Calgary. Department of Physics...
 This is the matrix class implementation. A short explanation is given for each function
 in the matrix.h file. Note that I am treating 1d aray as 2d array in my program in a column major routine.
  The main reason is that you can not pass 2d array in MPI library between nodes so this could be a good reason to have 1d and deal as 2d arrays.
 */
//=============================================================================//

matrix::matrix(unsigned int t1,unsigned int t2)
{
    rows=t1;
    cols=t2;
    
    array1d= new double[rows*cols];
    
}

//=============================================================================//
matrix::matrix(const matrix &t)
{
    rows=t.rows;
    cols=t.cols;
    int size=rows*cols;
    int one=1;
    array1d= new double[size];
    dcopy(&size,t.array1d,&one,array1d,&one);
}

//=============================================================================//

matrix::~matrix()
{
    delete [] this->array1d;
}

//=============================================================================//
matrix& matrix::operator=(const matrix& t)
{
    if (this!=&t) {  //be aware of self-assignement
        this->~matrix();
        this->rows=t.rows;
        this->cols=t.cols;
        int size=rows*cols;
        int one=1;
        this->array1d= new double[size];
        dcopy(&size,t.array1d,&one,this->array1d,&one);

       // std::copy(t.array1d,t.array1d+(t.rows*t.cols),this->array1d);
    }
    return *this;
}

//=============================================================================//
// I am chaning this to a column major one
double& matrix::operator()(unsigned int t1,unsigned int t2){
    unsigned int r=this->rows;
    return(this->array1d[r*t2+t1]);
   /* return(this->array1d[t1*c+t2]); ROW MAJOR*/
}

//=============================================================================//
matrix operator+(const matrix& m1,const matrix& m2)
{
    unsigned int size=m1.rows*m1.cols;
    matrix temp(m1.rows,m1.cols);
    for (int i=0; i<size; ++i) {
            temp.array1d[i]=m1.array1d[i]+m2.array1d[i];
    }
    return temp;
}

//=============================================================================//
matrix operator*(const matrix& t1,const matrix& t2)
{
    int m=t1.rows,n=m,k=m,lda=m,ka=m,ldc=m,ldb=m;
    double alpha=1.00,beta=0; //C<--alpha*A*B_beta*C
    
    matrix product(m,m);
    /*
    Since we have row-major array and blas does the product in column major we first tranpose each matrix and then
    pass them to the dgemm function. TSD and TSA with having the fist letter "T" sending this message to the dgemm. The other
    letters "S" and "A" do not matter */
    dgemm("NSA","NSA",&m,&n,&k,&alpha,t1.array1d,&lda,t2.array1d,&ldb,&beta,product.array1d,&ldc);
    return product;
}

//=============================================================================//
matrix operator*(const double c,const matrix& m)
{
    matrix product(m.rows,m.cols);
    for (int i=0; i<m.rows*m.cols; i++)
        product.array1d[i]=m.array1d[i]*c;

    return product;
}

//=============================================================================//
matrix operator*(const matrix& m,const double c)
{
    return (c*m);
}

//=============================================================================//
matrix matrix::rowi(unsigned int i)
{
    int c=this->cols;
    int r=this->rows;
    matrix tmp(1,c);
    
    for (int j=0; j<c; j++) {
        tmp.array1d[j]=this->array1d[r*j+i];
    }
    
    
    
    return tmp;
}

//=============================================================================//
void matrix::eye()
{
    int r=this->rows;
    int c=this->cols;
    
    for (int i=0 ; i<r*r; i++) {
        this->array1d[i]=0;
    }
    
    for (int j=0; j<r; j++) {
        this->array1d[c*j+j]=1.0;
    }
}
//=============================================================================//
void eig(matrix& eigevalue,matrix& eigenvector,matrix& t1)
{
    int n=t1.rows;
    int lda=n;
    double vl=0,vu=0;
    int il=1, iu=n;
    double abstol=1e-14;
    int ldz=n;
    int lwork=40*n;
    int liwork=20*n;
    int* iwork=new int[liwork];
    double* work=new double[lwork];
    int m;
    int info;
    int* isuppz=new int[2*n];
    
        dsyevr("VECTOR","ALL","LOWER",&n,t1.array1d,&lda,&vl,&vu,&il,&iu,&abstol,&m,eigevalue.array1d,eigenvector.array1d,&ldz,isuppz,work,&lwork,iwork,&liwork,&info);
    delete [] iwork;
    delete [] work;
    delete [] isuppz;
}

//=============================================================================//
matrix conjugate(const matrix& t1) {
    
    //first making the identity matrix
    int rows=t1.rows;
    int cols=t1.cols;
    matrix eye(rows,cols);
    for (int i=0 ; i<rows*rows; i++) {
        eye.array1d[i]=0;
    }
    
    for (int j=0; j<rows; j++) {
        eye.array1d[cols*j+j]=1;
    }
    
    matrix product(t1.rows,t1.cols);
    
    int m=t1.rows,n=m,k=m,lda=m,ka=m,ldc=m,ldb=m;
    double alpha=1,beta=0; //C<--alpha*A*B_beta*C
    dgemm("CSA","NSA",&m,&n,&k,&alpha,t1.array1d,&lda,eye.array1d,&ldb,&beta,product.array1d,&ldc);
    
    return product;
}

//=============================================================================//
matrix transpose(const matrix& t1) {
    
    //first making the identity matrix
    int rows=t1.rows;
    int cols=t1.cols;
    matrix eye(rows,cols);
    for (int i=0 ; i<rows*rows; i++) {
        eye.array1d[i]=0;
    }
    
    for (int j=0; j<rows; j++) {
        eye.array1d[cols*j+j]=1;
    }
    
    matrix product(t1.rows,t1.cols);
    int m=t1.rows,n=m,k=m,lda=m,ka=m,ldc=m,ldb=m;
    double alpha=1,beta=0; //C<--alpha*A*B_beta*C
    dgemm("TSA","NSA",&m,&n,&k,&alpha,t1.array1d,&lda,eye.array1d,&ldb,&beta,product.array1d,&ldc);
    return product;
}
//=============================================================================//
matrix matrix::min()
{
    int c=this->cols;
    int r=this->rows;
    
    matrix id_min(1,2);
    id_min.array1d[0]=this->array1d[0];
    id_min.array1d[1]=0;
    
    for (int i=0; i<r*c; i++) {
        if (this->array1d[i]<id_min.array1d[0]) {
            id_min.array1d[0]=this->array1d[i];
            id_min.array1d[1]=i;
        }
    }
    
    return id_min;
}

//=============================================================================//
matrix matrix::operator+(double a) {
    
    int c=this->cols;
    int r=this->rows;
    
    matrix tmp(c,r);

    for (int i=0; i<r*c; i++)
        tmp.array1d[i]=this->array1d[i]+a;
    
    return tmp;
    
}

//=============================================================================//
matrix operator+(double a,matrix& t1) {
    
    return(t1+a);
    
}

//=============================================================================//
matrix diagmat(matrix& t1) {
    
    int cols=t1.cols;
    matrix tmp(cols,cols);
    for (int i=0 ; i<cols*cols; i++) {
        tmp.array1d[i]=0;
    }
    
    for (int j=0; j<cols; j++) {
        tmp.array1d[cols*j+j]=t1.array1d[j];
    }
    return tmp;
}







