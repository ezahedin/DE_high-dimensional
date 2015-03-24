#include "kron.h"
//===================================================================================//

cmatrix kron(cmatrix k1,cmatrix k2)
{
    unsigned int r1=k1.rows;
    unsigned int c1=k1.cols;
    unsigned int r2=k2.rows;
    unsigned int c2=k2.rows;
    cmatrix K(r1*r2,r1*r2);
    for (int i=0; i<r1; i++)
        for (int j=0; j<c1; j++)
            for (int k=0; k<r2; k++)
                for (int l=0; l<c2; l++)
                    K(k+i*c2,l+j*c2)=k1(i,j)*k2(k,l);
    return K;
}

//===================================================================================//
cmatrix kron(matrix k1,cmatrix k2)
{
    unsigned int r1=k1.rows;
    unsigned int c1=k1.cols;
    unsigned int r2=k2.rows;
    unsigned int c2=k2.rows;
    cmatrix K(r1*r2,r1*r2);
    for (int i=0; i<r1; i++)
        for (int j=0; j<c1; j++)
            for (int k=0; k<r2; k++)
                for (int l=0; l<c2; l++)
                    K(k+i*c2,l+j*c2)=Complex(k1(i,j))*k2(k,l);
    
    return K;
}

//===================================================================================//
cmatrix kron(cmatrix k1,matrix k2)
{
    unsigned int r1=k1.rows;
    unsigned int c1=k1.cols;
    unsigned int r2=k2.rows;
    unsigned int c2=k2.rows;
    cmatrix K(r1*r2,r1*r2);

    for (int i=0; i<r1; i++)
        for (int j=0; j<c1; j++)
            for (int k=0; k<r2; k++)
                for (int l=0; l<c2; l++)
                    K(k+i*c2,l+j*c2)=k1(i,j)*Complex(k2(k,l));
    
    return K;
}
//===================================================================================//

matrix kron(matrix k1,matrix k2)
{
    unsigned int r1=k1.rows;
    unsigned int c1=k1.cols;
    unsigned int r2=k2.rows;
    unsigned int c2=k2.rows;
    matrix K(r1*r2,r1*r2);

    for (int i=0; i<r1; i++)
        for (int j=0; j<c1; j++)
            for (int k=0; k<r2; k++)
                for (int l=0; l<c2; l++)
                    K(k+i*c2,l+j*c2)=k1(i,j)*k2(k,l);

    return K;
}
//===================================================================================//



