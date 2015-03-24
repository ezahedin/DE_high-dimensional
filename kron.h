#ifndef KRONKER
#define KRONKER
#include "cmatrix.h"
#include "matrix.h"

cmatrix kron(cmatrix k1,cmatrix k2);
cmatrix kron(matrix k1,cmatrix k2);
cmatrix kron(cmatrix k1,matrix k2);
matrix  kron(matrix k1,matrix k2);

#endif