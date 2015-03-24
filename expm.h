#include "cmatrix.h"
#include "matrix.h"
#include "math.h"
#include "complex.h"
#include "algorithm"
#include <iostream>
#include "eigen/Dense"  //Eigen library to exploit the matrix left devision



cmatrix expm(cmatrix& A);
cmatrix PadeApproximantOfDegree(int m,cmatrix& A);
//int ceil(double m);

matrix getPadeCoefficients(int m);

