#ifndef UTILITIES_H
#define UTILITIES_H

#include "cmatrix.h"
#include <cstdlib>      // std::rand, std::srand
#include <time.h>
#include <iostream>
#include "problem.h"
#include <math.h>
#include "printi.h"
/* rand_matrix function explanation--This matrix generate random numbers according to the following definitions:
 a=type of distribution
    a=1: uniform (0,1)
    a=2: uniform (-1,1)
    a= 3: normal (0,1)
 seed--initial seed. its getting updated after each time we call this function
 matrix t1:matrix includes tha final results
 */
void rand_matrix(int a,int* seed,matrix& t1);


// random generator functions for index of the matrices:
void random_shuffle1(int* mat,int size,int mynode);

//constraining the elements of a matrix between a min and max value
void constrain(matrix& pop,double min,double max);

//Projecting from Hilbert space to Computational space
cmatrix project_null(cmatrix& t1);
matrix projection(matrix& t1);
matrix erfun(problem& Problem,matrix& f);

#endif

