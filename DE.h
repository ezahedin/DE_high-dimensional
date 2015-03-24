#ifndef DE_P
#define DE_P

//DE takes the problem as the input.
//problem include all the information that DE needs for optimization
#include "problem.h"
#include "algorithm"
#include "utilities.h"
#include "EvalUT.h"
#include "printi.h"
#include <string>
//void constrain(matrix* pop,float min,float max,unsigned int popsize,unsigned int vecsize);
//unsigned int uniquerand(trng::yarn2 *R, unsigned int i1,unsigned int i2);
void DE(problem& Problem,int mynode);
double DE(problem& Problem,matrix& _pop,matrix& single_pop,int mynode,int* seed,matrix& Fidelity,matrix& rates,int Gen);
#endif