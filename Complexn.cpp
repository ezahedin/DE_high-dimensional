#include "Complexn.h"

//Constructor
Complex::Complex(float t1,float t2)
{
    real=t1;
    imag=t2;
}



//Complex copy constructor
Complex::Complex(const Complex &t)
{
    real=t.real;
    imag=t.imag;
}

//Complex assignment copy
Complex& Complex::operator=(const Complex& t)
{
    if (this!=&t) {  //be aware of self-assignement
        real=t.real;
        imag=t.imag;
    }
    return *this;
    
}
//overloading += operator for complex numebrs
inline Complex& Complex::operator+=(Complex t1)
{
    real+=t1.real;
    imag+=t1.imag;
    return *this;
}

//getter and setter for the complex number

float Complex::getreal()
{
    return real;
}

float Complex::getimag()
{
    return imag;
}
//overloading + operator for complex numebrs, this is not a member function
Complex operator+(Complex t1,Complex t2)
{
    Complex tmp=t1;
    return t1+=t2;
}



