#ifndef COMPLEXh_H
#define COMPLEXh_H


class Complex {
    
private:
    float real,imag;
    
public:
    
    Complex(float real=0,float imag=0); // Complex Constructor
    Complex(float real) {imag=0;};      // Complex Constructor when imaganry part is zero
    Complex(const Complex &t);       //Complex copy constructor
    Complex& operator=(const Complex& t); //Complex copy assignment
    Complex& operator+=(Complex t1);  //Complex sum operator
    float getreal();
    float getimag();
    
};

Complex operator+(Complex t1,Complex t2);

#endif /* COMPLEXh_H */