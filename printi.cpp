#include "printi.h"

void printi(matrix t,string cfilename) {
    
    const char *filename = cfilename.c_str();
    std::ofstream file;
    int r=t.rows;
    int c=t.cols;
    
    if (!(file.is_open())){
        file.open (filename, ios::out | ios::trunc);
    }

    
    
//    file <<"rows "<<t.rows<<" columns "<<t.cols<<'\n';
//    file<<"======================================================================"<<'\n';
    
    
    for (int i=0; i<r; i++) {
        for (int j=0; j<c; j++) {
            file<<t(i,j)<<"  ";
        }
        file<<'\n';
    }

    
    file.close();
}


void printi(cmatrix t,string cfilename) {
    
    const char *filename = cfilename.c_str();
    std::ofstream file;
    int r=t.rows;
    int c=t.cols;
    
    if (!(file.is_open())){
        file.open (filename, ios::out | ios::trunc);
    }
    
    
    
    //    file <<"rows "<<t.rows<<" columns "<<t.cols<<'\n';
    //    file<<"======================================================================"<<'\n';
    
    
    for (int i=0; i<r; i++) {
        for (int j=0; j<c; j++) {
            file<<t(i,j)<<"  ";
        }
        file<<'\n';
    }
    
    
    file.close();
}





void prints(matrix t)
{
    
    for (int i=0; i<t.rows; i++) {
        for (int j=0; j<t.cols; j++) {
            std::cout<<t(i,j)<<"  ";
        }
        std::cout<<'\n';
    }
    
    
}


void printcs(cmatrix t)
{
    
    for (int i=0; i<t.rows; i++) {
        for (int j=0; j<t.cols; j++) {
            std::cout<<t(i,j)<<"  ";
        }
        std::cout<<'\n';
    }
    
    
}