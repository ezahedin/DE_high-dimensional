#include "expm.h"
using namespace Eigen;
using namespace std;

cmatrix expm(cmatrix& A)
{
    cmatrix F(A.row,A.row);
    
    
    //Check the class of input A and initialize M_VALS and THETA accordingly.
    int m_vals[5] = { 3, 5, 7, 9, 13};
    double theta[5]={0.01495585217958292,0.2539398330063230,0.9504178996162932,2.097847961257068,5.371920351148152};
    double normA=A.norm1();
    if (normA<=theta[4]){
        // no scaling and squaring is required.
        
        for (int i=0; i<5; i++) {
            if (normA<=theta[i]) {
                F=PadeApproximantOfDegree(m_vals[i],A); //I need to define this matrix lol F is a cmatrix type
                break;
            }
        }
        
        
    } else {
        int n;
        double t;
        int s;
        t = frexp(normA/theta[4],&s);
        s = s - (t == 0.5); //adjust s if normA/theta(end) is a power of 2.
        A=A*(1/pow(2,s)); //Scaling
        
        F = PadeApproximantOfDegree(m_vals[4],A);
        for (int i=0; i<s; i++) {
            F=F*F; //Squaring
        }
    }
    
    
    return F;
    
    
}


//This is an expensive function and a good target for GPU part.
cmatrix PadeApproximantOfDegree(int m,cmatrix& A)
{
    
    int n=A.row;
    matrix c=getPadeCoefficients(m); //have to define this function
    cmatrix I(n,n);
    cmatrix U(n,n),V(n,n),F(n,n);
   
    
    
    //cmatrix Apowers[(int)ceil((m+1)/2)];
    cmatrix *Apowers=new cmatrix[(int)ceil((m+1)/2)];
    //for (int j=0; j<(int)ceil((m+1)/2); j++) {
    //    Apowers[j]=new cmatrix(n,n);
    //}
    //Apowers=new cmatrix[(int)ceil((m+1)/2)];
    MatrixXcd _A(n,n),B(n,n);  //This type is from Eigen lib
    VectorXcd x(n);                //This type is from Eigen lib. This might lead to error "Segmentation error". Be hold when you test the code
    for (int i=0; i<n; i++) {
        I.carray2d[i][i]=std::complex<double>(1,0);
    }
    
    switch (m) {
        case 3: case 5: case 7: case 9:
            Apowers[0]=I;
            Apowers[1]=A*A;
            
            if (m==3) {
                for (int i=2; i>0; i--) {
                    Apowers[i]=(Apowers[i-1])*(Apowers[1]);
                }
            }
            else {
                for (int i=2; i<((int)ceil((m+1))/2)+1; i++) {
                    Apowers[i]=(Apowers[i-1])*(Apowers[1]);
                }
            }
            
            for (int i=m; i>0; i-=2) {
                U=U+c.array2d[0][i]*(Apowers[(i+1)/2-1]);
            }
                U=A*U;
            
                for (int j=m-1 ; j>-1; j-=2) {
                    V=V+c.array2d[0][j]*(Apowers[(j/2)]);
                }
            break;
            case 13:
        {
                cmatrix A2=A*A;cmatrix A4=A2*A2; cmatrix A6=A2*A4;
            
                U = A * (A6*(c.array2d[0][13]*A6 + c.array2d[0][11]*A4 + c.array2d[0][9]*A2) + c.array2d[0][7]*A6 + c.array2d[0][5]*A4 + c.array2d[0][3]*A2 + c.array2d[0][1]*I);
                V = A6*(c.array2d[0][12]*A6 + c.array2d[0][10]*A4 + c.array2d[0][8]*A2) + c.array2d[0][6]*A6 + c.array2d[0][4]*A4 + c.array2d[0][2]*A2 + c.array2d[0][0]*I;
        }
            break;
            
        default:
            break;
            }
            //This part is tricky. Here we need to solve a system of linear equations
            // I'll do this part in a type-casting method in the future. So we can speed up
            for (int i=0; i<n; i++) {
                for (int j=0; j<n; j++) {
                    _A(i,j)=V.carray2d[i][j]-U.carray2d[i][j];
                    B(i,j)=std::complex<double>(2,0)*U.carray2d[i][j];
                }
            }
            
            for (int i=0; i<n; i++) {
                x=_A.colPivHouseholderQr().solve(B.col(i));
                for (int j=0; j<n; j++) {
                    F.carray2d[i][j]=x(j);
                }
                
            }
    
    
    delete[] Apowers;
             return (F = F + I); // (-U+V)\(U+V);
    
    }

    
matrix getPadeCoefficients(int m)
    {
    // GETPADECOEFFICIENTS Coefficients of numerator P of Pade approximant
    //C = GETPADECOEFFICIENTS(M) returns coefficients of numerator
    //    of [M/M] Pade approximant, where M = 3,5,7,9,13.
        matrix c(1,14);
        switch (m) {
            case 3:
                c.array2d[0][0]=120;
                c.array2d[0][1]=60;
                c.array2d[0][2]=12;
                c.array2d[0][3]=1;
                break;
            case 5:
                c.array2d[0][0]=30240;
                c.array2d[0][1]=15120;
                c.array2d[0][2]=3360;
                c.array2d[0][3]=420;
                c.array2d[0][4]=30;
                c.array2d[0][5]=1;

                break;
            case 7:
                c.array2d[0][0]=17297280;
                c.array2d[0][1]=8648640;
                c.array2d[0][2]=1995840;
                c.array2d[0][3]=277200;
                c.array2d[0][4]=25200;
                c.array2d[0][5]=1512;
                c.array2d[0][6]=56;
                c.array2d[0][7]=1;
                
                break;
                
            case 9:
                c.array2d[0][0]=17643225600;
                c.array2d[0][1]=8821612800;
                c.array2d[0][2]=2075673600;
                c.array2d[0][3]=302702400;
                c.array2d[0][4]=30270240;
                c.array2d[0][5]=2162160;
                c.array2d[0][6]=110880;
                c.array2d[0][7]=3960;
                c.array2d[0][8]=90;
                c.array2d[0][9]=1;
                
                break;
            case 13:
                c.array2d[0][0]=64764752532480000;
                c.array2d[0][1]=32382376266240000;
                c.array2d[0][2]=7771770303897600;
                c.array2d[0][3]=1187353796428800;
                c.array2d[0][4]=129060195264000;
                c.array2d[0][5]=10559470521600;
                c.array2d[0][6]=670442572800;
                c.array2d[0][7]=33522128640;
                c.array2d[0][8]=1323241920;
                c.array2d[0][9]=40840800;
                c.array2d[0][10]=960960;
                c.array2d[0][11]=16380;
                c.array2d[0][12]=182;
                c.array2d[0][13]=1;
                
                break;
            default:
                break;
        }
        return c;
    }
    
   
    
    
/*
int ceil(int m)
{
    if (((int)m==m)||m<0) {
        return (int)m;
    } else  {
        return ((int)m+1);
    }
    
}
 */









