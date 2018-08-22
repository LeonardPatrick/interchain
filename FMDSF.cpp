#include<cmath>
#include<complex>
#include "/Users/leonardenglish/Desktop/Thoughts on Life/Two Spin Soliton/Numerical/armadillo-3.900.6/include/armadillo"
#include <Accelerate/Accelerate.h>

#include"FMcoeff.h"  // accessing coefficients


using namespace std;
using namespace arma;


const double PI = 4.0*atan(1.0);
complex<double>iii = complex<double>(0.0,1.0);               // Defining the purely imaginary unit vector.                                                                                                 



double SxxH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k)
{

  // return imag( G.submat(4*k, 0, (4*k + 3), 3)(1,0));

  //  return -2.0 * (1.0/PI) * imag( conj(ax(Jy, Jx, Jz, j, A ) ) * ax(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,0) + conj( ax(Jy, Jx, Jz, j, A ) ) * bx(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,1)  + conj(ax(Jy, Jx, Jz, j, A ) ) * cx(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,2) + conj(bx(Jz, B, j, A )) * ax(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,0)  + conj(cx(Jy, Jx, Jz, j, A )) * ax(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,0) );


    return -2.0 * (1.0/PI) * imag( conj(ax(Jy, Jx, Jz, j, A ) ) * ax(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,0) + conj( ax(Jy, Jx, Jz, j, A ) ) * bx(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,1)  + conj(ax(Jy, Jx, Jz, j, A ) ) * cx(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,2) + conj(bx(Jz, B, j, A )) * ax(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,0) + conj(bx(Jz, B, j, A ) ) * bx(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,1) + conj(bx(Jz, B, j, A )) * cx(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,2) + conj(cx(Jy, Jx, Jz, j, A )) * ax(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,0) + conj(cx(Jy, Jx, Jz, j, A )) * bx(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,1) + conj(cx(Jy, Jx, Jz, j, A )) * bx(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,2)  );

}





double SyyH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k)
{

  return -2.0 * (1.0/PI) * imag( conj(ay(Jy, Jx, Jz, j, A ) ) * ay(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,0) + conj( ay(Jy, Jx, Jz, j, A ) ) * by(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,1)  + conj(ay(Jy, Jx, Jz, j, A ) ) * cy(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,2) + conj(by(Jz, B, j, A )) * ay(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,0) + conj(by(Jz, B, j, A ) ) * by(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,1) + conj(by(Jz, B, j, A )) * cy(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,2) + conj(cy(Jy, Jx, Jz, j, A )) * ay(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,0) + conj(cy(Jy, Jx, Jz, j, A )) * by(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,1) + conj(cy(Jy, Jx, Jz, j, A )) * by(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,2)  );

}






double SzzH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k)
{

  return -2.0 * (1.0/PI) * imag( conj(az(Jz, B, j, A ) ) * az(Jz, B , j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,0) + conj( az(Jz, B, j, A ) ) * bz(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,1) + conj(bz(Jy, Jx, Jz, j, A )) * az(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,0) + conj(bz(Jy, Jx, Jz, j, A ) ) * bz(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,1) );


}









double SyzzyH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k)
{

  // minus is in place at front to account for i*i = -1. Syz - Szy is purely imaginary and total cross section multiplies this term by i. 

  return -2.0 * (1.0/PI) * imag(   iii*(       (  conj(ay(Jy, Jx, Jz, j, A))  * az(Jz, B , j, A) - conj(az(Jz, B, j, A))  * ay(Jy, Jx, Jz, j, A)   ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,0)  +  ( conj(ay(Jy, Jx, Jz, j, A)) * bz(Jy, Jx, Jz, j, A) - conj(az(Jz, B, j,  A )) * by(Jz, B, j, A) ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,1)  +  (  conj(by(Jz, B, j, A )) * az(Jz, B, j, A) - conj(bz(Jy, Jx, Jz, j, A )) * ay(Jy, Jx, Jz, j, A) ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,0)  +  ( conj(by(Jz, B, j, A ) ) * bz(Jy, Jx, Jz, j, A ) - conj(bz(Jy, Jx, Jz, j, A ) ) * by(Jz, B, j, A ) ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,1)  +  conj(cy(Jy, Jx, Jz, j, A ) ) * az(Jz, B, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,0)  -  conj(az(Jz, B, j, A ) ) * cy(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(0,2)  +  conj(cy(Jy, Jx, Jz, j, A ) ) * bz(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(2,1)  -  conj(bz(Jy, Jx, Jz, j, A ) ) * cy(Jy, Jx, Jz, j, A ) * G.submat(4*k, 0,(4*k+ 3), 3)(1,2)   )  );


}

