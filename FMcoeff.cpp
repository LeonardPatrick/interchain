#include<cmath>
#include<complex>

using namespace std;

complex<double> im = complex<double>(0.0,1.0);                                                                                                            
double Hgfac = -2.00231930436;                                                                       

double Hmu = 5.7883818066e-5 * 1000;                // meV T^-1                               





//Sx_Q|GS>


complex<double> ax(double Jy, double Jx, double Jz, int j, double C[])
{



  return 0.5*( exp(im*C[j]) - ( 1.0 + exp(im*2.0*C[j]) )*( (Jy - Jx)/(4.0*Jz) ) );


}


complex<double> bx(double Jz, double B,  int j, double C[])
{


 return -0.5 * ((Hgfac*Hmu*B)/(2.0*Jz)) * ( exp(im*C[j]) + exp(im*2.0*C[j]) );

}






 complex<double> cx(double Jy, double Jx, double Jz, int j, double C[])
{



  return -0.5 * ( (Jy - Jx)/(4.0*Jz) ) * ( exp(im*C[j]) + exp(im*3.0*C[j]) );

}










//Sy_Q|GS>



 complex<double> ay(double Jy, double Jx, double Jz, int j, double C[])
{



 return  -0.5*im*( -exp(im*C[j]) + ( -(Jy - Jx)/(4.0*Jz) ) * ( 1.0 + exp(im*2.0*C[j]) ) );

}





complex<double> by( double Jz, double B,  int j, double C[])
{



return  -0.5*im*( (-Hgfac*Hmu*B)/(2.0*Jz) )*(-1.0)*( exp(im*C[j])  +   exp(im*2.0*C[j]) ); 


}





complex<double> cy(double Jy, double Jx, double Jz, int j, double C[])
{




 return  -0.5*im*( -(Jy - Jx)/(4.0*Jz) ) * (-1.0) * ( exp(im*C[j]) + exp(im*3.0*C[j]) );


}








//Sz_Q|GS>


complex<double> az( double Jz, double B,  int j, double C[])
{



  return  -( (-Hgfac*Hmu*B)/(2.0*Jz) ) * exp(im*C[j]);


}





complex<double> bz(double Jy, double Jx, double Jz,  int j, double C[])
{

 


 return  ( -(Jy -  Jx)/(4.0*Jz) ) * (-1.0) * ( exp(im*C[j]) + exp(im*2.0*C[j]) );


}








