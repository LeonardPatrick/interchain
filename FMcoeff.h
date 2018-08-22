#ifndef FMcoeff_H
#define FMcoeff_H



#include<complex>

using namespace std;


// In the following j refers to the index for the array C, which represents the QQ values. Jy, Jx and Jz refer to the exchange energies and B the magnitude of the transverse magnetic field.



//Sx_Q|GS>

complex<double> ax(double Jy, double Jx, double Jz,  int j, double C[]);

complex<double> bx(double Jz, double B, int j, double C[]);

complex<double> cx(double Jy, double Jx, double Jz, int j, double C[]);



//Sy_Q|GS>

complex<double> ay(double Jy, double Jx, double Jz, int j, double C[]);

complex<double> by(double Jz, double B, int j, double C[]);

complex<double> cy(double Jy, double Jx, double Jz, int j, double C[]);





//Sz_Q|GS>

complex<double> az( double Jz, double B,  int j, double C[]);

complex<double> bz(double Jy, double Jx, double Jz, int j, double C[]);




#endif

