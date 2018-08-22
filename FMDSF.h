#ifndef FMDSF_H
#define FMDSF_H



#include<complex>
#include"FMcoeff.h"      // to access the coefficients to permit full computation of DSF.


using namespace std;
using namespace arma;    // to use submatrix operation, submat.




// These are the function declarations for the dynamical structure factors (DSF) [p. 894 in notes]. The input parameters should be an array (A[]) (for Q_||), an index to mark each element in the array(j), the three exhange energies (Jy, Jx, Jz), the transverese B field (B), the Green function matrix (G), an index to mark the different submatrices(k). A double will be returned since the DSF are real (Syz - Szy is treated as real as it is purely imaginary but is multiplied by an i in the expression for the total cross section).


// NOTE: Jy is placed before Jx since (Jy - Jx) naturally arises in XYZ calculations.


//(A[], j, Jy, Jx, Jz, B, G, k)

double SxxH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k);

double SyyH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k);

double SzzH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k);

double SyzzyH(double A[], int j, double Jy, double Jx, double Jz, double B, cx_mat G, int k);










#endif
