#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>       
#include "/Users/leonardenglish/Desktop/Staggered Field/armadillo-3.900.6/include/armadillo"
#include <Accelerate/Accelerate.h>
#include<time.h>

#include "FMcoeff.h"    // To access coefficients. I think the ordering does matter here. 
#include "FMDSF.h"      // To access the dynamical structure factors.
#include "FMBS.h"    // Header file for magnitude of external magnetic field in the z direction, no. of data points (H,J) and how many times the for loop should be divided (Z). Z MUST THE SAME VALUE AS Z IN THE BASH SCRIPT FM Copy GF Stag.sh 


 
using namespace std;

using namespace arma;


const double PI = 4.0*atan(1.0);





cx_mat HSF(int p, complex<double> t, complex<double> tDag, complex<double> L,  complex<double> LDag, complex<double> g, complex<double> gDag, double QQ);





// Initial matrices




cx_mat a0(int p,  complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag, double EnM, double QQ);

cx_mat b0(int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag,  cx_mat d, cx_mat dDag, double EnM, double QQ);

cx_mat e0 (int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ); 

cx_mat es0(int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ); 




// Recurrence relation






cx_mat ak(int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ); 

cx_mat bk (int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ); 

cx_mat ek(int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ); 

cx_mat esk(int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ); 






// NOTE: h, in the case of the FM configuration, is g mu_B B_z. Bz is declared and defined in the header file.


double Jz;   

double Jx;
double Jy;

double B;
double DM;



complex<double>ii = complex<double>(0.0,1.0);               // Defining the purely imaginary unit vector.                                                                                               
double  gfac = -2.00231930436;  
double  mu = 5.7883818066e-5 * 1000;                       // meV T^-1



int pp;

double eta =  0.0001;                                          // This is just a constant                                                                                                               


mat id = eye<mat>(4,4);









int main()
{

  clock_t time1, time2;

  time1 = clock();

 


 Jz = 13.8;                                         // meV.
 Jx = 3.96986 ;
 Jy = 1.89041;
  B = 0.0;                                           // T 
 DM = 0.0;//(Jy-Jx)/7.0; //                                        // meV


 Bz = 2.2;   // Uniform (in space and time) magnetic field in Ising direction.



double P = 1.0;


double eta =  0.0001;                                          // This is just a constant                                                                                                               
double a = 1.0;


complex<double>ii = complex<double>(0.0,1.0);               // Defining the purely imaginary unit vector.                                                                                               


cx_mat d;
cx_mat dDag;




 complex<double> g;
 complex<double> gDag;

 complex<double> t;
 complex<double> tDag;

 complex<double> L;
 complex<double> LDag;





 

 ofstream outSxx;
 outSxx.open("FM_Server_data_xx/Rep_xx.dat");



 ofstream outSyy;
 outSyy.open("FM_Server_data_yy/Rep_yy.dat");



 ofstream outSzz;
 outSzz.open("FM_Server_data_zz/Rep_zz.dat");



 ofstream outSyzzy;
 outSyzzy.open("FM_Server_data_yzzy/Rep_yzzy.dat");





 // Instead of creating an array with each element being a matrix, we create one large matrix with rows 0 - 3, the first matrix; rows 4 - 7, the second matrix; rows 8 - 11, the third matrix and so on, all the way up to 4M - 1, which will be specfied by the user  (useful to use the sub matrix method in armadillo to compare matrices in C++ and those in Mathematica). The dimension of this matrix will be 2M x 2. Another option is to overwrite the matrices in each subsequent iteration. This might the way to go.    




 //  int J = 300;



 double EnM[J];




// Below, the subatrix method is exploited .submatrix(first row, first column, last row, last column). Using this in conjunction with the large matrices defined above, one can perform the recursive algortihm simimlar to that in Mathemaitca. It will be something like the FIRST row will be in the set {0, 4, 8, 12, 16..} and the LAST row will be in the set {3, 7, 11, 15, ...}; the first column and last column will stay the same, 0 and 1 respectively. However, because the indexing in the recursive relation start at 1, the indices in submat must be as follows: submat(4*i,0, (4*i)+3 ,3).In the recursion relation, we wish to relate one submatrix with the submatrix just before it. If we have submat(12, 0, 15, 3) we want to relate to submat(8, 0, 11 , 3). In general, we are relatingsubmat (4*i, 0 , 4*i + 3, 3) to submat( 4*(i - 1), 0 , 4*(i -1), 3  )  


 cx_mat fineSM(4*J,4); // This matrix is used to store all the final value of epsilonSM for a given energy value. J is included instead of M since we are dealing with given energies.

 cx_mat G00tab(4*J,4); // This matrix is stores the value of the Green function for a given energy value.



 // Q_{||} or Q_z

 // int H = 300 ;

 double QQ[H];





 // Q_{perp} or Q_x                                                                                                                                                                                       

 int F = 1;
 double QP[F];






 // SxxNum3D: Before we had a 1D array, now we have a matrix. Each row represents a given wave vector value, while a given column represents a particular energy.  



 mat Sxx (H,J);

 mat Syy (H,J);

 mat Szz (H,J);

 mat Syzzy (H,J);





 // g = gamma, t = tau, L = lambda, d = delta 



 
 for( int u = 0; u < H; u++ )
{


  

  //   QQ[u] = 1.4;
    QQ[u] = (u*PI)/(H-1);


      

 
  
  g = ((Jy/4.0) - (Jx/4.0))*(1.0 + cos(2.0 * QQ[u] * a ) - ii * sin( 2.0 * QQ[u] * a ) );

  gDag = ((Jy/4.0) - (Jx/4.0))*(1 + cos(2.0 * QQ[u] * a ) + ii * sin( 2.0 * QQ[u] * a ));


  t = ((gfac * mu * B)/2.0) * (1.0 + cos( QQ[u] * a ) - ii * sin( QQ[u] * a ));

  tDag = ((gfac * mu * B)/2.0) * (1.0 + cos( QQ[u] * a ) + ii * sin( QQ[u] * a ));


  L =  ((ii*  DM)/2.0) *( cos(QQ[u]*a)  - ii*sin(QQ[u]*a) - 1.0 );

  LDag =  (( -ii*  DM)/2.0) *( cos(QQ[u]*a) +  ii*sin(QQ[u]*a) - 1.0 ); 








  d << 0 << 0 << 0 << 0 << endr
    << 0 << 0 << 0 << 0 << endr
    << g << 0 << 0 << 0 << endr
    << t + L  << g << 0 << 0 << endr; 



 dDag << 0 << 0 << gDag  << tDag + LDag << endr
      << 0 << 0 << 0 << gDag << endr
      << 0 << 0 << 0 << 0 << endr
      << 0 << 0 << 0 << 0 << endr;  


 









 for(int n = 0; n < J; n++)
   {
 
     
     // EnM[n] = 15.0;
   
      EnM[n] = 9.5 + 9.5*((double)n/ J);   // Initally I had brackets around n and J: (n/J), but that was not working with typecasting. I then removed the parathenses and it works. Note that EnM is in units of meV.



      // k = 9

      fineSM.submat(4*n , 0, (4*n + 3), 3) =  esk( 2, 0, t,  tDag, L,  LDag,  g,   gDag, d, dDag, EnM[n], QQ[u]);
	

      G00tab.submat(4*n, 0, (4*n + 3) , 3) = inv(EnM[n] * id - fineSM.submat(4*n, 0, (4*n + 3) ,3 ) + ii * eta * id  );


       cout << (n+1) + J*u << endl;





      





 

      
 


 // Note that the array QQ is not QQ[j] or even QQ[]. I think having QQ is related to the pointer aspect of the array and the function is using this interpretation.



      
       
 Sxx(u,n) = SxxH(QQ, u, Jy, Jx, Jz, B, G00tab, n);


 /*
 Syy(u,n) = SyyH(QQ, u, Jy, Jx, Jz, B, G00tab, n);

 Szz(u,n) = SzzH(QQ, u, Jy, Jx, Jz, B, G00tab, n);

 Syzzy(u,n) = SyzzyH(QQ, u, Jy, Jx, Jz, B, G00tab, n);
     

 */


 // Exporting data


  outSxx << QQ[u] << " "  << EnM[n] << " "  << Sxx(u,n) << endl;
   
 

   }


 


 }








 time2 = clock();

 double ttime;

 ttime = ((double)(time2 -time1))/CLOCKS_PER_SEC;

 cout << ttime << " sec  "   << endl;



 return 0;


}






// The difference now is that there is a new term in the HSF function, namely double QQ. This should be placed in all other functions. While it will not directly affect the recursion relation, it will allow one to avoid mistakes.

 
cx_mat HSF(int p, complex<double> t, complex<double> tDag, complex<double> L,  complex<double> LDag, complex<double> g, complex<double> gDag, double QQ)
{


  cx_mat HZS;
  cx_mat Hfill;


  if(p==0)
    {


      HZS << Jz + gfac * mu * Bz * (1.0 + 4.0*p) -( (Jx + Jy)/2.0 * cos(QQ) ) << t+L << g << 0 << endr
	 << tDag+ LDag << Jz + gfac * mu * Bz * (2.0 + 4.0*p) << t + L << g << endr
	 << gDag << tDag + LDag << Jz + gfac * mu * Bz * (3.0 + 4.0*p ) << t + L << endr
	 << 0 << gDag << tDag + LDag << Jz + gfac * mu * Bz * (4.0 + 4.0*p) << endr;

      return HZS;
 

    }




else
  {


  Hfill << Jz + gfac * mu * Bz * (1.0 + 4.0*p)  << t+L << g << 0 << endr
	 << tDag+ LDag << Jz + gfac * mu * Bz * (2.0 + 4.0*p) << t + L << g << endr
	 << gDag << tDag + LDag << Jz + gfac * mu * Bz * (3.0 + 4.0*p ) << t + L << endr
	 << 0 << gDag << tDag + LDag << Jz + gfac * mu * Bz * (4.0 + 4.0*p) << endr;


  return Hfill;
  }


  

}







// The following functions encapsulate the intial matrices for the algorithm.



// From p. 992 in notes. p corresponds to n, the input/ independent variable.




// Instead of passing an array into the function (EnM[n]), a double will passed instead


cx_mat es0(int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ) 
{



  return   HSF(p , t, tDag, L,  LDag,   g, gDag, QQ )  + d * inv( EnM * id  - HSF((p+1), t, tDag, L,  LDag,  g,  gDag, QQ )  + ii * eta * id ) * dDag; 

 


}



// Keep in mind that p should always be one less than the integer you wish to input as a0 is defined as a0(n+1) and not a0(n) (refer to 992). I do this for consistency.




cx_mat a0(int p,  complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag, double EnM, double QQ)
{

  return d * inv( EnM * id -  HSF( p + 1, t, tDag, L, LDag, g, gDag, QQ ) + ii * eta * id   ) * d;


}






 
cx_mat e0 (int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ) 
{


 

  return  HSF(p, t, tDag,  L,   LDag,  g, gDag, QQ ) + d * inv( EnM * id  -  HSF((p +1), t, tDag, L,  LDag,  g,  gDag, QQ )  + ii * eta * id ) * dDag +  dDag * inv( EnM * id  -  HSF((p - 1), t, tDag, L,  LDag,  g,  gDag, QQ )  + ii * eta * id ) * d; 




}







cx_mat b0(int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag,  cx_mat d, cx_mat dDag, double EnM, double QQ)
{


  return dDag * inv( EnM * id -  HSF((p - 1), t,   tDag,  L, LDag,  g,  gDag, QQ ) + ii * eta * id   ) * dDag;



}



















// n => p




cx_mat ak(int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ) 
{






  if(k == 0)
    {

      return a0(p, t, tDag, L, LDag, g, gDag,  d,  dDag,  EnM, QQ);

    }



  // p . 998 of notes

  else
    {

      return ak(k-1, p, t, tDag, L, LDag, g, gDag, d, dDag, EnM, QQ)   * inv( EnM * id - ek(k-1, p + pow(2.0,k) ,  t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ)  + ii * eta*id  ) * ak(k - 1, p + pow(2.0,k), t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ);

    }





}











cx_mat  bk (int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ) 
{


  if( k == 0  )
    {

      return b0(p, t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ);

    }



  else
    {


      return bk(k-1, p, t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ) * inv(EnM * id  -  ek(k-1, p - pow(2.0,k) , t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ)  + ii * eta*id ) * bk(k - 1, p - pow(2.0,k), t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ );    


    }






}









cx_mat ek(int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ) 
{


  if( k == 0)
    {


      return e0(p, t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ );    


    }




  else
    {



      return ek(k - 1, p, t, tDag, L, LDag, g, gDag, d, dDag, EnM, QQ ) + ak(k-1, p, t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ) * inv( EnM*id - ek( k - 1, p + pow(2.0,k), t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ ) + ii * eta * id ) * bk( k - 1, p + pow(2.0,k), t, tDag, L, LDag, g, gDag,  d,  dDag, EnM, QQ )  + bk( k - 1, p,  t, tDag, L, LDag, g, gDag,  d,  dDag,  EnM, QQ ) * inv( EnM * id - ek( k - 1, p - pow(2.0,k), t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ )  + ii * eta * id  ) * ak(k - 1, p - pow(2.0,k), t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ );     



    }






}














cx_mat esk(int k, int p, complex<double> t,  complex<double> tDag,  complex<double> L,  complex<double> LDag,  complex<double> g,  complex<double> gDag, cx_mat d, cx_mat dDag,  double EnM, double QQ) 
{


  if( k == 0)
    {


      return es0(p, t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ );    

    }




  else
    {


      return esk(k - 1, p, t, tDag, L, LDag, g, gDag, d,  dDag, EnM, QQ ) + ak(k - 1, p, t, tDag, L, LDag, g, gDag, d,  dDag,  EnM, QQ ) * inv(EnM *  id - ek(k - 1, pow(2.0,k),t, tDag, L, LDag, g, gDag, d, dDag, EnM, QQ) + ii * eta * id  ) * bk( k - 1, pow(2.0,k), t, tDag, L, LDag, g, gDag, d, dDag, EnM, QQ ); 


    }





}











								







   








