/*------------------------------------------------
 This program accompanies the book:
      C. Pozrikidis
  Numerical Computation in Science and Engineering
     Oxford University Press
------------------------------------------------*/
 
/*-------------------------------------------
   ell_int_cc
   Evaluation of complete elliptic integrals
   of the first and second kind
   using a recursion formula
 
   SYMBOLS:
   --------
 
   F:     first kind
   E:     second kind
   acc: specified accuracy
 -------------------------------------------*/
 
#include<cmath>
using namespace std;
 
void ell_int_cc (double RK2, double F, double E)
{
   double acc, pi, pih;
   double B, C, D, G, RK;
   acc=0.000000001;
 
     // constants
        pi  = 3.14159265358;
        pih = 0.5*pi;
 
     // launching
        RK = sqrt(RK2);
        F  = pih;
        E  = 1.0;
        G  = 1.0;
        B  = RK;
 
        D = 1; //initialize D for the while-loop
 
        while (abs(D) > acc)
          {
            C = sqrt(1.0-(B*B));
            B = (1.0-C)/(1.0+C);
            D = F*B;
            F = F+D;
            G = 0.5*G*B;
            E = E+G;
           }
 
        E=F*(1.0-0.5*RK2*E);
 
} //Done
