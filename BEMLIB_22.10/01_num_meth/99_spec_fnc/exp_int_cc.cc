
/*------------------------------------------------
 This program accompanies the book:
           C. Pozrikidis
 Numerical Computation in Science and Engineering
          Oxford University Press
------------------------------------------------*/

/*----------------------------------
 exp_int_cc(x)

 Polynomial and rational expression
 for the exponential integral

 Abromowitz & Stegun, p. 231
----------------------------------*/

#include<math.h>
using namespace std;

double exp_int_cc( double x)
{
  double a0, a1, a2, a3, a4, a5, x2, x3, x4, x5;
  double expint, b1, b2, b3, b4;

      if(x <= 1.0)
      {
       a0 = -0.57721566;
       a1 =  0.99999193;
       a2 = -0.24991055;
       a3 =  0.05519968;
       a4 = -0.00976004;
       a5 =  0.00107857;

       x2 = x*x;
       x3 = x*x2;
       x4 = x*x3;
       x5 = x*x4;

       expint = -log(x)+a0 + a1*x + a2*x2 +
                 a3*x3 + a4*x4 + a5*x5;

/*c--------
c lower-order accuracy
c
c        expint = -log(x);
c        sign = 1.0;
c        fact = 1.0;
c        power = 1.0;
c        for (i = 1; i<=10; i++)
c                                  {
c          fact = fact*i;
c          power = power*x;
c          expint = expint + sign*power/(i*fact);
c          sign = -sign;
c        }
c--------*/
          }
          else
          {
/*c--------
c lower-order accuracy
c
c        a1 = 2.334733;
c        a2 = 0.250621;
c        b1 = 3.330657;
c        b2 = 1.681534;
c        x2 = x*x;
c        expint = (x2+a1*x+a2)/(x2+b1*x+b2)/(x*exp(x));
c--------*/

          a1 = 8.5733287401;
          a2 =18.0590169730;
          a3 = 8.6347608925;
          a4 =  .2677737343;

          b1 = 9.5733223454;
          b2 =25.6329561486;
          b3 =21.0996530827;
          b4 = 3.9584969228;

          x2 = x*x;
          x3 = x*x2;
          x4 = x*x3;
 
          expint =(x4+a1*x3+a2*x2+a3*x+a4)
           /(x4+b1*x3+b2*x2+b3*x+b4)/(x*exp(x));
        }

      return expint;
}
