/*------------------------------
  roots of a quadratic equation
 ------------------------------*/

#include<iostream>
#include<cmath>

using namespace std;

int main()
{
   double a,b,c,d,a2,srd,x1,x2,imag,real,root;
   cout<<"\nRoots of: a*x**2 + b*x + c = 0\n";

   while(1)
   {
      cout<<"\n Please enter the coefficients a, b, c";
      cout<<"\n 999 for any one to quit";
      cout<<"\n------------------------\n";
      cin>>a;
      cin>>b;
      cin>>c;
      if( (a==999) || (b==999) || (c==999) )
      break;

      /*------
      c traps
      c------*/

      if(abs(a) < 0.0000001)
        {
         cout<<"\nThe coefficient a is too small";
         cout<<"\nand the equation looks linear\n";

           if(abs(b) > 0.0000001)
            {
             root = -c/b;
             cout<<"\nThe root is: "<<root<<"\n";
            }
        else
          cout<<"\nNo root for the constant binomial\n";
        }

/*-------------------
  analytical solution
 --------------------*/

     else
     {
       d = b*b-4.0*a*c;   //discriminant
       a2 = 2.0*a;
       if(d>=0)
        {
         srd = sqrt(d);
         x1=(-b+srd)/a2;
         x2=(-b-srd)/a2;
         cout<<"\nroots are:"<<x1<<", "<<x2<<"\n";
        }
         else
        {
          d = -d;
          imag = sqrt(d)/a2;
          real = -b/a2;
          cout<<"\n roots are complex conjugate";
          cout<<"\n real part:\t"<<real;
          cout<<"\n imaginary part:\t"<<imag<<"\n";
        }
       }

    } //end while-loop

 return 0;
 }
