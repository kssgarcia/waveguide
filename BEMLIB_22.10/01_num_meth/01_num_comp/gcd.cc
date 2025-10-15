/*------------------------------------
  This program finds the
  greatest common divisor
  using Euclid's algorithm
------------------------------------*/

/*-------------------------------------------------
  This program accompanies the book:
             C. Pozrikidis
  "Numerical Computation in Science and Engineering"
            Oxford University Press
---------------------------------------------------*/

#include <iostream>

using namespace std;

int main()
{
  int n,m,k,nave;

  cout<<"\nWill compute the Greatest Common Divisor";
  cout<<"\n\tof two positive integers\n";

  cout<<"\nPlease enter the two integers";
  cout<<"\n\t0 for either one to quit\n";

  cin>>n;
  cin>>m;

  while ( (n!=0) && (m!=0) )
  {
    if(n==m)
    { k=n;
      cout<<"\nThe Greatest Common Divisor is: "<<k<<"\n";
    }
    else while (n!=m)
    {  if(n>m)
       { nave=m;
         m=n;
         n=nave;
       }
       k=m-n;
       m=n;
       n=k;
       if(n==m)
       { k=n;
         cout<<"\nThe Greatest Common Divisor is: "<<k<<"\n";
       }
    }

    cout<<"\nPlease enter the two integers";
    cout<<"\n\t0 for either one to quit\n";

    cin>>n;
    cin>>m;

  }      //end of while-loop

cout<<"\n Thank you for your business.\n";

return 0;
}
