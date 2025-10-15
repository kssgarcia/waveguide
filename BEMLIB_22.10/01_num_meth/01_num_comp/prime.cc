/*---------------------------------------
* Check whether a given integer is prime
---------------------------------------*/

#include<iostream>

using namespace std;

int n;

/*------input-------------*/

int GetN()
{
   n=-1;
   while (n<0)
   {
     cout<<"\nPlease enter the integer to be tested";
     cout<<"\n0 to quit";
     cout<<"\n----------\n";
     cin>>n;
     if(n<0) cout<<"\nThe integer must be positive; try again\n";
   }
 return n;
}

/*------main-------------*/

int main()
{
  int k,l,m;
  n=1;

  while (n!=0)
  {
     n=GetN();
     for(m=2; m<=n-1; m++)
     {
       l=n/m; //testing for the remainder
       k=l*m;
       if(k==n)
       {
        cout<<"\n"<<n<<" is not a prime number";
        cout<<"\nIt's highest divisor is "<<l<<"\n";
        break;
       }
     }
     if( (k!=n) && (n!=0) || (n==2))
     cout<<"\n"<<n<<" is a prime number\n";
  }
return 0;
}
