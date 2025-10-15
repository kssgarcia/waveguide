#include<iostream>
#include<iomanip>

using namespace std;

/* sum a series  that decays like 1/i^2
using Aitken extrapolation */

//------- main-------------

int main()
{

const int kmax = 15;
double x[kmax+1], a[kmax+1];
double sum = 0.0;
int N=1, p=2, M=N, L=0, k, i;

for(k=0; k<=kmax; k++)
 {
  for(i=L+1;i<=M;i++)
   {
    sum = sum+1.0/(i*i);
   }
   x[k] = sum-0.5/(M*M);
   if(k>1)
   {
    a[k-1] = (x[k-2]*x[k]-x[k-1]*x[k-1])/(x[k]-2.0*x[k-1]+x[k-2]);
   }
   L=M;
   M=M*p;
 }

//--- print the sequences:

  double pi=3.14159265358979323846;
  double exact=pi*pi/6.0;

  a[0]=a[kmax]=0.0;
  cout << setprecision(8) << setw(12) << setiosflags(ios::fixed);
  cout << setw(3) << 0 << "  " << x[0]  << "  " << x[0]-exact << endl;

  for(k=1;k<=kmax-1;k++)
  {
   cout << setw(3) << k << "  " << x[k] << "  " << x[k]-exact 
        << "  " << a[k] << "  " << a[k]-exact << endl;
  }
  cout << setw(3) << kmax << "  " << x[kmax] << "  " << x[kmax]-exact << endl;

return 0;
}
