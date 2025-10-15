/* -------------------
  Newton interpolation
  through N+1 points
----------------------*/

#include <iostream>
#include <iomanip>

using namespace std;

int main()
{

const int N=3;
double x[N+2], f[N+2], c[N+2], u[N+2];
double A[N+2][N+2];
int i, j, m;

/* -------------------
formating
----------------------*/

cout << setiosflags(ios::fixed | ios::showpoint);

/* -------------------
define the N+1 data points
----------------------*/

x[1]=0.0;
x[2]=0.6;
x[3]=1.2;
x[4]=1.8;

f[1]=0.0;
f[2]=1.0;
f[3]=2.0;
f[4]=5.0;

/* -------------------
display the data
----------------------*/

cout<<" \nInterpolated data:\n\n";
for(i=1; i<=N+1; i++)
{
 cout << setprecision(6) << setw(12) << right << x[i]
      << setprecision(6) << setw(12) << right << f[i] <<"\n";
}

/* -------------------
Newton algorithm
----------------------*/

for(i=1; i<=N+1; i++)
{
 u[i] = f[i];
 A[i][1] = u[i];  // Optional, used for the Newton table
}

c[1] = u[1];

for(m=2; m<=N+1; m++)  // run over columns
{
  for(i=1; i<=N-m+2; i++)  // run over rows
  {
  u[i] = (u[i+1]-u[i])/(x[i+m-1]-x[i]);
  A[i][m] = u[i]; // Optional, used for the Newton table
  }
  c[m] = u[1];
}

/* -------------------
print the coefficients
----------------------*/

cout<<"\nNewton coefficients:\n\n";

for(i=1; i<=N+1; i++)
{
  cout << setprecision(6) << setw(12) << right << c[i] << "\n";
}

/* -------------------
print the table
----------------------*/

cout<<"\nNewton table:\n\n";

for(i=1; i<=N+1; i++)
{
  for(j=1; j<=N+2-i; j++)
  {
   cout<< setprecision(6) << setw(12) << right << A[i][j];
  }
  cout << "\n";
}

return 0;
}
