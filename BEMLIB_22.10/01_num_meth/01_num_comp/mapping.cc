/* ---------------------------------------------
Multiply a vector by a square matrix many times

  This program accompanies the book:
            C. Pozrikidis
  ``Numerical Computation in Science and Engineering''
          Oxford University Press
------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

/*-------------------------------------------------
                mat_vec
 function to perform matrix-vector multiplication:
	c_i = a_ij b_j
--------------------------------------------------*/

void mat_vec
	(int n
	,double a[][50]
	,double b[]
	,double c[]
	)
{
int i; int j;
  for (i=1;i<=n;i++)
    {c[i] = 0;
     for (j=1;j<=n;j++)
     {c[i] = c[i] + a[i][j] * b[j];
     }
    }
}

/*-------------------------------------------------
                  main
--------------------------------------------------*/

int main()
{
int norm; int n; int i; int j;
double a[50][50];
double b[50];
double c[50];

cout << "\n";
cout << " Normalize the vector after each projection?\n";
cout << " Enter 1 for yes, 0 for no\n";
cout << " -------------------------\n";
cin >> norm;

/*--- READ THE MATRIX AND THE VECTOR */

ifstream input_data;
input_data.open("matrix_v.dat");

input_data >> n;

for (i=1;i<=n;i++)
	{for (j=1;j<=n;j++)
		{input_data >> a[i][j];
		}
	}

for (i=1;i<=n;i++)
	{input_data >> b[i];
	}

input_data.close();

/*--- DONE READING */

/*--- DISPLAY*/

cout << "\n";
cout << " Matrix - initial vector:";
cout << "\n\n";

for (i=1;i<=n;i++)
	{for (j=1;j<=n;j++)
		{cout << setw(8) << a[i][j];
		}
	cout << " " << setw(8) << b[i] << "\n";
	}
cout << "\n\n";

/*--- END DISPLAYING*/

/*--- MAPPING*/

int icount = 0;
int more = 1;

while (more!=0)
  {
	mat_vec (n,a,b,c);

	for (i=1;i<=n;i++)
		{b[i]=c[i];
		}	

	if(norm == 1)
		{double rnorm = 0;
		for (i=1;i<=n;i++)
			{rnorm = rnorm + b[i]*b[i];
			}
		rnorm = sqrt(rnorm);
		for (i=1;i<=n;i++)
			{b[i]=b[i]/rnorm;
			}
		}

	cout << " Projected vector at stage: " << icount;
	cout << "\n\n";

	for (i=1;i<=n;i++)
		{
		cout << setprecision(5) << setw(10);
		cout << b[i] << "\n";
		}

	icount = icount+1;

	cout << " One more projection?\n" ;
	cin >> more;
  }

return 0;
}
