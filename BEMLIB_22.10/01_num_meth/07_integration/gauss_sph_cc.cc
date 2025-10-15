/*c----------------------------------------
c Abscissae and weights for integration
c over the surface of a sphere
c----------------------------------------*/

/*c-----------------------------------------
c Copyright by C. Pozrikidis, 1999
c All rights reserved.
c
c This program is to be used only under the
c stipulations of the licensing agreement.
c----------------------------------------*/

/*c------------------------------------------------
c This program accompanies the book:
c
c C. Pozrikidis
c
c Numerical Computation in Science and Engineering
c
c Oxford University Press
c
c Fortran to cc translation by:
c Rhodalynn Degracia
c------------------------------------------------*/

#include <iostream>
#include <math.h>

using namespace std;

void gauss_sph_cc(long int Nint, double xx[], double yy[], double zz[], double ww[])
{
		double srti, two;

		if( (Nint!=6)  && (Nint!=18) )
		{
		 cout<<"\n Requested number of points";
		 cout<<"\n not available ";
	         cout<<"\n Will take Ninit = 18"<<"\n";
		 Nint = 18;
		}

		/*-----------------------*/
		if(Nint==6)
		/*-----------------------*/
		{
			xx[0] = 1.0;
			yy[0] = 0.0;
			zz[0] = 0.0;
			ww[0] = 1/6.0;

			xx[1] = -1.0;
			yy[1] = 0.0;
			zz[1] = 0.0;
			ww[1] = ww[0];

			xx[2] = 0.0;
			yy[2] = 1.0;
			zz[2] = 0.0;
			ww[2] = ww[0];

			xx[3] = 0.0;
			yy[3] = -1.0;
			zz[3] = 0.0;
			ww[3] = ww[0];

			xx[4] = 0.0;
			yy[4] = 0.0;
			zz[4] = 1.0;
			ww[4] = ww[0];

			xx[5] = 0.0;
			yy[5] = 0.0;
			zz[5] = -1.0;
			ww[5] = ww[0];
		}

		/*-----------------------*/
		else if (Nint==18)
		/*-----------------------*/
		{
			two = 2.0;
			srti = 1.0/sqrt(two);

			xx[0] = 1.0;
			yy[0] = 0.0;
			zz[0] = 0.0;
			ww[0] = 1.0/30.0;

			xx[1] = -1.0;
			yy[1] = 0.0;
			zz[1] = 0.0;
			ww[1] = ww[0];

			xx[2] = 0.0;
			yy[2] = 1.0;
			zz[2] = 0.0;
			ww[2] = ww[0];

			xx[3] = 0.0;
			yy[3] = -1.0;
			zz[3] = 0.0;
			ww[3] = ww[0];

			xx[4] = 0.0;
			yy[4] = 0.0;
			zz[4] = 1.0;
			ww[4] = ww[0];

			xx[5] = 0.0;
			yy[5] = 0.0;
			zz[5] = -1.0;
			ww[5] = ww[0];

			xx[6] = srti;
			yy[6] = srti;
			zz[6] = 0.0;
			ww[6] = 1.0/15.0;

			xx[7] = -srti;
			yy[7] = srti;
			zz[7] = 0.0;
			ww[7] = ww[6];

			xx[8] = srti;
			yy[8] = -srti;
			zz[8] = 0.0;
			ww[8] = ww[6];

			xx[9] = -srti;
			yy[9] = -srti;
			zz[9] = 0.0;
			ww[9] = ww[6];

			xx[10] = 0.0;
			yy[10] = srti;
			zz[10] = srti;
			ww[10] = ww[6];

			xx[11] = 0.0;
			yy[11] = -srti;
			zz[11] = srti;
			ww[11] = ww[6];

			xx[12] = 0.0;
			yy[12] = srti;
			zz[12] = -srti;
			ww[12] = ww[6];

			xx[13] = 0.0;
			yy[13] = -srti;
			zz[13] = -srti;
			ww[13] = ww[6];

			xx[14] = srti;
			yy[14] = 0.0;
			zz[14] = srti;
			ww[14] = ww[6];

			xx[15] = -srti;
			yy[15] = 0.0;
			zz[15] = srti;
			ww[15] = ww[6];

			xx[16] = srti;
			yy[16] = 0.0;
			zz[16] = -srti;
			ww[16] = ww[6];

			xx[17] = -srti;
			yy[17] = 0.0;
			zz[17] = -srti;
			ww[17] = ww[6];
		}
}
