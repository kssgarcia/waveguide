/* gauss_sph_cc_dr*/

#include <iostream>
#include "gauss_sph_cc.h"

using namespace std;

int main()
{
  long int nn;
  double xxi[50], yyi[50], zzi[50], wwi[50];

  nn = 4;

  gauss_sph_cc (nn, xxi, yyi, zzi, wwi);
}
