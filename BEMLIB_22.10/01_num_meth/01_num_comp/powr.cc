/*---------------powr-----------------
  compute the integral power a**n or a^n
--------------------------------------*/

int powi(double a, int n)
{
  int k;
  double accum = 1.0;
  for(k=1; k<=n; k++)
   {
    accum = accum * a;
   }
return accum;
}
