/*--------- function powi -----------------
  computes the integral power i**j or i^j
--------------------------------------*/

int powi(int i, int j)
{
  int k;
  int accum = 1;
  for(k=1; k<=j; k++)
   {
    accum = accum * i;
   }
return accum;
}
