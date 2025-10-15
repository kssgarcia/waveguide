#include <stdio.h>
#include <math.h>
#include <stdlib.h>
 
void quadratic (double r, double s);
 
int main (void)
{
 
   double *a, *b, *c, r=1.0, s=1.0, dr=1.0, ds=0.0, det, detr, dets;
   int deg,i,cycle=1;
 
   printf("Please enter the degree of the polynomial\n");
   scanf("%d", &deg);
 
   a=(double*)calloc(deg+1,sizeof(double));
   b=(double*)calloc(deg+1,sizeof(double));
   c=(double*)calloc(deg+1,sizeof(double));
                                                                                
   printf("Please enter the coeffecients of the polynomial starting from the
highest order:\n");

        for (i=0;i<deg+1;i++)
                scanf("%lf", &a[i]);
 
        while(deg>2)
        {
 
                b[0]=a[0];
                c[0]=0.0;
                c[1]=b[0];
 
                        while ((fabs(dr)+fabs(ds))>1e-16)
                        {
 
                                b[1]=a[1]+r*b[0];
                                c[2]=b[1]+r*c[1];
                                for (i=2;i<deg+1;i++)
                                        b[i]=a[i]+r*b[i-1]+s*b[i-2];
                                for (i=3;i<deg+1;i++)
                                        c[i]=b[i-1]+r*c[i-1]+s*c[i-2];
 
                                det=c[deg-1]*c[deg-1]-c[deg]*c[deg-2];
                                detr=c[deg-2]*b[deg]-b[deg-1]*c[deg-1];
                                dets=c[deg]*b[deg-1]-b[deg]*c[deg-1];
 
                                if (fabs(det) < 1e-16)
                                {
                                        det = 1.0;
                                        detr = 1.0;
                                        dets = 1.0;
                                }
 
                                dr=detr/det;
                                ds=dets/det;
                                r=r+1.0*dr;
                                s=s+1.0*ds;
                                cycle++;
                        }
 
                for (i=0;i<deg-1;i++)
                        a[i]=b[i];
                quadratic(r,s);
 
                deg-=2;
                cycle=1;
        }
        //add last case
 
        return 0;
}
 

void quadratic (double r, double s)

{
  double s1,s2;
 
  if((r*r-4*(-s))>0.0)
   {
     s1=(r+sqrt(r*r-4.0*(-s)))/2.0;
     s2=(r-sqrt(r*r-4.0*(-s)))/2.0;
     printf("The roots are %.3f and %.3f\n", s1,s2);
   } /* if */
  else
   {
     s1=r/2.0;
     s2=sqrt(-1.0*(r*r-4.0*(-s)))/2.0;
     printf("The roots are %.3f + %.3fi and %.3f - %.3fi\n", s1,s2,s1,s2);

   } /* else */
 
  return;
}
