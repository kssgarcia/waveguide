/*---------------------------------------------
  Compute the binomial coefficient in two ways
 --------------------------------------------*/
 
#include<iostream>
using namespace std;

int n, m;

void GetNM();
 
/*---------- main --------------*/

int main()
{
     int i, nf, mf,nmmf, l, k, nmc;
 
     while(1)
     {
       cout<<"\nWill compute the (n/m) binomial coefficient\n";
       GetNM();
 
       if ( (n==0) || (m==0) )
       break;
 
        // direct method using factorials
 
          nf = 1;      //compute n!
          for (i = 1; i<=n; i++)
               nf = nf*i;

          mf = 1;      //compute m!
          for (i = 1; i<=m; i++)
               mf = mf*i;
 
          nmmf = 1;    //compute (n-m)!
          for (i = 1; i<=(n-m); i++)
               nmmf = nmmf*i;
 
          nmc = nf/(mf*nmmf);

          cout<<"\nBinomial coefficient by the direct method:\t";
          cout<<nmc<<"\n";
 
   // better method
   // find the minimum of m and n-m
 
        l = m;
        if((n-m) < l)
        l = n-m;
 
        //some magic
 
        nmc = 1;
        for(k = 1; k<=l;k++)
          nmc = nmc*(n-k+1)/k;
 
        cout<<"\nBinomial coefficient by the improved method:\t";
        cout<<nmc<<"\n";
 
    }     //pack it in
 
 return 0;
}

/*---------- GetNM --------------*/

void GetNM()
{
  cout<<"\nPlease enter n and m (n>=m)";
  cout<<"\n   0 for either one to quit";
  cout<<"\n---------------------------\n";
  cin>>n; cin>>m;
 
  if( (m > n) && (m!=0) && (n!=0) )
    {
     cout<<"\nInvalid input; please try again\n";
     GetNM();
    }
}
