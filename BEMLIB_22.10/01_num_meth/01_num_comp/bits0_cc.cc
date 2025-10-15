/* ---------------- ----------------------------
  Greatest integer that can be described by n bits
------------------------------------------------*/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

/*---------------- main program -------- */

int main()
{
 int n=1;
 int i;
 const int two = 2;
 cout << " Will compute the greatest integer\n";
 cout << " that can be described with n bits\n";

while(n!=0)
  {
     cout << "\n";
     cout << " Please enter the number of bits\n";
     cout << "        (should be less than 32)\n";
     cout << " q to quit\n";
     cout << " -----------\n";

     if(!(cin >> n)) break;

     cout << setw(13) << "   bits  " << "   "
          << setw(10) << " increment" << "  "
          << setw(16) << "largest integer" << "\n";

     int q = 0;

     for(i=0; i<=n-1; i++)
     {
      int p = pow(two, i);
      q = q + p;
      cout << setw(10) << i+1 << "  "
	   << setw(10) << p   << "  "
	   << setw(10) << q   << "\n";
     };
  };

return 0;
}
