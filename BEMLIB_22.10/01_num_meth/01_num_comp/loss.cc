/* gradual loss of significance */

#include <iostream>
#include <iomanip>

using namespace std;

int main()
{
double x = 1.0/9.0;
for (int i=0; i<30; i++)
  {cout << setw(20) <<  setprecision(18) << x << endl;
   x = (x-1)*10; // C treats as (x-1.0)*10.0
  }
return 0;
}
