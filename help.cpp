#include <iostream>
using namespace std;

void growth_help()
{
 cout << "\n USAGE                              src - C \n"
      << " ---------------------------------------------\n";
 cout << " ostwald -Options\n"
      << " ---------------------------------------------\n\n"
      << " CALCULATION MODE              -m NGG (or AGG)\n"
      << " SCREW-DISLOCATION EFFECT      -w\n"
      << " TEMPERATURE      (1100C - 1700C)      -t 1500\n"
      << " STEP FREE ENERGY (0.0hG - 1.0hG)      -s 0.33\n"
      << " LIQUID FRACTION  (0.00L - 1.00L)      -l 0.46\n"
      << " ---------------------------------------------\n"
      << " CALCULATION TIME STEP (slow, 0.01 s)  -z\n"
      << " CALCULATION TIME STEP (fast, 1 s)     -q\n"      
      << " ---------------------------------------------\n"
      << " GENERATION OF INITIAL GRAIN SET       -i\n"
      << " RE-CALCULATION OF HISTOGRAM BINS      -k\n"
      << " ---------------------------------------------\n"
      << "              Edward Yang-il Jung @ 2008-12-28\n"
      << endl;
}
