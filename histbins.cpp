// setting the grain size data
#include <iostream>
#include <fstream>
using namespace std;
#include "convert.h"

int histo_bin()
{
 int cts, range;
 char file[15] = "h_0000000.txt";

 cout << "|                                           |" << endl;
 cout << "|    -> Histogram Re-Counter                |" << endl;
 cout << "+-------------------------------------------+" << endl;
 cout << "Calculation Time Step to process:_______\b\b\b\b\b\b\b";
 cin >> cts;
 cout << "Grain size range to merge(1=0.01um):___\b\b\b";
 cin >> range;
 if (range < 1)
 {
  cout << "range was set to 1" << endl;
  range = 1;
 }
 
 string name;
 name = to_string<int>(cts, dec);
 int h = name.length() - 1;
 int k = 8;
 while(h>=0)
  file[k--] = name[h--];        // step into file name

 ifstream fin;
 fin.open(file,ios_base::in);
 int bins;
 fin >> bins;
 int * count;
 count = new int[bins];
 k = 0;
 while (fin >> count[k] && k<bins)
  k++;
 fin.close();

 file[0] = 'e';
 ofstream hout(file,ios_base::out);
 int s = 0;
 double sum = 0;
 while (s<bins)
 {       
  for (k=0; k<range; k++)
  {
   if (s>=bins)
    break;
   sum += count[s++];
  }
  hout << sum << endl;
  sum = 0.0;
 }
 hout.close();
 
 delete [] count;
 return 0;
} 
