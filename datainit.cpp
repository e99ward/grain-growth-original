// setting the grain size data
#include <iostream>
#include <fstream>
using namespace std;
#include <cmath>
#include <ctime>

const double PI = 4.0* atan(1.0);

double Normal_Distribution(int av, int st);
void Sort_Quick(double num[], int left, int right);

//-------------------------------------

int data_init()
{
 cout << "|                                           |" << endl;
 cout << "|    -> Grain Size Generatior               |" << endl;
 cout << "+-------------------------------------------+" << endl;

 int MAX;
 cout << "Number of Grains:_______\b\b\b\b\b\b\b";
 cin >> MAX;
 double * sizes;
 sizes = new double[MAX];
 char file[15] = "d_0000000.txt";
 srand((unsigned)time(NULL));
 int av, st;
 cout << "Standard Gaussian Distribution" << endl;
 cout << "Average (100=1um):___\b\b\b";
 cin >> av;
 cout << "StandardDeviation:___\b\b\b";
 cin >> st;

 // make normal distribution
 for (int i=0; i<MAX; i++)
  sizes[i] = Normal_Distribution(av,st);
 cout << "--------------------------------------------" << endl;
 cout << "sizes are generated.    now sorting..." << endl;
 // sort the sizes
 Sort_Quick(sizes, 0, MAX-1);
 cout << "sorting finished.    now saving file..." << endl;
 // write to file
 ofstream fout(file,ios_base::out);
 fout.precision(10);
 fout.setf(ios_base::fixed,ios_base::floatfield);
 int k=0;
 for (int j=0; j<MAX; j++)
 {
  fout << sizes[j] << " ";
  k++;
  if ((k%10)==0)
   fout << '\n';
 }
 fout.close();
 cout << "process end." << endl;
 cout << "--------------------------------------------" << endl;
 return 0;
}

double Normal_Distribution(int av, int st)
{
 double r1, r2, z1, z2;
 r1 = (1.0+rand())/(RAND_MAX+2);
 r2 = 1.0*rand()/RAND_MAX;
// r1 = (double)(rand()%RAND_MAX)/RAND_MAX;
// r2 = (double)(rand()%RAND_MAX)/RAND_MAX;
 z1 = sqrt(-2.0 * log(r1)) *  cos(2.0 * PI * r2);
 z2 = z1 * st + av;
 if(z2<0) z2 *= -1.0;
 return z2;
}

void Sort_Quick(double num[], int left, int right)
{
 int l_hold, r_hold;
 double pivot;
 l_hold = left;
 r_hold = right;
 pivot = num[left];
 while (left < right)
 {
  while ((num[right] >= pivot) && (left < right))
   right--;
  if (left != right)
  {
   num[left] = num[right];
   left++;
  }
  while ((num[left] <= pivot) && (left < right))
   left++;
  if (left != right)
  {
   num[right] = num[left];
   right--;
  }
 }
 num[left] = pivot;
 int n_p = left;
 left = l_hold;
 right = r_hold;
 if (left < n_p)
  Sort_Quick(num, left, n_p-1);
 if (right > n_p)
  Sort_Quick(num, n_p+1, right);
}
