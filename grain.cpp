// specification GRAIN class
#include <iostream>
#include <fstream>
using namespace std;
#include <string>
#include <cmath>
#include "grain.h"
#include "convert.h"
#include "constant.h"
using namespace constant;


Grain::Grain()
{
 MAX = 1;
 g_size = new double[MAX];
 g_size[0] = 0.0;
 h_time = 0.1;
 rs = 0.0;
 rav = 0.0;
 rst = 0.0;
 s_l = 0.0;
 s_e = 0.0;
 s_s = 0.0;
}

Grain::Grain(int step)
{
 char file[15] = "d_0000000.txt";
 string name;
 name = to_string<int>(step, dec);
 int h = name.length() - 1;
 int k = 8;
 while(h>=0)
  file[k--] = name[h--];        // step into file name

 MAX = 1000000;
 g_size = new double[MAX];
 ifstream fin;
 fin.open(file,ios_base::in);
 k = 0;
 while (fin >> g_size[k] && k<MAX)
  k++;
 fin.close();
 MAX = k;
 Grain::avst();
 rs = 0.0;
}

Grain::~Grain()
{
 delete [] g_size;
}

//---------------------------------------------------------

void Grain::spec(bool mode1, bool mode2, int temp, double sfe, double liq, double speed)
{
 s_e = 1.0;
 s_s = 1.0;
 f_A = 1.0;
 f_B = 1.0;
 f_C = 1.0;
 if (mode1==false)                  // grain growth kinetics
  modes = 'n';                      // n = normal
 else if (mode2==false)
  modes = 'a';                      // a = abnormal 
 else
  modes = 's';                      // s = screw-disl. assisted
 if (temp != def_temp)                  // change in temperature
 {
  double t_a, t_A, new_sfe;
  t_A = 1.0 * exp(-1.0*def_act/(def_temp+273));
  t_a = 1.0 * exp(-1.0*def_act/(temp+273));
  f_B = t_a / t_A;
  f_C = 1.0 * (def_temp+273) /(temp+273);
  f_A = f_B * f_C;
  new_sfe = 1.0 * exp(-1.0*def_rta/sqrt(1.0*def_rtt-temp-273));
  s_e = new_sfe * new_sfe / def_sfe / def_sfe; 
  s_s = def_sfe / new_sfe;
 }
 if (sfe != def_sfe)                   // change in step free energy
 {
  s_e = sfe * sfe / def_sfe / def_sfe;          // system step free energy: squared
  s_s = def_sfe / sfe;                 // system step free energy: inverse
 }
 if (liq == 1.0)                    // system liquid fraction: Ardell's
  s_l = 0.0;
 else
  s_l = 26.0 - 25.0*liq;
 h_time = speed;
}

void Grain::avst()
{
 t_grain = 0;
 double t_radius=0, t_log=0;
 for (int i=0; i<MAX; i++)
 {
  if (g_size[i] == 0) continue;
  t_radius += g_size[i];
  t_grain++;
 }
 rav = t_radius / t_grain;
 for (int j=0; j<MAX; j++)
 {
  if (g_size[j] == 0) continue;
  t_log += log10(g_size[j]/rav) * log10(g_size[j]/rav);
 }
 rst = sqrt(t_log / (t_grain-1));
}

void Grain::rstar()
{
 double v_r, v_mass;
 v_r = rav;
 v_mass = Grain::mass(v_r);
 cout << "TRACE: " << v_mass << '\t';  // check for mass
 double d_r = 5.0;
 while(v_mass > 0.1 || v_mass < -0.1)
 {
  if (v_mass > 0)
  {
   v_r += d_r;
   v_mass = Grain::mass(v_r);
   if (v_mass < 0)
    d_r /= 2;
  }
  else
  {
   v_r -= d_r;
   v_mass = Grain::mass(v_r);
   if (v_mass > 0)
    d_r /= 2;
  }
 }
 rs = v_r;
 cout << "rs=" << rs << endl;          // check for r_Star
}

void Grain::calc()
{
 double * g_vary;
 g_vary = new double[MAX];
 Grain::rstar();
 Grain::vary(g_vary, rs);
 for (int j=0; j<MAX; j++)
 {
  g_size[j] += g_vary[j];
  if (g_size[j] < 5.0e-3)       // r<=0.05nm means no grain
   g_size[j] = 0;
 }
 delete [] g_vary;
 Grain::avst();
}

void Grain::writeData(int step)
{
 char file[15] = "d_0000000.txt";
 string name;
 name = to_string<int>(step, dec);
 int h = name.length() - 1;
 int k = 8;
 while(h>=0)
  file[k--] = name[h--];        // step into file name

 int t_zero, redus;
 t_zero = MAX - t_grain;
 if (t_zero > 100)
  redus = int(t_zero / 100) * 100;
 else
  redus = 0;

 ofstream fout(file,ios_base::out);
 fout.precision(10);
 fout.setf(ios_base::fixed,ios_base::floatfield);
 k = 0;
 for (int i=redus; i<MAX; i++)
 {
  fout << g_size[i] << " ";
  if (redus==0)
   k++;
  else
   g_size[k++] = g_size[i];
  if ((k%10)==0)
   fout << '\n';
 }
 fout.close();

 MAX -= redus;
}

void Grain::writeHist(int step)
{
 char file[15] = "h_0000000.txt";
 string name;
 name = to_string<int>(step, dec);
 int h = name.length() - 1;
 int k = 8;
 while(h>=0)
  file[k--] = name[h--];        // step into file name

 int bins = int(g_size[MAX-1])+1;
 int * m_dist;
 m_dist = new int[bins];
 for (int j=0; j<bins; j++)
  m_dist[j] = 0;
 int rels;
 for (int t=0; t<MAX; t++)
 {
  if (g_size[t]>0)
  {
   rels = int(g_size[t]);
   m_dist[rels]++;
  }
 }

 ofstream hout(file,ios_base::out);
 hout << bins;
 for (int s=0; s<bins; s++)
 {
  if ((s%15)==0)
   hout << '\n';
  hout << m_dist[s] << " ";
 }
 hout.close();

 delete [] m_dist;
}

void Grain::writeStat(int step)
{
 static int sett=0;
 ofstream fout("rs_N_avst.txt",ios_base::out | ios_base::app);
 if (sett==0)
 {
  fout << "step\tr_star\tr_ave\tr_std\tr_max\tr_num" << '\n';
  sett++;
 }
 fout << step << '\t' << rs << '\t' << rav << '\t' << rst
      << '\t' << g_size[MAX-1] << '\t' << t_grain <<  '\n';
 fout.close();
}

void Grain::sort(double num[], int left, int right)
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
  sort(num, left, n_p-1);
 if (right > n_p)
  sort(num, n_p+1, right);
}

void Grain::vary(double g_vary[], double r)
{
 double * k1;
 k1 = new double[MAX];
 double * k2;
 k2 = new double[MAX];
 double * k3;
 k3 = new double[MAX];
 double * k4;
 k4 = new double[MAX];
 if (modes == 'n')                // normal grain growth (NGG)
 {
  for (int i=0; i<MAX; i++)
  {
   if (g_size[i] == 0)
    g_vary[i] = 0.0;
   else if (g_size[i] < r)        // NGG dissolution
   {
    g_vary[i] = h_time * c_A * f_A * 1.0/g_size[i] * (1.0/r - 1.0/g_size[i]) * (1.0+ s_l * g_size[i]/r);
    if ((g_vary[i] + g_size[i]) < 0)
     g_vary[i] = -1.0 * g_size[i];
   }
   else                           // NGG growth
   {
    k1[i] = h_time * c_A * f_A * 1.0/g_size[i] * (1.0/r - 1.0/g_size[i]) * (1.0+ s_l * g_size[i]/r);
    k2[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 0.5*k1[i]) * (1.0/r - 1.0/(g_size[i]+ 0.5*k1[i])) * (1.0+ s_l * (g_size[i]+ 0.5*k1[i])/r);
    k3[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 0.5*k2[i]) * (1.0/r - 1.0/(g_size[i]+ 0.5*k2[i])) * (1.0+ s_l * (g_size[i]+ 0.5*k2[i])/r);
    k4[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 1.0*k3[i]) * (1.0/r - 1.0/(g_size[i]+ 1.0*k3[i])) * (1.0+ s_l * (g_size[i]+ 1.0*k3[i])/r);
    g_vary[i] = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
   }
  }
 }
 else if (modes == 'a')           // abnormal grain growth
 {
  double dg, g_vn, g_va;
  for (int i=0; i<MAX; i++)
  {
   if (g_size[i] == 0 || g_size[i] == r)
    g_vary[i] = 0.0;
   else if (g_size[i] < r)
   {
    g_vary[i] = h_time * c_A * f_A * 1.0/g_size[i] * (1.0/r - 1.0/g_size[i]) * (1.0+ s_l * g_size[i]/r);
    if ((g_vary[i] + g_size[i]) < 0)
     g_vary[i] = -1.0 * g_size[i];
   }
   else
   {
    k1[i] = h_time * c_A * f_A * 1.0/g_size[i] * (1.0/r - 1.0/g_size[i]) * (1.0+ s_l * g_size[i]/r);
    k2[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 0.5*k1[i]) * (1.0/r - 1.0/(g_size[i]+ 0.5*k1[i])) * (1.0+ s_l * (g_size[i]+ 0.5*k1[i])/r);
    k3[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 0.5*k2[i]) * (1.0/r - 1.0/(g_size[i]+ 0.5*k2[i])) * (1.0+ s_l * (g_size[i]+ 0.5*k2[i])/r);
    k4[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 1.0*k3[i]) * (1.0/r - 1.0/(g_size[i]+ 1.0*k3[i])) * (1.0+ s_l * (g_size[i]+ 1.0*k3[i])/r);
    g_vn = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    k1[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/g_size[i]));
    k2[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/(g_size[i]+ 0.5*k1[i])));
    k3[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/(g_size[i]+ 0.5*k2[i])));
    k4[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/(g_size[i]+ 1.0*k3[i])));
    g_va = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    g_vary[i] = 1.0/(1.0/g_va + 1.0/g_vn);
   }
  }
 }
 else                             // abnormal grain growth
 {                                // with assistance of screw-disl.
  double dg, g_vn, g_va, g_vs;
  for (int i=0; i<MAX; i++)
  {
   if (g_size[i] == 0 || g_size[i] == r)
    g_vary[i] = 0.0;
   else if (g_size[i] < r)
   {
    g_vary[i] = h_time * c_A * f_A * 1.0/g_size[i] * (1.0/r - 1.0/g_size[i]) * (1.0+ s_l * g_size[i]/r);
    if ((g_vary[i] + g_size[i]) < 0)
     g_vary[i] = -1.0 * g_size[i];
   }
   else
   {
    k1[i] = h_time * c_A * f_A * 1.0/g_size[i] * (1.0/r - 1.0/g_size[i]) * (1.0+ s_l * g_size[i]/r);
    k2[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 0.5*k1[i]) * (1.0/r - 1.0/(g_size[i]+ 0.5*k1[i])) * (1.0+ s_l * (g_size[i]+ 0.5*k1[i])/r);
    k3[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 0.5*k2[i]) * (1.0/r - 1.0/(g_size[i]+ 0.5*k2[i])) * (1.0+ s_l * (g_size[i]+ 0.5*k2[i])/r);
    k4[i] = h_time * c_A * f_A * 1.0/(g_size[i]+ 1.0*k3[i]) * (1.0/r - 1.0/(g_size[i]+ 1.0*k3[i])) * (1.0+ s_l * (g_size[i]+ 1.0*k3[i])/r);
    g_vn = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    k1[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/g_size[i]));
    k2[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/(g_size[i]+ 0.5*k1[i])));
    k3[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/(g_size[i]+ 0.5*k2[i])));
    k4[i] = h_time * c_B * f_B * exp(-1.0* s_e * f_C * c_C /(1.0/r - 1.0/(g_size[i]+ 1.0*k3[i])));
    g_va = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    k1[i] = h_time * c_D * f_A * s_s * 1.0/g_size[i] * pow(1.0/r - 1.0/g_size[i],2.0);
    k2[i] = h_time * c_D * f_A * s_s * 1.0/(g_size[i]+ 0.5*k1[i]) * pow(1.0/r - 1.0/(g_size[i]+ 0.5*k1[i]),2.0);
    k3[i] = h_time * c_D * f_A * s_s * 1.0/(g_size[i]+ 0.5*k2[i]) * pow(1.0/r - 1.0/(g_size[i]+ 0.5*k2[i]),2.0);
    k4[i] = h_time * c_D * f_A * s_s * 1.0/(g_size[i]+ 1.0*k3[i]) * pow(1.0/r - 1.0/(g_size[i]+ 1.0*k3[i]),2.0);
    g_vs = (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i])/6.0;
    g_vary[i] = 1.0/(1.0/g_va + 1.0/g_vn) + g_vs;
    if (g_vary[i] > g_vn)
     g_vary[i] = g_vn;
   }
  }
 }
 delete [] k1;
 delete [] k2;
 delete [] k3;
 delete [] k4;
}

double Grain::mass(double r)
{
 double * gm_v;
 gm_v = new double[MAX];
 double v_mass = 0.0;
 Grain::vary(gm_v, r);
 for (int j=0; j<MAX; j++)
 {
  v_mass += c_M * g_size[j] * g_size[j] * gm_v[j];
 }
 delete [] gm_v;
 return v_mass;
}

//----------------------------------------------------------
