// define GRAIN class
#ifndef GRAIN_H_
#define GRAIN_H_

class Grain
{
 private:
  int t_grain;           // number of grain in system
  int MAX;               // number of grain in calculation
  char modes;            // calculation modes
  double h_time;         // calculation step
  double * g_size;       // grain Radius [MAX]
  double rs;             // critical grain size also in Radius
  double rav;            // average Radius
  double rst;            // standard deviation
  double s_l;            // liquid volume fraction
  double s_e;            // step free energy (0.62hG)
  double s_s;            // step free energy (0.62hG)
  double f_A;            // temperature on diffusion (1500C=1773K)
  double f_B;            // temperature on diffusion (1500C=1773K)
  double f_C;            // temperature on exponants (1500C=1773K)
 public:
  Grain();
  Grain(int step);
  ~Grain();
  void spec(bool mode1, bool mode2, int temp, double sfe, double liq, double speed);
  void calc();               // calc diffusion
  void writeData(int step);  // save grain data to file (step)
  void writeHist(int step);  // save histograms to file (step)
  void writeStat(int step);  // save Av, St (every step)
 private:
  void avst();		         // average size, standard deviation
  void rstar();              // critical grain size
  void vary(double g_vary[], double r);  // (dr/dt)
  double mass(double r);                 // r^2(dr/dt)
  void sort(double num[], int left, int right);
};

#endif
