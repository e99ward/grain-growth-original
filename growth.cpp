// main function using GRAIN class
#include <iostream>
#include <fstream>
using namespace std;
#include <ctime>
#include "grain.h"
#include "constant.h"
using constant::NGG;
using constant::AGG;

void growth_help();
int data_init();
int histo_bin();

int main(int argc, char * argv[])
{
 cout << "+___________________________________________+" << endl;
 cout << "|                                           |" << endl;
 cout << "|    Grain Growth Calculator  < C >         |" << endl;
 cout << "|                                           |" << endl;
 cout << "|                 by Edward Yang-il Jung    |" << endl;
 cout << "|___________________________________________|" << endl;

 bool mode = NGG;               // NGG or AGG
 bool m_spiral = false;         // screw-disl. assisted
 int v_temper = 1500;           // 1500 celsius
 double v_liquid = 0.46;        // 0.46 (random packing)
 double v_stepfree = 0.33;      // 0.33 hGamma
 double m_speed = 0.1;          // calculation step (CTS)

 if (argc>1)
 {
  for (int i=1; i<argc; i++)
  {
   if (argv[i][0] == '-')
   {
    if (argv[i][1] == 'm')        // setting MODE
    {
     if (toupper(argv[i+1][0]) == 'A')
      mode = AGG;
     else
      mode = NGG;
    }
    else if (argv[i][1] == 'w')   // setting SCREW-DISL.
    {
     m_spiral = true;
    }
    else if (argv[i][1] == 't')   // setting TEMPERATURE
     v_temper = atoi(argv[i+1]);
    else if (argv[i][1] == 'l')   // setting LIQUID FRACTION
     v_liquid = atof(argv[i+1]);
    else if (argv[i][1] == 's')   // setting STEP FREE ENERGY
     v_stepfree = atof(argv[i+1]);
    else if (argv[i][1] == 'q')   // setting Faster
     m_speed = 1.0;
    else if (argv[i][1] == 'z')   // setting Slower
     m_speed = 0.01;
    else if (argv[i][1] == 'i')   // generate new INITIAL set
     exit(data_init());
    else if (argv[i][1] == 'k')   // generate HISTOGRAM bins
     exit(histo_bin());
    else                          // show HELP
    {
     growth_help();
     exit(0);
    }
   }
  }
 }
 else
 {
  growth_help();
  exit(0);
 }

 cout << "****** configure the simulation ******" << endl;
 cout << "---------------------------------------------" << endl;
 int InitialStep = 0;
 int FinalStep = 5000;
 int SaveModG = 100;
 int SaveModH = 100;
 char temp[10];
 // setting INITIAL TIME STEP and load the FILE
 cout << "initial step [0]:_____\b\b\b\b\b";
 cin.getline(temp,10);
 char file[15] = "d_0000000.txt";
 if (strlen(temp)!=0)
 {
  InitialStep = atoi(temp);
  int h = strlen(temp) - 1;
  int k = 8;
  while(h>=0)
   file[k--] = temp[h--];
 }
 ifstream fin;                  // can load initial file?
 fin.open(file,ios_base::in);
 if (!fin.is_open())
 {
  cerr << "can't open the file: " << file << endl;
  cerr << "restart the program" << endl;
  fin.close();
  exit(1);
 }
 fin.close();
 // setting FINAL TIME STEP
 cout << "final step [5000]:_____\b\b\b\b\b";
 cin.getline(temp,10);
 if (strlen(temp)!=0)
  FinalStep = atoi(temp);
 // setting SAVE TIME STEP for Grain
 cout << "save step for grain data [100]:_____\b\b\b\b\b";
 cin.getline(temp,10);
 if (strlen(temp)!=0)
  SaveModG = atoi(temp);
 // setting SAVE TIME STEP for Histogram
 cout << "save step for histograms [100]:_____\b\b\b\b\b";
 cin.getline(temp,10);
 if (strlen(temp)!=0)
  SaveModH = atoi(temp);

 cout << "---------------------------------------------" << endl;
 if (mode==AGG)
 {
  cout << "AGG simulation " << endl;
  if (m_spiral==true)
   cout << "with screw-disl. assisted growth" << endl;
 }
 else
 {
  cout << "NGG simulation " << endl;
  v_stepfree = 0.0;
  if (m_spiral==true)
   cout << "WARNING: NGG cannot incorporate spiral growth" << endl;
 }
 cout << "initial_step=" << InitialStep
      << ", finial_step=" << FinalStep << endl;
 cout << "data_save_step=" << SaveModG
      << ", histogram_step=" << SaveModH << endl;
 cout << "temperature is " << v_temper << "C" << endl;
 cout << "step free energy is " << v_stepfree << " hG" << endl;
 cout << "liquid volume fraction is " << v_liquid << endl;
 cout << "calculation speed is " << m_speed << " sec/CTS" << endl;
 cout << "---------------------------------------------" << endl;
 cout << "** continue? (Y)es, ENTER! / (Q)uit, ^C" << flush;
 char ch = cin.get();
 if (toupper(ch)=='Q')
  exit(2);
 
 time_t time1, time2;
 time(&time1);
 cout << "---------------------------------------------" << endl;
 cout << "****** start ******" << endl;
 cout << "****** " << ctime(&time1);
 cout << "---------------------------------------------" << endl;
 int cts = InitialStep;
 cout << "read from file: step " << cts << endl;
 Grain set = Grain(cts);                                //initialize
 set.spec(mode, m_spiral, v_temper, v_stepfree, v_liquid, m_speed); //specify
 set.writeStat(cts);
 set.writeHist(cts);
 cout << "now calculation is going on..." << endl;
 while (cts<FinalStep)
 {
  set.calc();                                   // calculation start
  cts++;
  cout << "calculation: step" << cts << endl;
  if ((cts%SaveModH)==0)                       // save the histogram
   set.writeHist(cts);
  if ((cts%SaveModG)==0)                       // save the grain data
   set.writeData(cts);
  set.writeStat(cts);                          // save the statistics
 }
 time(&time2); 
 cout << "---------------------------------------------" << endl;
 cout << "****** end ******" << endl;
 cout << "****** " << ctime(&time2);
 cout << "---------------------------------------------" << endl;
 
 int hour, mins;
 double secs = difftime(time2,time1);
 hour = int(secs/3600);
 secs -= 3600*hour;
 mins = int(secs/60);
 secs -= 60*mins;
 cout << "** times consumed ";
 cout.width(2);
 cout << hour << ":";
 cout.width(2);
 cout << mins << ":";
 cout.width(2);
 cout << secs << " hours **" << endl;
 cout << "---------------------------------------------" << endl;

 return 0;
}
