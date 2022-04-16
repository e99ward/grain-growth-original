// define CONVERT function 
#include <iostream>
#include <sstream>
#include <string>
using namespace std;


template <typename T>
string to_string(T t, ios_base & (*f)(ios_base&))
{
 ostringstream oss;
 oss << f << t;
 return oss.str();
}
// numeric type to string                        //
//-----------------------------------------------//
// cout << to_string<long>(123456, hex) << endl; //
// cout << to_string<double>(12.3, oct) << endl; //


template <typename T>
bool from_string(T &t, const string &s, ios_base & (*f)(ios_base&))
{
 istringstream iss(s);
 return !(iss >> f >> t).fail();
}
// string to numeric type                        //
//-----------------------------------------------//
// if (from_string<int>(i, string("ff"), hex))   //
// { cout << i << endl;}  f=float		 //


//----------------------------------------------------------
