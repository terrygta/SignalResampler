#include <iostream>
#include <vector>
#include <chrono>
#include "resample.h"

using namespace std;

int main ()
{
  const int kUpfactor = 2048;
  const int kDownfactor = 100;

  vector<double> input(10), output;
  std::iota(input.begin(), input.end(), 1.0);
  auto start = std::chrono::high_resolution_clock::now();
  resample<double> ( kUpfactor, kDownfactor, input, output );
  auto stop = std::chrono::high_resolution_clock::now();
  cout<<"double: ";
  for ( int i = 0, n = output.size(); i < n; i++ )
    cout << output[i] << " ";
  cout<<endl;
  auto dur_d = std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
  cout<<"double elapsed: "<< dur_d.count()/1e9l<<" secs";
  cout<<endl<<endl;

  vector<float> input_f(10), output_f;
  std::iota(input_f.begin(), input_f.end(), 1.0);
  start = std::chrono::high_resolution_clock::now();
  resample<float> ( kUpfactor, kDownfactor, input_f, output_f );
  stop = std::chrono::high_resolution_clock::now();
  cout<<"float : ";
  for ( int i = 0, n = output_f.size(); i < n; i++ )
    cout << output_f[i] << " ";
  cout<<endl;
  auto dur_f = std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start);
  cout<<"float  elapsed: "<< dur_f.count()/1e9l<<" secs";
  cout<<endl<<endl;

  return 0;
}
