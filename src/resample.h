//RESAMPLE  Change the sampling rate of a signal.
//   Y = RESAMPLE(UpFactor, DownFactor, InputSignal, OutputSignal) resamples the sequence in 
//   vector InputSignal at UpFactor/DownFactor times and stores the resampled data to OutputSignal.
//   OutputSignal is UpFactor/DownFactor times the length of InputSignal. UpFactor and DownFactor must be 
//   positive integers.

//This function is translated from Matlab's Resample funtion. 

//Author: Haoqi Bai

#pragma once

#include <boost/math/special_functions/bessel.hpp>

#include <vector>

using namespace std;

void resample ( int upFactor, int downFactor, 
  vector<double>& inputSignal, vector<double>& outputSignal );