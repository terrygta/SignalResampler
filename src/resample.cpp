#include <boost/math/special_functions/bessel.hpp>

#include <vector>
#include <math.h>
#include "upfirdn.h"
#include "resample.h"

using namespace std;

int getGCD ( int num1, int num2 )
{
  int tmp = 0;
  while ( num1 > 0 )
  {
    tmp = num1;
    num1 = num2 % num1;
    num2 = tmp;
  }
  return num2;
}

int quotientCeil ( int num1, int num2 )
{
  if ( num1 % num2 != 0 )
    return num1 / num2 + 1;
  return num1 / num2;
}

double sinc ( double x )
{
  if ( fabs ( x - 0.0 ) < 0.000001 )
    return 1;
  return sin ( M_PI * x ) / ( M_PI * x );
}

void firls ( int length, vector<double> freq, 
  const vector<double>& amplitude, vector<double>& result )
{
  vector<double> weight;
  int freqSize = freq.size ();
  int weightSize = freqSize / 2;

  weight.reserve ( weightSize );
  for ( int i = 0; i < weightSize; i++ )
    weight.push_back ( 1.0 );

  int filterLength = length + 1;

  for ( int i = 0; i < freqSize; i++ )
    freq[i] /= 2.0;

  vector<double> dFreq;
  for ( int i = 1; i < freqSize; i++ )
    dFreq.push_back ( freq[i] - freq[i - 1] );

  length = ( filterLength - 1 ) / 2;
  int Nodd = filterLength % 2;
  double b0 = 0.0;
  vector<double> k;
  if ( Nodd == 0 )
  {
    for ( int i = 0; i <= length; i++)
      k.push_back ( i + 0.5 );
  }
  else
  {
    for ( int i = 0; i <= length; i++)
      k.push_back ( i );
  }

  vector<double> b;
  int kSize = k.size();
  for ( int i = 0; i < kSize; i++ )
    b.push_back( 0.0 );
  for ( int i = 0; i < freqSize; i += 2 )
  {
    double slope = ( amplitude[i + 1] - amplitude[i] ) / ( freq[i + 1] - freq[i] );
    double b1 = amplitude[i] - slope * freq[i];
    if ( Nodd == 1 )
    {
      b0 += ( b1 * ( freq[i + 1] - freq[i] ) ) + 
        slope / 2.0 * ( freq[i + 1] * freq[i + 1] - freq[i] * freq[i] ) * 
          fabs( weight[(i + 1) / 2] * weight[(i + 1) / 2] );
    }
    for ( int j = 0; j < kSize; j++ )
    {
      b[j] += ( slope / ( 4 * M_PI * M_PI ) * 
        (cos ( 2 * M_PI * k[j] * freq[i + 1] ) - cos ( 2 * M_PI * k[j] * freq[i] )) / ( k[j] * k[j] )) *
          fabs( weight[(i + 1) / 2] * weight[(i + 1) / 2] );
      b[j] += ( freq[i + 1] * ( slope * freq[i + 1] + b1 ) * sinc( 2 * k[j] * freq[i + 1] ) - 
        freq[i] * ( slope * freq[i] + b1 ) * sinc( 2 * k[j] * freq[i] ) ) *
          fabs( weight[(i + 1) / 2] * weight[(i + 1) / 2] );
    }
  }
  if ( Nodd == 1 )
    b[0] = b0;
  vector<double> a;
  double w0 = weight[0];
  for ( int i = 0; i < kSize; i++ )
    a.push_back(( w0 * w0 ) * 4 * b[i]);
  if ( Nodd == 1 )
  {
    a[0] /= 2;
    for ( int i = length; i >= 1; i-- )
      result.push_back( a[i] / 2.0 );
    result.push_back( a[0] );
    for ( int i = 1; i <= length; i++ )
      result.push_back( a[i] / 2.0 );
  }
  else
  {
    for ( int i = length; i >= 0; i-- )
      result.push_back( a[i] / 2.0 );
    for ( int i = 0; i <= length; i++ )
      result.push_back( a[i] / 2.0 );
  }
}

void kaiser ( const int order, const double bta, vector<double>& window )
{
  double bes = fabs ( boost::math::cyl_bessel_i ( 0, bta ) );
  int odd = order % 2;
  double xind = ( order - 1 ) * ( order - 1 );
  int n = ( order + 1 ) / 2;
  vector<double> xi;
  xi.reserve ( n );
  for ( int i = 0; i < n; i++ )
  {
    double val = static_cast<double>( i ) + 0.5 * ( 1 - static_cast<double>( odd ) );
    xi.push_back ( 4 * val * val );
  }
  vector<double> w;
  w.reserve ( n );
  for ( int i = 0; i < n; i++ )
    w.push_back ( boost::math::cyl_bessel_i ( 0, bta * sqrt( 1 - xi[i] / xind ) ) / bes );
  for ( int i = n - 1; i >= odd; i-- )
    window.push_back ( fabs ( w[i] ) );
  for ( int i = 0; i < n; i++ )
    window.push_back ( fabs ( w[i] ) );
}

void resample ( int upFactor, int downFactor, 
  vector<double>& inputSignal, vector<double>& outputSignal )
{
  const int n = 10;
  const double bta = 5.0;
  if ( upFactor <= 0 || downFactor <= 0 )
    throw std::runtime_error ( "factors must be positive integer" );
  int gcd = getGCD ( upFactor, downFactor );
  upFactor /= gcd;
  downFactor /= gcd;

  if ( upFactor == downFactor )
  {
    outputSignal = inputSignal;
    return;
  }

  int inputSize = inputSignal.size();
  outputSignal.clear ();
  int outputSize =  quotientCeil ( inputSize * upFactor, downFactor );
  outputSignal.reserve ( outputSize );

  int maxFactor = max ( upFactor, downFactor );
  double firlsFreq = 1.0 / 2.0 / static_cast<double> ( maxFactor );
  int length = 2 * n * maxFactor + 1;
  double firlsFreqs[] = { 0.0, 2.0 * firlsFreq, 2.0 * firlsFreq, 1.0 };
  vector<double> firlsFreqsV;
  firlsFreqsV.assign ( firlsFreqs, firlsFreqs + 4 );
  double firlsAmplitude[] = { 1.0, 1.0, 0.0, 0.0 };
  vector<double> firlsAmplitudeV;
  firlsAmplitudeV.assign ( firlsAmplitude, firlsAmplitude + 4 );
  vector<double> coefficients;
  firls ( length - 1, firlsFreqsV, firlsAmplitudeV, coefficients );
  vector<double> window;
  kaiser ( length, bta, window );
  int coefficientsSize = coefficients.size();
  for( int i = 0; i < coefficientsSize; i++ )
    coefficients[i] *= upFactor * window[i];

  int lengthHalf = ( length - 1 ) / 2;
  int nz = downFactor - lengthHalf % downFactor;
  vector<double> h;
  h.reserve ( coefficientsSize + nz );
  for ( int i = 0; i < nz; i++ )
    h.push_back ( 0.0 );
  for ( int i = 0; i < coefficientsSize; i++ )
    h.push_back ( coefficients[i] );
  int hSize = h.size();
  lengthHalf += nz;
  int delay = lengthHalf / downFactor;
  nz = 0;
  while ( quotientCeil( ( inputSize - 1 ) * upFactor + hSize + nz, downFactor ) - delay < outputSize )
    nz++;
  for ( int i = 0; i < nz; i++ )
    h.push_back ( 0.0 );
  vector<double> y;
  upfirdn ( upFactor, downFactor, inputSignal, h, y );
  for ( int i = delay; i < outputSize + delay; i++ )
  {
    outputSignal.push_back ( y[i] );
  }
}