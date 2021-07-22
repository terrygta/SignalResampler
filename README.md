## Signal Resampler for C++

A C++ version of Matlab's resample function (y = resample(x, p, q), see http://www.mathworks.com.au/help/signal/ref/resample.html). The implementation of upfirdn function is based on Motorola's open source project (see http://sourceforge.net/motorola/upfirdn/home/Home/).

To compile:
- install g++ (g++7 or newer) and cmake
    - Ubuntu/Debian: *sudo apt install g++ cmake*
- *cmake .*
- *make*

To run example:
- *./example*
