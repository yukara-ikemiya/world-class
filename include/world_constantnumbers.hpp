//-----------------------------------------------------------------------------
// Copyright 2012 Masanori Morise,
// Copyright 2018 Yukara Ikemiya
//
// This header file only defines constant numbers used for several function.
//-----------------------------------------------------------------------------
#ifndef WORLD_CONSTANT_NUMBERS_HPP
#define WORLD_CONSTANT_NUMBERS_HPP

namespace world {

constexpr double kPi = 3.1415926535897932384;
constexpr double kMySafeGuardMinimum = 0.000000000001;
constexpr double kEps = 0.00000000000000022204460492503131;
constexpr double kFloorF0 = 71.0;
constexpr double kCeilF0 = 800.0;
constexpr double kDefaultF0 = 500.0;
constexpr double kLog2 = 0.69314718055994529;
// Maximum standard deviation not to be selected as a best f0.
constexpr double kMaximumValue = 100000.0;
// Note to me (fs: 48000)
// 71 Hz is the limit to maintain the FFT size at 2048.
// If we use 70 Hz as FLOOR_F0, the FFT size of 4096 is required.

// for D4C()
constexpr int kHanning = 1;
constexpr int kBlackman = 2;
constexpr double kFrequencyInterval = 3000.0;
constexpr double kUpperLimit = 15000.0;
constexpr double kThreshold = 0.85;
constexpr double kFloorF0D4C = 47.0;

// for Codec (Mel scale)
// S. Stevens & J. Volkmann,
// The Relation of Pitch to Frequency: A Revised Scale,
// American Journal of Psychology, vol. 53, no. 3, pp. 329-353, 1940.
constexpr double kM0 = 1127.01048;
constexpr double kF0 = 700.0;
constexpr double kFloorFrequency = 40.0;
constexpr double kCeilFrequency = 20000.0;

}  // namespace world

#endif
