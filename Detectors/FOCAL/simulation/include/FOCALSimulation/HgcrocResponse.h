// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_FOCAL_HGCROCRESPONSE_H_
#define ALICEO2_FOCAL_HGCROCRESPONSE_H_

#include <array>
#include <string>
#include <cstdint>

namespace o2
{
namespace focal
{

class HgcrocResponse
{
 public:
  static constexpr int NASICS = 18, NCHANNELS = 72;
  static constexpr int ADCMAX = 1023; // 10-bit ADC

  // sanity windows for the calibration fits (from the test-beam analysis)
  static constexpr double ADCSLOPE_MIN = 2500., ADCSLOPE_MAX = 4500.;
  static constexpr double TOTSLOPE_MIN = 200., TOTSLOPE_MAX = 550.;
  static constexpr double TOTSLOPEHIGH_MIN = 150., TOTSLOPEHIGH_MAX = 300.;

  // container for the calibration parameters of each pad
  struct ChannelCalib {
    double adcConst = 0, adcSlope = 0, adcLowFit = 0, adcHighFit = 0;
    double totConst = 0, totSlope = 0, totLowFit = 0, totHighFit = 0;
    double totConstHigh = 0, totSlopeHigh = 0;

    double G = 0, Qlin = 0, D = 0, Qsat = 0, Qtot = 0; // derived params (computed from the above)

    bool seenLow = false, seenHigh = false; // calibration found in each file
    uint16_t badMask = 0;                   // which groups of fit parameters are bad (replaced by the half mean)
  };

  // bits used to set ChannelCalib::badMask, one per group of fit parameters
  enum CalibBits : uint16_t {
    kBadADC = 1 << 0,
    kBadTOTLow = 1 << 1,
    kBadTOTHigh = 1 << 2
  };

  void Init();
  const std::array<std::array<ChannelCalib, NCHANNELS>, NASICS>& getCalibParams() const { return mCalib; }

 private:
  void LoadCalibLow(const std::string& filepath);
  void LoadCalibHigh(const std::string& filepath);
  void ValidateCalibParams();
  void RegularizeCalibParams();
  void ComputeDerivedParams();

  std::array<std::array<ChannelCalib, NCHANNELS>, NASICS> mCalib;
  bool mInitialized = false;
};

} // namespace focal
} // namespace o2

#endif