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
#include "Rtypes.h"

namespace o2{
namespace focal{

class HgcrocResponse{
 public:
  static constexpr int NASICS = 18, NCHANNELS = 72;

  // Container for the calibration parameters of each pad
  struct ChannelCalib {
  // raw, from the files
  double adcConst = 0, adcSlope = 0, adcLowFit = 0, adcHighFit = 0;
  double totConst = 0, totSlope = 0, totLowFit = 0, totHighFit = 0;
  double totConstHigh = 0, totSlopeHigh = 0;
  // derived, computed in init()
  double G = 0, Qlin = 0, D = 0, Qsat = 0, Qtot = 0;
  // bookkeeping:
  bool seenLow = false, seenHigh = false;   // line found in each file
  uint16_t replacedMask = 0;                // which constants were replaced by the half mean
  };

  void Init();

 private:
  std::array<std::array<ChannelCalib, NCHANNELS>, NASICS> mCalib;
  bool mInitialized = false;
};

} // namespace focal
} // namespace o2

#endif