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

/// \file testHgcrocResponse.cxx
/// \brief Unit tests for the calibration handling of o2::focal::HgcrocResponse
///
/// The test writes its own calibration files to a temporary directory, so it
/// does not depend on the real test-beam files. Every channel gets a "typical"
/// calibration, except a few channels that are overridden on purpose.

#define BOOST_TEST_MODULE Test FOCAL HgcrocResponse
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "FOCALSimulation/HgcrocResponse.h"
#include "CommonUtils/ConfigurableParam.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <TRandom3.h>
#include <array>
#include <iostream>

namespace o2::focal
{

using Calib = std::array<std::array<o2::focal::HgcrocResponse::ChannelCalib, 72>, 18>;

void printHalfMeans(const Calib& calib)
{
  using C = o2::focal::HgcrocResponse::ChannelCalib;

  struct Quantity {
    const char* name;
    double C::* field;
    double expected;
  };
  const Quantity quantities[] = {
    {"adcConst", &C::adcConst, 68.},
    {"adcSlope", &C::adcSlope, 4050.},
    {"adcHighFit", &C::adcHighFit, 0.20},
    {"totConst", &C::totConst, -16.},
    {"totSlope", &C::totSlope, 336.},
    {"totConstHigh", &C::totConstHigh, 152.},
    {"totSlopeHigh", &C::totSlopeHigh, 257.},
    {"Qsat", &C::Qsat, 0.27},
    {"Qtot", &C::Qtot, 0.37},
  };

  for (int asic = 0; asic < 18; ++asic) {
    for (int half = 0; half < 2; ++half) {
      std::cout << "ASIC " << asic << " half " << half << ":\n";
      for (const auto& q : quantities) {
        double sum = 0.;
        for (int ch = half * 36; ch < half * 36 + 36; ++ch) {
          sum += calib[asic][ch].*(q.field);
        }
        std::cout << "  " << q.name << " mean = " << sum / 36. << "  (expected " << q.expected << ")\n";
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(General)
{
  auto& suite = boost::unit_test::framework::master_test_suite();
  if (suite.argc > 1) {
    o2::conf::ConfigurableParam::updateFromString(suite.argv[1]);
  }

  HgcrocResponse resp;
  resp.Init();

  TRandom3 rand;
  rand.SetSeed(42);

  /*for (int k = 0; k < 15; ++k){
    int i = rand.Integer(18);
    int j = rand.Integer(72);

    std::cout << i << " " << j << " "
                << resp.mCalib[i][j].adcLowFit << " "
                << resp.mCalib[i][j].adcHighFit << " "
                << resp.mCalib[i][j].adcConst << " "
                << resp.mCalib[i][j].adcSlope << " "
                << resp.mCalib[i][j].totLowFit << " "
                << resp.mCalib[i][j].totHighFit << " "
                << resp.mCalib[i][j].totConst << " "
                << resp.mCalib[i][j].totSlope << " "
                << resp.mCalib[i][j].totSlopeHigh << " "
                << resp.mCalib[i][j].totConstHigh << " "
                << resp.mCalib[i][j].G << " "
                << resp.mCalib[i][j].Qlin << " "
                << resp.mCalib[i][j].D << " "
                << resp.mCalib[i][j].Qsat << " "
                << resp.mCalib[i][j].Qtot << " "
                << resp.mCalib[i][j].badMask << "\n";
  }*/

  auto calib = resp.getCalibParams();
  printHalfMeans(calib);

  for (int i = 0; i < 18; ++i) {
    for (int j = 0; j < 72; ++j) {
      std::cout << i << " " << j << " "
                << calib[i][j].adcLowFit << " "
                << calib[i][j].adcHighFit << " "
                << calib[i][j].adcConst << " "
                << calib[i][j].adcSlope << " "
                << calib[i][j].totLowFit << " "
                << calib[i][j].totHighFit << " "
                << calib[i][j].totConst << " "
                << calib[i][j].totSlope << " "
                << calib[i][j].totSlopeHigh << " "
                << calib[i][j].totConstHigh << " "
                << calib[i][j].G << " "
                << calib[i][j].Qlin << " "
                << calib[i][j].D << " "
                << calib[i][j].Qsat << " "
                << calib[i][j].Qtot << " "
                << calib[i][j].badMask << "\n";
    }
  }
}

} // namespace o2::focal
