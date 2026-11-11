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

#include "FOCALSimulation/HgcrocResponse.h"
#include "FOCALSimulation/PadSimParam.h"
#include <fairlogger/Logger.h>
#include <TSystem.h>
#include <TString.h>

#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <filesystem>

namespace o2::focal
{

void HgcrocResponse::Init()
{
  mInitialized = false;
  mCalib = {};

  const PadSimParam& pars = PadSimParam::Instance();

  TString lowInjCalibData = pars.getCalibFileLow();
  gSystem->ExpandPathName(lowInjCalibData);
  TString highInjCalibData = pars.getCalibFileHigh();
  gSystem->ExpandPathName(highInjCalibData);

  LoadCalibLow(lowInjCalibData.Data());
  LoadCalibHigh(highInjCalibData.Data());
  ValidateCalibParams();
  RegularizeCalibParams();

  mInitialized = true;
}

void HgcrocResponse::LoadCalibLow(const std::string& filepath)
{
  std::ifstream in(filepath);
  std::string line;

  if (!in) {
    LOG(fatal) << "focal::HgcrocResponse: could not read calibration file: " << filepath;
  }

  LOG(info) << "focal::HgcrocResponse: reading calibration file: " << filepath;

  int asic, ch;
  double adcLowFit, adcHighFit, adcConst, adcSlope, totLowFit, totHighFit, totConst, totSlope;
  int lineNum = 0;

  while (std::getline(in, line)) {
    ++lineNum;
    if (line.empty())
      continue;

    std::istringstream ss(line);

    if (!(ss >> asic >> ch >> adcLowFit >> adcHighFit >> adcConst >> adcSlope >> totLowFit >> totHighFit >> totConst >> totSlope)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":" << lineNum << ": error parsing line";
    }
    if ((asic < 0) || (asic > NASICS - 1) || (ch < 0) || (ch > NCHANNELS - 1)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":" << lineNum << ": out of range entry for ASIC " << asic << " channel " << ch;
    }
    if (mCalib[asic][ch].seenLow) {
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << " has received more than one calibration entry from " << std::filesystem::path(filepath).filename().string();
    }

    mCalib[asic][ch].adcLowFit = adcLowFit;
    mCalib[asic][ch].adcHighFit = adcHighFit;
    mCalib[asic][ch].adcConst = adcConst;
    mCalib[asic][ch].adcSlope = adcSlope;
    mCalib[asic][ch].totLowFit = totLowFit;
    mCalib[asic][ch].totHighFit = totHighFit;
    mCalib[asic][ch].totConst = totConst;
    mCalib[asic][ch].totSlope = totSlope;
    mCalib[asic][ch].seenLow = true;
  }
}

void HgcrocResponse::LoadCalibHigh(const std::string& filepath)
{
  std::ifstream in(filepath);
  std::string line;

  if (!in) {
    LOG(fatal) << "focal::HgcrocResponse: cannot read file " << filepath;
  }

  LOG(info) << "focal::HgcrocResponse: reading calibration file: " << filepath;

  int asic, ch;
  double totConstHigh, totSlopeHigh;
  int lineNum = 0;

  while (std::getline(in, line)) {
    ++lineNum;
    if (line.empty())
      continue;

    std::istringstream ss(line);

    if (!(ss >> asic >> ch >> totSlopeHigh >> totConstHigh)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":" << lineNum << ": error parsing line";
    }
    if ((asic < 0) || (asic > NASICS - 1) || (ch < 0) || (ch > NCHANNELS - 1)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":" << lineNum << ": out of range entry for ASIC " << asic << " channel " << ch;
    }
    if (mCalib[asic][ch].seenHigh) {
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << " has received more than one calibration entry from " << std::filesystem::path(filepath).filename().string();
    }

    mCalib[asic][ch].totConstHigh = totConstHigh;
    mCalib[asic][ch].totSlopeHigh = totSlopeHigh;
    mCalib[asic][ch].seenHigh = true;
  }
}

// validates whether mCalibration parameters are valid or nonsense and sets badMask accordingly
void HgcrocResponse::ValidateCalibParams()
{
  // sanity windows
  const double adcSlopeMin = 2500., adcSlopeMax = 4500.;
  const double totSlopeMin = 200., totSlopeMax = 550.;
  const double totSlopeHighMin = 150., totSlopeHighMax = 300.;

  for (int asic = 0; asic < NASICS; ++asic) {
    for (int ch = 0; ch < NCHANNELS; ++ch) {
      auto& c = mCalib[asic][ch];
      c.badMask = 0;

      if (!c.seenLow) {
        c.badMask |= kBadADC | kBadTOTLow;
        LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << " did not receive calibration parameters for low charge injection fit, half mean will be used";
      } else {
        const double D = 1023. - c.adcConst - c.adcSlope * c.adcHighFit;
        if (c.adcSlope < adcSlopeMin || c.adcSlope > adcSlopeMax || D <= 0) {
          c.badMask |= kBadADC;
          LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << ": bad ADC fit (slope = " << c.adcSlope << ", D = " << D << "), will be replaced by the half mean";
        }
        if (c.totSlope < totSlopeMin || c.totSlope > totSlopeMax) {
          c.badMask |= kBadTOTLow;
          LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << ": bad ToT low charge injection fit parameters (slope = " << c.totSlope << "), will be replaced by the half mean";
        }
      }

      if (!c.seenHigh) {
        c.badMask |= kBadTOTHigh;
        LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << " did not receive calibration parameters for high charge injection fit, half mean will be used";
      } else {
        if (c.totSlopeHigh < totSlopeHighMin || c.totSlopeHigh > totSlopeHighMax) {
          c.badMask |= kBadTOTHigh;
          LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch << ": bad ToT high charge injection fit parameters (slope = " << c.totSlopeHigh << "), will be replaced by the half mean";
        }
      }
    }
  }
}

// replace the flagged groups of bad fit parameters by the mean of the good channels in the same ASIC half
void HgcrocResponse::RegularizeCalibParams()
{
  using Param = double ChannelCalib::*;

  const std::vector<Param> adcParams{&ChannelCalib::adcLowFit, &ChannelCalib::adcHighFit, &ChannelCalib::adcConst, &ChannelCalib::adcSlope};
  const std::vector<Param> totLowParams{&ChannelCalib::totLowFit, &ChannelCalib::totHighFit, &ChannelCalib::totConst, &ChannelCalib::totSlope};
  const std::vector<Param> totHighParams{&ChannelCalib::totConstHigh, &ChannelCalib::totSlopeHigh};

  auto regularizeGroup = [this](int asic, int half, uint16_t bit, const std::vector<Param>& params) {
    const int firstChannel = half * NCHANNELS / 2;
    const int lastChannel = firstChannel + NCHANNELS / 2;

    int nGood = 0;
    std::vector<double> meanParams(params.size(), 0.);
    std::vector<int> badChannels;

    // calculating means over the good channels
    for (int ch = firstChannel; ch < lastChannel; ++ch) {
      const auto& c = mCalib[asic][ch];

      if (c.badMask & bit) {
        badChannels.push_back(ch);
        continue;
      }
      for (size_t i = 0; i < params.size(); ++i) {
        meanParams[i] += c.*params[i];
      }

      nGood++;
    }

    if (badChannels.empty())
      return; // nothing to replace
    if (nGood == 0) {
      std::string group = "unknown fit group";
      if (bit == kBadADC) {
        group = "ADC fit";
      } else if (bit == kBadTOTLow) {
        group = "TOT low charge injection fit";
      } else if (bit == kBadTOTHigh) {
        group = "TOT high charge injection fit";
      }
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " half " << half << " has no good channels to average for the " << group;
    }

    for (double& mean : meanParams) {
      mean /= nGood;
    }

    // replacing bad parameter with the mean for the ASIC half
    for (int ch : badChannels) {
      auto& c = mCalib[asic][ch];
      for (size_t i = 0; i < params.size(); ++i) {
        c.*params[i] = meanParams[i];
      }
      if (bit == kBadADC) { // the half mean does not guarantee D > 0, so check it again after the substitution
        const double D = 1023. - c.adcConst - c.adcSlope * c.adcHighFit;
        if (D <= 0) {
          LOG(fatal) << "focal::HgcrocResponse: the ADC fit parameters of ASIC " << asic << " half " << half << " are such that the dynamic range of the ADC is D <= 0 (D = " << D << ") after substitution of invalid fit parameters by half mean";
        }
      }
    }
  };

  for (int asic = 0; asic < NASICS; ++asic) {
    for (int half = 0; half < 2; ++half) {
      regularizeGroup(asic, half, kBadADC, adcParams);
      regularizeGroup(asic, half, kBadTOTLow, totLowParams);
      regularizeGroup(asic, half, kBadTOTHigh, totHighParams);
    }
  }
}

} // namespace o2::focal