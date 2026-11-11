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

} // namespace o2::focal