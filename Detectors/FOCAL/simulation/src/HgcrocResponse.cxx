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
#include <filesystem>

namespace o2
{
namespace focal
{

template <typename SimParam>
void HgcrocResponse::init(const SimParam& pars)
{
  mInitialized = false;

  // loading config from SimParam class for pads or HCAL
  mConfig.mNAsics = pars.mNAsics;
  mConfig.mNChannels = pars.mNChannels;
  mConfig.mADCMax = pars.mADCMax;
  mConfig.mADCSlopeMin = pars.mADCSlopeMin;
  mConfig.mADCSlopeMax = pars.mADCSlopeMax;
  mConfig.mTOTSlopeMin = pars.mTOTSlopeMin;
  mConfig.mTOTSlopeMax = pars.mTOTSlopeMax;
  mConfig.mTOTSlopeHighRangeMin = pars.mTOTSlopeHighRangeMin;
  mConfig.mTOTSlopeHighRangeMax = pars.mTOTSlopeHighRangeMax;
  mConfig.mTOTLowFitMin = pars.mTOTLowFitMin;
  mConfig.mTOTHighRangeDefaultSlope = pars.mTOTHighRangeDefaultSlope;
  mConfig.mTOTHighRangeDefaultConst = pars.mTOTHighRangeDefaultConst;
  mConfig.mCalibFileLowRange = pars.mCalibFileLowRange;
  mConfig.mCalibFileHighRange = pars.mCalibFileHighRange;
  mConfig.mTOTBadHalfFraction = pars.mTOTBadHalfFraction;

  if (mConfig.mNAsics <= 0 || mConfig.mNChannels <= 0) {
    LOG(fatal) << "focal::HgcrocResponse: number of ASICs and channels must be positive";
  }
  if (mConfig.mNChannels % 2 != 0) {
    LOG(fatal) << "focal::HgcrocResponse: number of channels per ASIC must be even (the ASIC is split in two halves)";
  }

  mCalibParams.assign(mConfig.mNAsics, std::vector<ChannelCalib>(mConfig.mNChannels));

  TString lowRangeCalibData = mConfig.mCalibFileLowRange;
  gSystem->ExpandPathName(lowRangeCalibData);
  TString highRangeCalibData = mConfig.mCalibFileHighRange;
  gSystem->ExpandPathName(highRangeCalibData);

  loadCalibLowRange(lowRangeCalibData.Data());
  loadCalibHighRange(highRangeCalibData.Data());
  validateCalibParams();
  flagDeadTOTHalves();
  regularizeBadCalibParams();
  computeDerivedParams();

  mInitialized = true;
}

template void HgcrocResponse::init<PadSimParam>(const PadSimParam&);
// add HCAL here once it provides the same getters

void HgcrocResponse::loadCalibLowRange(const std::string& filepath)
{
  std::ifstream in(filepath);
  std::string line;

  if (!in) {
    LOG(fatal) << "focal::HgcrocResponse: cannot read calibration file: " << filepath;
  }

  LOG(debug) << "focal::HgcrocResponse: reading calibration file: " << filepath;

  int asic, ch;
  double adcLowerFitRange, adcHigherFitRange, adcConst, adcSlope, totLowerFitRange, totHigherFitRange, totConst, totSlope;
  int lineNum = 0;

  while (std::getline(in, line)) {
    ++lineNum;
    if (line.empty()) {
      continue;
    }

    std::istringstream ss(line);

    if (!(ss >> asic >> ch >> adcLowerFitRange >> adcHigherFitRange >> adcConst >> adcSlope >> totLowerFitRange >> totHigherFitRange >> totConst >> totSlope)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": error parsing line";
    }
    if ((asic < 0) || (asic > mConfig.mNAsics - 1) || (ch < 0) || (ch > mConfig.mNChannels - 1)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": out of range entry for ASIC " << asic << " channel " << ch;
    }
    if (mCalibParams[asic][ch].seenLow) {
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                 << " has received more than one calibration entry from " << std::filesystem::path(filepath).filename().string();
    }

    mCalibParams[asic][ch].adcLowerFitRange = adcLowerFitRange;
    mCalibParams[asic][ch].adcHigherFitRange = adcHigherFitRange;
    mCalibParams[asic][ch].adcConst = adcConst;
    mCalibParams[asic][ch].adcSlope = adcSlope;
    mCalibParams[asic][ch].totLowerFitRange = totLowerFitRange;
    mCalibParams[asic][ch].totHigherFitRange = totHigherFitRange;
    mCalibParams[asic][ch].totConst = totConst;
    mCalibParams[asic][ch].totSlope = totSlope;
    mCalibParams[asic][ch].seenLow = true;
  }
}

void HgcrocResponse::loadCalibHighRange(const std::string& filepath)
{
  std::ifstream in(filepath);
  std::string line;

  if (!in) {
    LOG(fatal) << "focal::HgcrocResponse: cannot read calibration file: " << filepath;
  }

  LOG(debug) << "focal::HgcrocResponse: reading calibration file: " << filepath;

  int asic, ch;
  double totConstHigh, totSlopeHigh;
  int lineNum = 0;

  while (std::getline(in, line)) {
    ++lineNum;
    if (line.empty()) {
      continue;
    }

    std::istringstream ss(line);

    if (!(ss >> asic >> ch >> totSlopeHigh >> totConstHigh)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": error parsing line";
    }
    if ((asic < 0) || (asic > mConfig.mNAsics - 1) || (ch < 0) || (ch > mConfig.mNChannels - 1)) {
      LOG(fatal) << "focal::HgcrocResponse: " << std::filesystem::path(filepath).filename().string() << ":"
                 << lineNum << ": out of range entry for ASIC " << asic << " channel " << ch;
    }
    if (mCalibParams[asic][ch].seenHigh) {
      LOG(fatal) << "focal::HgcrocResponse: ASIC " << asic << " channel " << ch
                 << " has received more than one calibration entry from " << std::filesystem::path(filepath).filename().string();
    }

    mCalibParams[asic][ch].totConstHigh = totConstHigh;
    mCalibParams[asic][ch].totSlopeHigh = totSlopeHigh;
    mCalibParams[asic][ch].seenHigh = true;
  }
}

// validates whether calibration parameters are valid or nonsense and sets badFitMask accordingly
void HgcrocResponse::validateCalibParams()
{
  int nMissingLow = 0, nMissingHigh = 0, nBadADC = 0, nBadTOTLow = 0, nBadTOTHigh = 0;
  std::ostringstream logMissingLow, logMissingHigh, logBadADC, logBadTOTLow, logBadTOTHigh;

  for (int asic = 0; asic < mConfig.mNAsics; ++asic) {
    for (int ch = 0; ch < mConfig.mNChannels; ++ch) {
      auto& c = mCalibParams[asic][ch];
      c.badFitMask = 0;

      if (!c.seenLow) {
        c.badFitMask |= BadADCFit | BadTOTLowFit;
        nMissingLow++;
        logMissingLow << " (" << asic << ", " << ch << ")";
      } else {
        const double adcMargin = mConfig.mADCMax - c.adcConst - c.adcSlope * c.adcHigherFitRange;
        const bool badADCSlope = c.adcSlope < mConfig.mADCSlopeMin || c.adcSlope > mConfig.mADCSlopeMax;
        const bool badADCMargin = adcMargin <= 0;
        const bool badTOTSlope = c.totSlope < mConfig.mTOTSlopeMin || c.totSlope > mConfig.mTOTSlopeMax;
        const bool badTOTFitRange = c.totLowerFitRange < mConfig.mTOTLowFitMin; // ToT never turned on

        if (badADCSlope || badADCMargin) {
          c.badFitMask |= BadADCFit;
          nBadADC++;
          logBadADC << " (" << asic << ", " << ch << ")";
        }
        if (badTOTSlope || badTOTFitRange) {
          c.badFitMask |= BadTOTLowFit;
          nBadTOTLow++;
          logBadTOTLow << " (" << asic << ", " << ch << ")";
        }
      }

      if (!c.seenHigh) {
        c.badFitMask |= BadTOTHighFit;
        nMissingHigh++;
        logMissingHigh << " (" << asic << ", " << ch << ")";
      } else {
        const bool badTOTSlopeHigh = c.totSlopeHigh < mConfig.mTOTSlopeHighRangeMin || c.totSlopeHigh > mConfig.mTOTSlopeHighRangeMax;
        // in test-beam 2023, default parameters = bad-channel, so we flag it:
        const bool hasDefaultSlope = c.totSlopeHigh == mConfig.mTOTHighRangeDefaultSlope;
        const bool hasDefaultConst = c.totConstHigh == mConfig.mTOTHighRangeDefaultConst;
        const bool hasDefaultParameters = hasDefaultSlope && hasDefaultConst;

        if (badTOTSlopeHigh || hasDefaultParameters) {
          c.badFitMask |= BadTOTHighFit;
          nBadTOTHigh++;
          logBadTOTHigh << " (" << asic << ", " << ch << ")";
        }
      }
    }
  }

  if (nMissingLow > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nMissingLow << " channels with no low charge injection calibration, "
                 << "replaced by the ASIC-half mean where available";
    LOG(debug) << "focal::HgcrocResponse: channels with no low charge injection calibration (asic, ch):" << logMissingLow.str();
  }
  if (nMissingHigh > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nMissingHigh << " channels with no high charge injection calibration, "
                 << "replaced by the ASIC-half mean where available";
    LOG(debug) << "focal::HgcrocResponse: channels with no high charge injection calibration (asic, ch):" << logMissingHigh.str();
  }
  if (nBadADC > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nBadADC << " channels with a bad ADC fit, "
                 << "replaced by the ASIC-half mean where available";
    LOG(debug) << "focal::HgcrocResponse: channels with a bad ADC fit (asic, ch):" << logBadADC.str();
  }
  if (nBadTOTLow > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nBadTOTLow << " channels with a bad low charge injection ToT fit, "
                 << "replaced by the ASIC-half mean where available";
    LOG(debug) << "focal::HgcrocResponse: channels with a bad low charge injection ToT fit (asic, ch):" << logBadTOTLow.str();
  }
  if (nBadTOTHigh > 0) {
    LOG(warning) << "focal::HgcrocResponse: " << nBadTOTHigh << " channels with a bad high charge injection ToT fit, "
                 << "replaced by the ASIC-half mean where available";
    LOG(debug) << "focal::HgcrocResponse: channels with a bad high charge injection ToT fit (asic, ch):" << logBadTOTHigh.str();
  }
}

void HgcrocResponse::flagDeadTOTHalves()
{
  const int nChannelsHalf = mConfig.mNChannels / 2;
  const uint16_t badTOTBit = BadTOTLowFit | BadTOTHighFit;

  for (int asic = 0; asic < mConfig.mNAsics; ++asic) {
    for (int half = 0; half < 2; ++half) {
      const int firstChannel = half * nChannelsHalf;
      const int lastChannel = firstChannel + nChannelsHalf;

      int nBadTOT = 0;
      for (int ch = firstChannel; ch < lastChannel; ++ch) {
        if (mCalibParams[asic][ch].badFitMask & badTOTBit) {
          nBadTOT++;
        }
      }

      if (nBadTOT < mConfig.mTOTBadHalfFraction * nChannelsHalf) {
        continue; // not a bad half
      }

      // for dead ToT half, flag all channels for a bad ToT fit
      for (int ch = firstChannel; ch < lastChannel; ++ch) {
        mCalibParams[asic][ch].badFitMask |= badTOTBit;
      }

      LOG(debug) << "focal::HgcrocResponse: ASIC " << asic << " half " << half << " has " << nBadTOT << "/" << nChannelsHalf << " channels with a bad ToT fit, "
                 << "all its ToT fits are flagged as bad (dead ToT half)";
    }
  }
}

// replace the flagged groups of bad fit parameters by the mean of the good channels:
// same ASIC half -> other half of the same ASIC -> whole detector
void HgcrocResponse::regularizeBadCalibParams()
{
  using Param = double ChannelCalib::*;

  const std::vector<Param> adcParams{&ChannelCalib::adcLowerFitRange, &ChannelCalib::adcHigherFitRange, &ChannelCalib::adcConst, &ChannelCalib::adcSlope};
  const std::vector<Param> totLowParams{&ChannelCalib::totLowerFitRange, &ChannelCalib::totHigherFitRange, &ChannelCalib::totConst, &ChannelCalib::totSlope};
  const std::vector<Param> totHighParams{&ChannelCalib::totConstHigh, &ChannelCalib::totSlopeHigh};

  const int nChannelsHalf = mConfig.mNChannels / 2;

  // fills 'mean' vector with mean of fit parameters over the good channels in [firstAsic, lastAsic) x [firstChannel, lastChannel); returns nGood
  auto meanOverGood = [this](uint16_t bit, const std::vector<Param>& fitParams, int firstAsic, int lastAsic,
                             int firstChannel, int lastChannel, std::vector<double>& mean) {
    mean.assign(fitParams.size(), 0.);
    int nGood = 0;

    for (int asic = firstAsic; asic < lastAsic; ++asic) {
      for (int ch = firstChannel; ch < lastChannel; ++ch) {
        const auto& c = mCalibParams[asic][ch];

        if (c.badFitMask & bit) {
          continue;
        }
        for (size_t i = 0; i < fitParams.size(); ++i) {
          mean[i] += c.*fitParams[i];
        }

        nGood++;
      }
    }

    if (nGood > 0) {
      for (double& m : mean) {
        m /= nGood;
      }
    }

    return nGood;
  };

  auto regularizeGroup = [&](int asic, int half, uint16_t bit, const std::vector<Param>& fitParams, const std::string& fitGroup) {
    const int firstChannel = half * nChannelsHalf;
    const int lastChannel = firstChannel + nChannelsHalf;

    std::vector<int> badFitChannels;

    for (int ch = firstChannel; ch < lastChannel; ++ch) {
      if (mCalibParams[asic][ch].badFitMask & bit) {
        badFitChannels.push_back(ch);
      }
    }

    if (badFitChannels.empty()) {
      return; // nothing to replace
    }

    std::vector<double> meanParams;

    const int otherHalf = 1 - half;
    // mean of the same half
    int nGood = meanOverGood(bit, fitParams, asic, asic + 1, firstChannel, lastChannel, meanParams);
    if (nGood == 0) {
      // no good channels in half, so take the mean of the other half
      nGood = meanOverGood(bit, fitParams, asic, asic + 1, otherHalf * nChannelsHalf, (otherHalf + 1) * nChannelsHalf, meanParams);
      if (nGood > 0) {
        LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " half " << half << " has no good channels for the " << fitGroup
                     << ", so its bad fits are replaced by the mean of the other half (half " << otherHalf << ") of the ASIC";
      }
    }
    if (nGood == 0) {
      // no good channels in the other half, so take the mean of the whole detector
      nGood = meanOverGood(bit, fitParams, 0, mConfig.mNAsics, 0, mConfig.mNChannels, meanParams);
      if (nGood > 0) {
        LOG(warning) << "focal::HgcrocResponse: ASIC " << asic << " has no good channels for the " << fitGroup
                     << ", so the bad fits of half " << half << " are replaced by the mean of the whole detector instead of the ASIC-half mean";
      }
    }
    if (nGood == 0) {
      // no good channels in the detector
      LOG(fatal) << "focal::HgcrocResponse: no good channels in the whole detector for the " << fitGroup;
    }

    // replacing bad parameters for the mean values
    for (int ch : badFitChannels) {
      auto& c = mCalibParams[asic][ch];
      for (size_t i = 0; i < fitParams.size(); ++i) {
        c.*fitParams[i] = meanParams[i];
      }
      if (bit == BadADCFit) { // the mean does not guarantee adcMargin > 0, so check it again after the substitution
        const double adcMargin = mConfig.mADCMax - c.adcConst - c.adcSlope * c.adcHigherFitRange;
        if (adcMargin <= 0) {
          LOG(fatal) << "focal::HgcrocResponse: the ADC fit parameters of ASIC " << asic << " half " << half
                     << " leave no ADC counts above the linear range (adcMargin = "
                     << adcMargin << " <= 0) after substitution of invalid fit parameters by the mean";
        }
      }
    }
  };

  for (int asic = 0; asic < mConfig.mNAsics; ++asic) {
    for (int half = 0; half < 2; ++half) {
      regularizeGroup(asic, half, BadADCFit, adcParams, "ADC fit");
      regularizeGroup(asic, half, BadTOTLowFit, totLowParams, "ToT low charge injection fit");
      regularizeGroup(asic, half, BadTOTHighFit, totHighParams, "ToT high charge injection fit");
    }
  }
}

void HgcrocResponse::computeDerivedParams()
{
  for (auto& row : mCalibParams) {
    for (auto& c : row) {
      c.gain = c.adcSlope;
      c.qLin = c.adcHigherFitRange;
      c.adcMargin = mConfig.mADCMax - c.adcConst - c.gain * c.qLin;
      c.qSat = c.qLin + 2 * c.adcMargin / c.gain;
      c.qTOT = c.totLowerFitRange;
    }
  }
}

} // namespace focal
} // namespace o2
