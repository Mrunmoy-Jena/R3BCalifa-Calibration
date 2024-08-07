/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BCALIFAMAPPED2CRYSTALCALPAR_H
#define R3BCALIFAMAPPED2CRYSTALCALPAR_H

#include "FairTask.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"

class TClonesArray;
class R3BCalifaMappingPar;
class R3BCalifaCrystalCalPar;
class R3BEventHeader;

class R3BCalifaMapped2CrystalCalPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BCalifaMapped2CrystalCalPar();

    /** Standard constructor **/
    R3BCalifaMapped2CrystalCalPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BCalifaMapped2CrystalCalPar();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Search peaks and calibrate **/
    virtual void SearchPeaks();

    /** Virtual method Fit peaks **/
    virtual void FitPeaks();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Accessor functions **/
    const Int_t GetNumCrystals() { return fNumCrystals; }
    const Int_t GetCalRange_left() { return fMapHistos_left; }
    const Int_t GetCalRange_right() { return fMapHistos_right; }
    const Int_t GetCalRange_bins() { return fMapHistos_bins; }
    const Int_t GetNumPeaks() { return fNumPeaks; }
    const Double_t GetSigma() { return fSigma; }
    const Double_t GetThreshold() { return fThreshold; }
    const Int_t GetNumParameterFit() { return fNumParam; }
    const Int_t GetMinStadistics() { return fMinStadistics; }
    const TString GetSourceName() { return fSourceName; }

    TArrayF* GetEnergyPeaks() { return fEnergyPeaks; }

    void SetNumCrystals(Int_t numberCry) { fNumCrystals = numberCry; }
    void SetCalRange_left(Int_t Histos_left) { fMapHistos_left = Histos_left; }
    void SetCalRange_right(Int_t Histos_right) { fMapHistos_right = Histos_right; }
    void SetCalRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    void SetCalRangeP_left(Int_t Histos_left) { fMapHistos_leftp = Histos_left; }
    void SetCalRangeP_right(Int_t Histos_right) { fMapHistos_rightp = Histos_right; }
    void SetCalRangeP_bins(Int_t Histos_bins) { fMapHistos_binsp = Histos_bins; }
    void SetNumPeaks(Int_t numberpeaks) { fNumPeaks = numberpeaks; }
    void SetMaxPeaks(Int_t maxPeaks) { fMaxPeaks = maxPeaks; } 
    void SetMinSlope(Double_t minSlope) { fMinSlope = minSlope; }
    void SetMaxSlope(Double_t maxSlope) { fMaxSlope = maxSlope; }
    void SetMinSlopeP(Double_t minSlopeP) { fMinSlopeP = minSlopeP; }
    void SetMaxSlopeP(Double_t maxSlopeP) { fMaxSlopeP = maxSlopeP;}
    void SetSigma(Double_t sigma) { fSigma = sigma; }
    void SetThreshold(Double_t threshold) { fThreshold = threshold; }
    void SetChi2Threshold(Double_t chi2Threshold) { fChi2Threshold = chi2Threshold; }
    void SetSigLowThreshold(Double_t sigLowThreshold) { fSigLowThreshold = sigLowThreshold; }
    void SetSigHighThreshold(Double_t sigHighThreshold) { fSigHighThreshold = sigHighThreshold; }
    void SetMinWidth(Double_t minWidth) { fMinWidth = minWidth; }
    void SetGausRange(Double_t gausRange) { fGausRange = gausRange; }
    void SetGausRangeP(Double_t gausRangeP) { fGausRangeP = gausRangeP; }
    void SetGausBaseEnergy(Double_t gausBaseEnergy) { fGausBaseEnergy = gausBaseEnergy; }
    void SetNumParameterFit(Int_t numberParFit) { fNumParam = numberParFit; }
    void SetMinStadistics(Int_t minstad) { fMinStadistics = minstad; }
    void SetSourceName(TString sourceName) { fSourceName = sourceName; }

    void SetDebugMode(Bool_t debug) { fDebugMode = debug; }

    void SetEnergyPeaks(TArrayF* thePeaks)
    {
        fEnergyPeaks = thePeaks;
        fNumPeaks = thePeaks->GetSize();
    }

  private:
    void SetParameter(); 
    Double_t FindChisquare(Double_t* X, Double_t* Y, Double_t* eX, Int_t ndf, TF1* f);
    
    Bool_t fDebugMode;
    Int_t fNumCrystals;
    Int_t fMapHistos_left; // gamma range
    Int_t fMapHistos_right;
    Int_t fMapHistos_bins;
    Int_t fMapHistos_leftp; // particle range
    Int_t fMapHistos_rightp;
    Int_t fMapHistos_binsp;

    Int_t fNumParam;
    Int_t fMinStadistics;

    Int_t fNumPeaks;
    Int_t fMaxPeaks;
    Double_t fSigma;
    Double_t fThreshold;
    Double_t fChi2Threshold;
    Double_t fSigLowThreshold;
    Double_t fSigHighThreshold;
    Double_t fMinWidth;

    Double_t fMinSlope;
    Double_t fMinSlopeP;
    Double_t fMaxSlope;
    Double_t fMaxSlopeP;

    Double_t fGausRange;
    Double_t fGausRangeP;
    Double_t fGausBaseEnergy;

    TString fSourceName;

    TArrayF* fEnergyPeaks;
    Double_t* fChannelPeaks;

    R3BCalifaMappingPar* fMap_Par;     /**< Parameter container with mapping. >*/
    R3BCalifaCrystalCalPar* fCal_Par;  /**< Container for Cal parameters. >*/
    TClonesArray* fCalifaMappedDataCA; /**< Array with CALIFA Mapped-input data. >*/

    //FILE *outfile;
    TFile *outrootfile;
    TTree *outroottree;

    TH1F** fh_Map_energy_crystal;
    TH2F* fh2_Map_crystal_gamma;
    TH2F* fh2_Map_crystal_proton;
    TH2F* fh_peak_crystal_gamma;
    TH2F* fh_peak_crystal_proton;
    TH2F* fh_sigma_crystal_gamma;
    TH2F* fh_sigma_crystal_proton;
    TH2F** fh2_peak_cal;
    TH2F** fh2_residual_energy;
    TH2F* fh2_slope_crystalID;
    TH2F* fh2_intercept_crystalID;
    TH2F** fh2_sig_crystal;
    TH2F* fh2_chi2_crystal;
    TH1F* fh_numPeak;
    TH2F* fh2_resolution_crystalID;
    TH1F** fh_Map_fit;
    TF1* f1 = nullptr;

  public:
    ClassDef(R3BCalifaMapped2CrystalCalPar, 2);
};

#endif
