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


    /** Virtual method Pulser Calibration **/
    virtual void PulserCalibration();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Accessor functions **/
    const Int_t GetNumCrystals() { return fNumCrystals; }
    const Int_t GetCalRange_left() { return fMapHistos_left; }
    const Int_t GetCalRange_right() { return fMapHistos_right; }
    const Int_t GetCalRange_bins() { return fMapHistos_bins; }
    const Int_t GetNumPeaks() { return fNumPeaks; }
    const Int_t GetNumVoltages_gamma() { return fNumVoltages_gamma; }
    const Int_t GetNumVoltages_proton() { return fNumVoltages_proton; }
        
    const Double_t GetSigma() { return fSigma; }
    const Double_t GetThreshold() { return fThreshold; }
    const Int_t GetNumParameterFit() { return fNumParam; }
    const Int_t GetMinStadistics() { return fMinStadistics; }
    const Int_t GetRunId() { return fRunId; }
    
    const TString GetSourceName() { return fSourceName; }
    const TString GetSpectrumName() { return fSpectrumName; }	
    const TString GetSpectrumName_gamma() { return fSpectrumName_gamma; }	
    const TString GetSpectrumName_proton() { return fSpectrumName_proton; }	
    
    const TString GetoutputName() { return foutputName; }
    const TString GetCalName() { return fCalName; }
    const TString Getcalifamapfilename() { return fcalifamapfilename; }
    const TString Getanglesfilename() { return fanglesfilename; }
    const TString GetPeaksCalibrated() { return fPeaksCalibrated; }
    const TString GetPeakErrors() { return fPeakErrors; }
    
    const Double_t GetMaxSigma() { return fMaxSigma; }
    const Int_t GetMinPeakEvents() { return fMinPeakEvents; }
    const Int_t GetPulserNumber() { return fPulserNumber; }

    TArrayF* GetEnergyPeaks() { return fEnergyPeaks; }
    TArrayF* GetPulserVoltages_gamma() { return fPulserVoltages_gamma; }
    TArrayF* GetPulserVoltages_proton() { return fPulserVoltages_proton; }
    

    void SetNumCrystals(Int_t numberCry) { fNumCrystals = numberCry; }
    void SetCalRange_left(Int_t Histos_left) { fMapHistos_left = Histos_left; }
    void SetCalRange_right(Int_t Histos_right) { fMapHistos_right = Histos_right; }
    void SetCalRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    void SetCalRangeP_left(Int_t Histos_left) { fMapHistos_leftp = Histos_left; }
    void SetCalRangeP_right(Int_t Histos_right) { fMapHistos_rightp = Histos_right; }
    void SetCalRangeP_bins(Int_t Histos_bins) { fMapHistos_binsp = Histos_bins; }
    
    void SetSourcePeaks_left(Int_t SourcePeaks_left) { fSourcePeaks_left = SourcePeaks_left; }
    void SetSourcePeaks_right(Int_t SourcePeaks_right) { fSourcePeaks_right = SourcePeaks_right; }
    void SetPulserPeaks_left(Int_t PulserPeaks_left) { fPulserPeaks_left = PulserPeaks_left; }
    void SetPulserPeaks_right(Int_t PulserPeaks_right) { fPulserPeaks_right = PulserPeaks_right; }
    
    void SetSourcePeaksP_left(Int_t SourcePeaksP_left) { fSourcePeaksP_left = SourcePeaksP_left; }
    void SetSourcePeaksP_right(Int_t SourcePeaksP_right) { fSourcePeaksP_right = SourcePeaksP_right; }
    void SetPulserPeaksP_left(Int_t PulserPeaksP_left) { fPulserPeaksP_left = PulserPeaksP_left; }
    void SetPulserPeaksP_right(Int_t PulserPeaksP_right) { fPulserPeaksP_right = PulserPeaksP_right; }    

    void SetNumPeaks(Int_t numberpeaks) { fNumPeaks = numberpeaks; }
    void SetNumVoltages_gamma(Int_t numbervoltages_gamma) { fNumVoltages_gamma = numbervoltages_gamma; }
    void SetNumVoltages_proton(Int_t numbervoltages_proton) { fNumVoltages_proton = numbervoltages_proton; }
    
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
    void SetSourceName(TString SourceName) { fSourceName = SourceName; }
    void SetSpectrumName(TString SpectrumName) { fSpectrumName = SpectrumName; }
    void SetSpectrumName_gamma(TString SpectrumName_gamma) { fSpectrumName_gamma = SpectrumName_gamma; }
    void SetSpectrumName_proton(TString SpectrumName_proton) { fSpectrumName_proton = SpectrumName_proton; }
    
    void SetoutputName(TString outputName) { foutputName = outputName; }
    void SetCalName(TString CalName) { fCalName = CalName; }
    void SetRunId(Int_t RunId) { fRunId = RunId; }
    void Setcalifamapfilename(TString califamapfilename) { fcalifamapfilename = califamapfilename; }
    void Setanglesfilename(TString anglesfilename) { fanglesfilename = anglesfilename; }
    void SetPeaksCalibrated(TString PeaksCalibrated) { fPeaksCalibrated = PeaksCalibrated; }
    void SetPeakErrors(TString PeakErrors) { fPeakErrors = PeakErrors; }
    
    void SetMaxSigma(Double_t MaxSigma) { fMaxSigma = MaxSigma; }
    void SetMinPeakEvents(Int_t MinPeakEvents) { fMinPeakEvents = MinPeakEvents; }
    void SetPulserNumber(Int_t PulserNumber) { fPulserNumber = PulserNumber; }
    
    void SetOffsetCalibration(Double_t OffsetCalibration) { fOffsetCalibration = OffsetCalibration; }
    
    void SetDebugMode(Bool_t debug) { fDebugMode = debug; }

    void SetEnergyPeaks(TArrayF* thePeaks)
    {
        fEnergyPeaks = thePeaks;
        fNumPeaks = thePeaks->GetSize();
    }
    
    void SetPulserVoltages_gamma(TArrayF* thePulserVoltages_gamma)
    {
        fPulserVoltages_gamma = thePulserVoltages_gamma;
        fNumVoltages_gamma = thePulserVoltages_gamma->GetSize();
    }
    
    void SetPulserVoltages_proton(TArrayF* thePulserVoltages_proton)
    {
        fPulserVoltages_proton = thePulserVoltages_proton;
        fNumVoltages_proton = thePulserVoltages_proton->GetSize();
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
    
    Int_t fSourcePeaks_left;
    Int_t fSourcePeaks_right;
    Int_t fPulserPeaks_left;
    Int_t fPulserPeaks_right;
    Int_t fSourcePeaksP_left;
    Int_t fSourcePeaksP_right;
    Int_t fPulserPeaksP_left;
    Int_t fPulserPeaksP_right;

    Int_t fNumParam;
    Int_t fMinStadistics;
    Int_t fNumParam_Source;
    
    Int_t fNumPeaks;
    Int_t fNumVoltages_gamma;
    Int_t fNumVoltages_proton;
    
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

    Double_t fMaxSigma;
    Int_t fMinPeakEvents;
    Int_t fPulserNumber;
    Int_t fRunId;
    
    Double_t fOffsetCalibration;
    
    TString fSourceName;
    TString fSpectrumName;
    TString fSpectrumName_gamma;
    TString fSpectrumName_proton;
    
    TString foutputName;
    TString fCalName;
    TString fcalifamapfilename;
    TString fanglesfilename;
    TString fPeaksCalibrated;
    TString fPeakErrors;
    
    TArrayF* fEnergyPeaks;
    TArrayF* fPulserVoltages_gamma;
    TArrayF* fPulserVoltages_proton;
    
    
    Double_t* fChannelPeaks;
    Double_t* fChannelPeaks_source;
    Double_t* fChannelPeaks_pulser;

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
