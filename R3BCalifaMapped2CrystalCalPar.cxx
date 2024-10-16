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

#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TVector3.h"
#include "TPolyMarker.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TKey.h"
#include "TLine.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairParAsciiFileIo.h"

#include "R3BCalifaCrystalCalPar.h"
#include "R3BCalifaMapped2CrystalCalPar.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaMappingPar.h"

#include <vector>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <TCanvas.h>
using namespace std;

R3BCalifaMapped2CrystalCalPar::R3BCalifaMapped2CrystalCalPar()
    : R3BCalifaMapped2CrystalCalPar("R3B CALIFA Calibration Parameters Finder ", 1)
{
}

R3BCalifaMapped2CrystalCalPar::R3BCalifaMapped2CrystalCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMap_Par(NULL)
    , fCal_Par(NULL)
    , fCalifaMappedDataCA(NULL)
    , fNumCrystals(1)
    , fNumParam(0)
    , fMinStadistics(100)
    , fMapHistos_left(0)
    , fMapHistos_right(0)
    , fMapHistos_bins(0)
    , fMapHistos_leftp(0)
    , fMapHistos_rightp(0)
    , fMapHistos_binsp(0)
    , fNumPeaks(0)
    , fMaxPeaks(20)
    , fSigma(0)
    , fThreshold(0)
    , fChi2Threshold(0.0)
    , fSigLowThreshold(5.0)
    , fSigHighThreshold(500.0)
    , fMinWidth(1.0)
    , fGausRange(6.0)
    , fGausRangeP(30.0)
    , fGausBaseEnergy(511.0)
    , fMinSlope(0.)
    , fMaxSlope(100.)
    , fMinSlopeP(0.)
    , fMaxSlopeP(100.)
    , fEnergyPeaks(NULL)
    , fDebugMode(0)
    , fMaxSigma(50.0)
    , fMinPeakEvents(100)
{
}

R3BCalifaMapped2CrystalCalPar::~R3BCalifaMapped2CrystalCalPar()
{
    LOG(info) << "R3BCalifaMapped2CrystalCalPar: Delete instance";
    if (fCalifaMappedDataCA)
        delete fCalifaMappedDataCA;
    if (fEnergyPeaks)
        delete fEnergyPeaks;
}

void R3BCalifaMapped2CrystalCalPar::SetParContainers()
{
    cout<<"SetParContainers() called"<<endl;
    // Parameter Container
    // Reading califaMappingPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }

    fMap_Par = dynamic_cast<R3BCalifaMappingPar*>(rtdb->getContainer("califaMappingPar"));
    if (!fMap_Par)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaMappingPar container";
    }
    else
    {
        LOG(info) << "R3BCalifaMapped2CrystalCalPar:: califaMappingPar container open";
    }
}

void R3BCalifaMapped2CrystalCalPar::SetParameter()
{
    cout<<"SetParameter() called"<<endl;
    if (!fMap_Par)
    {
        LOG(warn) << "R3BCalifaMapped2CrystalCalPar::Container califaMappingPar not found.";
    }
    //--- Parameter Container ---
    fNumCrystals = fMap_Par->GetNumCrystals(); // Number of crystals x 2
    
    LOG(info) << "R3BCalifaMapped2CrystalCalPar::NumCry " << fNumCrystals;
    // fMap_Par->printParams();
}

InitStatus R3BCalifaMapped2CrystalCalPar::Init()
{
    LOG(info) << "R3BCalifaMapped2CrystalCalPar::Init()";

    if (!fEnergyPeaks)
    {
        fEnergyPeaks = new TArrayF;
        fEnergyPeaks->Set(fNumPeaks);
    }

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() FairRootManager not found";
        return kFATAL;
    }

    fCalifaMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("CalifaMappedData"));
    if (!fCalifaMappedDataCA)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() CalifaMappedData not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() FairRuntimeDb not found";
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BCalifaCrystalCalPar*>(rtdb->getContainer("califaCrystalCalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaCrystalCalPar container";
        return kFATAL;
    }

    // Initiate output file
    outrootfile = TFile::Open(fSpectrumName,"UPDATE");
    if (!outrootfile) {
        cout << "no spectrum.root found, create new file" << endl;
    	outrootfile = new TFile(fSpectrumName,"RECREATE");
    outrootfile->cd();
    outroottree = new TTree("genT","General Tree");
    }
    cout << "spectrum.root file name: " << fSpectrumName << endl;

    // Set container with mapping parameters
    SetParameter();

    // Create histograms for crystal calibration
    char name1[100];
    char name2[100];
    char name3[100];
    char name4[100];
    char name5[100];
    char name6[100];
    Int_t fright, fleft, fbins;
    if (fSourceName == "spectrum") 
    {
        fh_Map_energy_crystal = new TH1F*[fNumCrystals];
    }
    for (Int_t i = 0; i < fNumCrystals; i++)
        //if (fMap_Par->GetInUse(i + 1) == 1)
        {

            sprintf(name1, "fh_Map_energy_crystal_"+fSourceName+"_%i", i + 1);
	        sprintf(name2, "fh2_peak_cal_%i", i + 1);
	        sprintf(name3, "fh_map_crystal_fit_22Na_%i", i + 1);
	        sprintf(name4, "fh_map_crystal_fit_60Co_%i", i + 1);
	        sprintf(name5, "fh_map_crystal_fit_AmBe_%i", i + 1);
	        sprintf(name6, "fh2_residual_energy_%i", i + 1);
	        
            if (i < fMap_Par->GetNumCrystals() / 2)
            {
                fright = fMapHistos_right;
                fleft = fMapHistos_left;
                fbins = fMapHistos_bins;
            }
            else
            {
                fright = fMapHistos_rightp;
                fleft = fMapHistos_leftp;
                fbins = fMapHistos_binsp;
            }
            
	        if (fSourceName == "spectrum")
	        {
                    fh_Map_energy_crystal[i] = new TH1F(name1, name1, fbins, fleft, fright);
	        }
        }

    

    if (fSourceName == "spectrum")
    {
        fh2_Map_crystal_gamma = new TH2F("fh2_Map_crystal_gamma","fh2_Map_crystal_gamma;crystal ID;Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh2_Map_crystal_proton = new TH2F("fh2_Map_crystal_proton","fh2_Map_crystal_proton;crystal ID;Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
    }

    return kSUCCESS;
}


InitStatus R3BCalifaMapped2CrystalCalPar::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BCalifaMapped2CrystalCalPar::Exec(Option_t* opt)
{
    if (fSourceName != "spectrum") return;

    Int_t nHits = fCalifaMappedDataCA->GetEntries();
    if (!nHits) return;
        
    R3BCalifaMappedData** MapHit = new R3BCalifaMappedData*[nHits];
    Int_t crystalId = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        MapHit[i] = dynamic_cast<R3BCalifaMappedData*>(fCalifaMappedDataCA->At(i));
        crystalId = MapHit[i]->GetCrystalId();
        
        // Fill histograms
        //if (fMap_Par->GetInUse(crystalId) == 1)
        //{
	        Double_t fleft, fright;
	        if (crystalId<=fNumCrystals/2)
	        {
		        fleft = fMapHistos_left;
		        fright = fMapHistos_right;
	        } 
	        else
	        {
		        fleft = fMapHistos_leftp;
		        fright = fMapHistos_rightp;
	        }

	        if (MapHit[i]->GetEnergy()>=fleft and MapHit[i]->GetEnergy()<=fright)
	        {
                    fh_Map_energy_crystal[crystalId - 1]->Fill(MapHit[i]->GetEnergy());
	            if (crystalId < fNumCrystals/2)
	            {
                    fh2_Map_crystal_gamma->Fill(crystalId-1,MapHit[i]->GetEnergy());
	            } 
	            else
	            {
		            fh2_Map_crystal_proton->Fill(crystalId-1,MapHit[i]->GetEnergy());
		        }
	        }
	    //}	
	
    }

    if (MapHit)
        delete[] MapHit;
    return;
}

void R3BCalifaMapped2CrystalCalPar::Reset() {}

void R3BCalifaMapped2CrystalCalPar::FinishEvent() {}

void R3BCalifaMapped2CrystalCalPar::FinishTask()
{
    cout << "finishTask() called" << endl;
    
    outrootfile->Write();
    outrootfile->Close();
}


//_______________________________________________________//

//________________Pulser_Calibration_____________________//

//_______________________________________________________//

void R3BCalifaMapped2CrystalCalPar::PulserCalibration()  
{
    cout << "Pulser_Calibration() called" <<endl;
  
    //Initialization
    SetParContainers();
    SetParameter();
  
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() FairRuntimeDb not found";
        return;
    }
  
    //_____________________________________open .par file and add to FairRuntimeDb__________________________________________________//
    FairParAsciiFileIo* parIo = new FairParAsciiFileIo();  
    parIo->open(fcalifamapfilename, "in");
    rtdb->setFirstInput(parIo);
    rtdb->addRun(fRunId);
	rtdb->initContainers(fRunId);
	rtdb->print();

    fCal_Par = dynamic_cast<R3BCalifaCrystalCalPar*>(rtdb->getContainer("califaCrystalCalPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BCalifaMapped2CrystalCalPar::Init() Couldn't get handle on califaCrystalCalPar container";
        return;
    }
    
    //_________________________________________________open_files___________________________________________________________________//
    //open spectrum.root file
    TFile* spectrum_file_gamma = TFile::Open(fSpectrumName_gamma, "READ");
    if (!spectrum_file_gamma) 
    {
    
        cerr << "No gamma_spectrum.root file found. Please perform histogram creation first!" << endl;
        return;
    }
    else
    {
        cout << "use gamma_spectrum file: " << fSpectrumName_gamma << endl;
    }
    
    TFile* spectrum_file_proton = TFile::Open(fSpectrumName_proton, "READ");
    if (!spectrum_file_proton) 
    {
    
        cerr << "No proton_spectrum.root file found. Please perform histogram creation first!" << endl;
        return;
    }
    else
    {
        cout << "use proton_spectrum file: " << fSpectrumName_proton << endl;
    }
  
    //open calibrated.root file
    TFile* local_outputFile = TFile::Open(foutputName, "RECREATE");
  
    // open file for calibrated pulser peaks in gamma range
    vector<int> crystalnumber_gamma;
    vector<vector<double>> pulserpeaks_gamma;
    
    
    ofstream crystals_errors;    
    crystals_errors.open(fPeakErrors, ios::out);
    if (!crystals_errors.is_open()) 
    {
        cerr << "Error: The file '" << fPeakErrors << "' could not be created!" << endl;
    }
    

    //___________________________________________________calibration____________________________________________________________________//
    
    //fNumPeaks corresponds to the number of source energy values ​​specified in the macro, fNumVoltages corresponds to the number of different pulser signals
    Double_t Num_all_peaks_gamma = fNumVoltages_gamma+fNumPeaks;
    Double_t Num_all_peaks_proton = fNumVoltages_proton+fNumPeaks; 
    
    // create histograms
    TH1F* local_fh_Map_energy_crystal[fNumCrystals];
    
    TH2F* PeaksFoundVsCrystal = new TH2F("Peaks_vs_Crystal", "Peaks vs Crystal ID; Crystal ID; Bin number", fNumCrystals, 0, fNumCrystals, 10000, 0, 30000);
    TH2F* SourcePeaksFoundGammaVsCrystal = new TH2F("gamma_source_peaks_vs_Crystal", "SourcePeaksGamma vs Crystal ID; Crystal ID; Bin number", fNumCrystals, 0, fNumCrystals, 2000, 0, 2000);
    TH2F* SourcePeaksFoundProtonVsCrystal = new TH2F("proton_source_peaks_vs_Crystal", "SourcePeaksProton vs Crystal ID; Crystal ID; Bin number", fNumCrystals, 0, fNumCrystals, 2000, 0, 200);
    TH2F* PeaksCalibratedGammaVsCrystal = new TH2F("gamma_peaks_calibrated_vs_Crystal", "PeaksCalibratedGamma vs Crystal ID; Crystal ID; Energy in keV", fNumCrystals, 0, fNumCrystals, 5000, 0, 30000);
    TH2F* PeaksCalibratedProtonVsCrystal = new TH2F("proton_peaks_calibrated_vs_Crystal", "PeaksCalibratedProton vs Crystal ID; Crystal ID; Energy in keV", fNumCrystals, 0, fNumCrystals, 5000, 0, 400000);
    TH2F* RangeFactorVsCrystal = new TH2F("Range_Factor_vs_Crystal", "RangeFactor vs Crystal ID; Crystal ID; RangeFactor", fNumCrystals, 0, fNumCrystals, 3000, 0, 50);
    TH2F* CalibratedEnergies_VoltageVsCrystal = new TH2F("proton_E/V_vs_Crystal", "E/V vs Crystal ID; Crystal ID; E/V", fNumCrystals, 0, fNumCrystals, 2000, 0, 500);
    
    vector<TH2F*> CalibratedEnergies_differenceVsCrystal(fNumVoltages_proton-fNumVoltages_gamma);
    
    for (int j = 0; j < (fNumVoltages_proton-fNumVoltages_gamma); ++j) {
        string histName = "proton_pulser_Peak_" + to_string(j + 1 + fNumVoltages_gamma) + "_E_fit/E_V_vs_CrystalID";
        string histTitle = "proton pulser Peak " + to_string(j + 1 + fNumVoltages_gamma) + " E_fit/E_V vs Crystal ID; Crystal ID; percentage deviation";

        CalibratedEnergies_differenceVsCrystal[j] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 3000, -1, 5);
    }
    
    // Arrays to store the x and y values ​​for the graphs
    vector<double> xValues, xValues_gamma, xValues_proton, BinRangePeak_gamma, yOffsetValues, ySlopeValues, yOffsetPulserValues, ySlopePulserValues, ychi2, ypvalue;

    vector<vector<double>> gamma_ySigmaValues(Num_all_peaks_gamma, vector<double>());
    vector<vector<double>> proton_ySigmaValues(Num_all_peaks_proton, vector<double>());
    
    // 2D Array to store calibrated Pulser Energies from Gamma Range
    vector<vector<double>> PulserEnergiesCalibrated(fNumVoltages_gamma+1, vector<double>());
    

    Int_t numPars = 2; // Number of parameters=2 by default
    Int_t nfound = 0;
    Int_t nfound_source = 0;
    Int_t nfound_pulser = 0;
    
    if (fNumParam)  // if num of parameters is explicitly specified then get this value
    {
        numPars = fNumParam;
    }
    
    //get  pulser voltages
    Double_t V_gamma[fNumVoltages_gamma];
    Double_t V_proton[fNumVoltages_proton];
    
    
    for (Int_t p=0; p<fNumVoltages_gamma; p++)
    {
        V_gamma[p] = fPulserVoltages_gamma->GetAt(p);
    }
    
    for (Int_t p=0; p<fNumVoltages_proton; p++)
    {
        V_proton[p] = fPulserVoltages_proton->GetAt(p);
    }
    
    fCal_Par->SetNumCrystals(fNumCrystals);
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(numPars * fNumCrystals);
    Int_t fright, fleft;

    for (Int_t i=0; i<fNumCrystals; i++)
    {    

        Int_t Num_all_peaks = (i < fNumCrystals / 2) ? Num_all_peaks_gamma : Num_all_peaks_proton;

        TSpectrum* ss = new TSpectrum(fNumPeaks);
        TSpectrum* ss_pulser = new TSpectrum(Num_all_peaks-fNumPeaks);
        
        char histName[100];
        sprintf(histName, "fh_Map_energy_crystal_spectrum_%i", i + 1);

        if (i<fNumCrystals/2)
            local_fh_Map_energy_crystal[i] = (TH1F*)spectrum_file_gamma->Get(histName);
        else
            local_fh_Map_energy_crystal[i] = (TH1F*)spectrum_file_proton->Get(histName);
          
          
        if (local_fh_Map_energy_crystal[i] && local_fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics) 
        {
        
            TH1F* subHist_source = (TH1F*)local_fh_Map_energy_crystal[i]->Clone();
            TH1F* subHist_pulser = (TH1F*)local_fh_Map_energy_crystal[i]->Clone();
                
            if (i<fNumCrystals/2)
            {
                subHist_source->GetXaxis()->SetRangeUser(fSourcePeaks_left, fSourcePeaks_right);
                subHist_pulser->GetXaxis()->SetRangeUser(fPulserPeaks_left, fPulserPeaks_right);
            }
            else
            {
                subHist_source->GetXaxis()->SetRangeUser(fSourcePeaksP_left, fSourcePeaksP_right);
                subHist_pulser->GetXaxis()->SetRangeUser(fPulserPeaksP_left, fPulserPeaksP_right);
            }
            
            gErrorIgnoreLevel = kError;     // suppress Warning in <TSpectrum::SearchHighRes>: Peak buffer full
            
            nfound_source = ss->Search(subHist_source, fSigma, "", fThreshold);
            nfound_pulser = ss_pulser->Search(subHist_pulser, fSigma, "", fThreshold);
            
            gErrorIgnoreLevel = kInfo;
            
            nfound = nfound_source + nfound_pulser;

            delete subHist_source;
            delete subHist_pulser;
        
            //safe histograms in outputFile    
            local_outputFile->cd();
            
            fChannelPeaks_source = (Double_t*) ss->GetPositionX();

            fChannelPeaks_pulser = (Double_t*) ss_pulser->GetPositionX();

            fChannelPeaks = new Double_t[nfound_source + nfound_pulser];

            for (int j = 0; j < nfound_source; j++) 
            {
                fChannelPeaks[j] = fChannelPeaks_source[j];
            }

            for (int j = 0; j < nfound_pulser; j++) 
            {
                fChannelPeaks[nfound_source + j] = fChannelPeaks_pulser[j];
            }


            Int_t idx[nfound];
            TMath::Sort(nfound, fChannelPeaks, idx, kFALSE);	//kFALSE: sort in ascending order (from lowest to highest) -> sorted indexes: idx (kTRUE: in decscending)

            // Calibrated Spectrum
            Double_t X[nfound + 1];	//array with size nfound+1 elements
            Double_t Y[nfound + 1];

            cout << "Crystal " << i + 1 << endl;  
                
            if (nfound == Num_all_peaks)	//Check that the number of source peaks found matches the expected number
            {          
                xValues.push_back(i + 1);
            
                if (i < fNumCrystals / 2)
                {
                    xValues_gamma.push_back(i + 1);
                }
                else
                {
                    xValues_proton.push_back(i + 1);
                }
            
                crystals_errors << i + 1;
                
        	    //gauss fit each found peak and check sigma and events
                for (Int_t j=0; j<Num_all_peaks; j++)
                {
                    Double_t posX = fChannelPeaks[idx[j]];
                    
                    TF1 * gaussfit;
                    
                    if (i<fNumCrystals/2) 
                    {
		                gaussfit = new TF1(Form("gaussfit_%d_%d", i, j),"gaus",posX-fSigma*15,posX+fSigma*15);
		            }
                    else 
                    {
                    	gaussfit = new TF1(Form("gaussfit_%d_%d", i, j),"gaus",posX-fSigma*1.5,posX+fSigma*1.5);
                    }

	                local_fh_Map_energy_crystal[i]->Fit(gaussfit, "RQ+");
	                
	                Double_t peakHeight = gaussfit->GetParameter(0);
	                Double_t mean = gaussfit->GetParameter(1);
	                Double_t sigma = gaussfit->GetParameter(2);

	                if (i < fNumCrystals / 2)
	                {
	                    gamma_ySigmaValues[j].push_back(sigma);
                    }
                    else
                    {
                        proton_ySigmaValues[j].push_back(sigma);
                    }
                    
                    X[j] = mean;    
                        
                    PeaksFoundVsCrystal->Fill(i + 1, X[j]); 
                    
                    if (j<fNumPeaks) 
                    {
                    	Y[j] = fEnergyPeaks->GetAt(j);	//fills the Y array with reference energy values from fEnergyPeaks
                    	
                    	cout << "Source Peak " << j + 1 << ", bin number: " << X[j] << endl;
                    	
                        if (i<fNumCrystals/2) 
                        {
		                    SourcePeaksFoundGammaVsCrystal->Fill(i + 1, X[j]); 
		                }
                        else 
                        {
                        	SourcePeaksFoundProtonVsCrystal->Fill(i + 1, X[j]); 
                        }
                    	
	                }
	                else
	                {
	                    cout << "Pulser Peak " << j + 1 - fNumPeaks << ", bin number: " << X[j] << endl;
	                }
	                
	                
                    if (j==(fNumPeaks-1) && i<fNumCrystals/2) 
                    {
                        BinRangePeak_gamma.push_back(X[j]);
                    }
	                
                    Double_t peakArea = peakHeight * sigma * sqrt(2 * M_PI);

	                if (sigma > fMaxSigma) // Check if the standard deviation is too large
	                {
                        LOG(error)  << "Peak " << j + 1 << ": Sigma too large ( " << sigma << " )!";
                            
                        crystals_errors << ", " << "Peak " << j + 1 << ": Sigma";
                    }
                    
                    if (peakArea < fMinPeakEvents) // Check if the number of events at the peak is too small
                    {
                        LOG(error)  << "Peak " << j + 1 << ": Too few events at peak ( " << peakArea << " )!";
                        
                        crystals_errors << ", " << "Peak " << j + 1 << ": events";
                    }
                }
            
                local_fh_Map_energy_crystal[i]->Write();
            
                X[nfound]=0.;	//set last value in array equal to 0 (sentinel value)
                Y[nfound]=0.;
                
                if (i < fMap_Par->GetNumCrystals() / 2)
                {
                    fright = fMapHistos_right;
                    fleft = fMapHistos_left;
                }
                else
                {
                    fright = fMapHistos_rightp;
                    fleft = fMapHistos_leftp;
                }
                
                Double_t CalParams[fNumParam];  //standard, [0]: offset, [1]: slope
                
                if (i<fNumCrystals/2)
                { 
                    TF1* f1fit = nullptr;

                    if (fOffsetCalibration != 0) 
                    {
                        f1fit = new TF1("f1fit", Form("%f + [0]*x", fOffsetCalibration), fleft, fright);
                    }
                    else
                    {
                        f1fit = new TF1("f1fit", "[0]+[1]*x", fleft, fright);
                    }

                    TGraph* graph = new TGraph(fNumPeaks, X, Y);
                    graph->Fit("f1fit", "Q"); // Quiet mode (minimum printing)
                    
                    cout << "Fit X values: ";
                    for (int a = 0; a < fNumPeaks; a++) {
                        cout << X[a];
                        if (a < fNumPeaks - a) 
                        {
                            cout << ", ";
                        }
                    }
                    cout << endl;
                    
                    cout << "Fit Y values: ";
                    for (int b = 0; b < fNumPeaks; b++) {
                        cout << Y[b];
                        if (b < fNumPeaks - b) 
                        {
                            cout << ", ";
                        }
                    }
                    cout << endl;
                    
                    //get offset and slope
                    if (fOffsetCalibration != 0) 
                    {
                        CalParams[0] = fOffsetCalibration;
                        
                        for (Int_t h = 1; h < fNumParam; h++)
                        {
                            CalParams[h] = f1fit->GetParameter(h-1);
                        }
                    }
                    else
                    {
                        for (Int_t h = 0; h < fNumParam; h++)
                        {
                            CalParams[h] = f1fit->GetParameter(h);
                        }
                    }
                    
                    cout << "Offset: " << CalParams[0] << endl;                
                    cout << "Slope: " << CalParams[1] << endl;

                    // Fill the histograms
                    for (Int_t c = 0; c < fNumPeaks; c++)
                    {
                        PeaksCalibratedGammaVsCrystal->Fill(i + 1, CalParams[0] + X[c] * CalParams[1]);
                    }
                    
                    yOffsetValues.push_back(CalParams[0]);
                    ySlopeValues.push_back(CalParams[1]);
                }
                
                //_____________________________pulser_calibration_____________________________________//

                if (fNumVoltages_gamma != 0 && fNumVoltages_proton != 0)
                {

                    //_________________pulser_gamma_________________________//
                    
                    if (i<fNumCrystals/2)
                    {  
                    
                        PulserEnergiesCalibrated[0].push_back(i+1);
                        
                        cout << "Index: " << PulserEnergiesCalibrated[0].size() - 1 << endl;

                        Double_t PulserEnergycalibrated[fNumVoltages_gamma];
                        
                        Double_t PulserEnergyexpected[fNumVoltages_gamma];
                        Double_t pulser_sigma = fSigma/100;
                        
                        
                        for (Int_t v = 0; v < fNumVoltages_gamma; v++)
                        {
                            PulserEnergycalibrated[v] = CalParams[0] + X[fNumPeaks + v] * CalParams[1];
                            
                            PeaksCalibratedGammaVsCrystal->Fill(i + 1, PulserEnergycalibrated[v]);
                            
                            PulserEnergiesCalibrated[v+1].push_back(PulserEnergycalibrated[v]);
                            
                            cout << v + 1 << " Pulser energy calibrated: " << PulserEnergycalibrated[v] << endl;
                        }
                        
                        
                        PulserEnergyexpected[0] = PulserEnergycalibrated[0];
                        
                        for (Int_t e = 1; e < fNumVoltages_gamma; e++)
                        {
                            PulserEnergyexpected[e] = PulserEnergyexpected[0] * V_gamma[e]/V_gamma[0];
                          
                            if( PulserEnergyexpected[e] < (PulserEnergycalibrated[e] * (1 - pulser_sigma)) || PulserEnergyexpected[e] > (PulserEnergycalibrated[e] * (1 + pulser_sigma)) )
                            {
                                cout << "Pulser energy: " << PulserEnergycalibrated[e] << " does not correspond to expected energy by voltage_gamma: " << PulserEnergyexpected[e] << endl;
                            }
                            else
                            {
                                cout << e + 1 << " Pulser energy by voltage_gamma: " << PulserEnergyexpected[e] << endl;
                            }
                        }
                        
                        
                        //pass slope and offset from simple source calibration
                        for (Int_t p = 0; p < fNumParam; p++)
                        {
                            fCal_Par->SetCryCalParams(CalParams[p], fNumParam*i+p);
                        }   

                    }
                    
                    //_________________pulser_proton_________________________//
                    
                    else
                    {
                        
                        Double_t PulserEnergycalibrated[fNumVoltages_proton];
                    
                        // pulser bin numbers              
                        Double_t X_Pulser[fNumVoltages_proton];
                        
                        for (Int_t r = 0; r < fNumVoltages_proton; r++)
                        {
                            X_Pulser[r] = X[fNumPeaks + r];
                        }
                        
                        Double_t PulserParams[fNumParam];
                    
                        Double_t range_factor;
                        
                        // find corresponding Calibrated Pulser Voltages
                        bool found = false;

                        for (int col = 0; col < PulserEnergiesCalibrated[0].size(); col++) 
                        {
                            if (PulserEnergiesCalibrated[0][col] == (i + 1 - fNumCrystals/2)) 
                            {
                                found = true;
                                
                                cout << "Index " << i + 1 << " found at column: " << col << endl;
                                
                                for (int num = 0; num < fNumVoltages_gamma; num++)
                                {
                                    PulserEnergycalibrated[num] = PulserEnergiesCalibrated[num+1][col];
                                    
                                    PeaksCalibratedProtonVsCrystal->Fill(i + 1, PulserEnergycalibrated[num]);
                                    
                                    cout << num + 1 << " Pulser energy calibrated: " << PulserEnergycalibrated[num] << endl; 
                                    
                                    CalibratedEnergies_VoltageVsCrystal->Fill(i + 1, PulserEnergycalibrated[num]/V_proton[num]);
                                }
                                
                                range_factor = BinRangePeak_gamma[col]/X[fNumPeaks-1];
                                
                                RangeFactorVsCrystal->Fill(i + 1, range_factor);
                                
                                cout << "Range factor: " << range_factor << " (" << BinRangePeak_gamma[col] << ", " << X[fNumPeaks-1] << ")" << endl;
                                break;
                            }
                        }

                        if (!found) 
                        {
                            cout << "Warning: no corresponding calibrated pulser voltages from gamma range found!" << endl;
                        }
                        
                        TF1* f2fit = new TF1("f2fit", "[0]+[1]*x", fleft, fright);
                        
                        TGraphErrors* graph2 = new TGraphErrors(fNumVoltages_gamma);

                        for (int g = 0; g < fNumVoltages_gamma; g++) 
                        {
                            Double_t x_value = X_Pulser[g];
                            Double_t y_value = PulserEnergycalibrated[g]; 
                            Double_t error_y = 0.02 * y_value;
                            
                            graph2->SetPoint(g, x_value, y_value);
                            graph2->SetPointError(g, 0, error_y);
                        }
                        
                        graph2->Fit("f2fit", "Q");

                        Double_t chi2 = f2fit->GetChisquare();
                        Int_t ndf = f2fit->GetNDF();
                        Double_t pvalue = f2fit->GetProb();

                        for (Int_t f = 0; f < fNumParam; f++)
                        {
                            PulserParams[f] = f2fit->GetParameter(f);
                        }
                        
                        cout << "Pulser offset: " << PulserParams[0] << endl;
                        cout << "Pulser slope: " << PulserParams[1] << endl;
                        cout << "Chi²: " << chi2 << endl;
                        cout << "NDF: " << ndf << endl;
                        cout << "p-Wert: " << pvalue << endl;
                        
                        yOffsetPulserValues.push_back(PulserParams[0]);
                        ySlopePulserValues.push_back(PulserParams[1]);
                        ychi2.push_back(chi2);
                        ypvalue.push_back(pvalue);
                        
                        //pass slope and offset from pulser calibration
                        for (Int_t p = 0; p < fNumParam; p++)
                        {
                            fCal_Par->SetCryCalParams(PulserParams[p], fNumParam*i+p);
                        }  
                        
                        for (Int_t e = fNumVoltages_gamma; e < fNumVoltages_proton; e++)
                        {
                        
                            PulserEnergycalibrated[e] = PulserParams[0] + X_Pulser[e] * PulserParams[1];

                            Double_t PulserEnergycalibratedbyVoltage = PulserEnergycalibrated[0] * V_proton[e]/V_proton[0];
                            Double_t Pulsercalibrateddifference = (PulserEnergycalibrated[e]/PulserEnergycalibratedbyVoltage-1) * 100;
                            
                            PeaksCalibratedProtonVsCrystal->Fill(i + 1, PulserEnergycalibrated[e]);
                            CalibratedEnergies_VoltageVsCrystal->Fill(i + 1, PulserEnergycalibrated[e]/V_proton[e]);
                            CalibratedEnergies_differenceVsCrystal[e-fNumVoltages_gamma]->Fill(i + 1, Pulsercalibrateddifference);
                            
                            cout << e + 1 << " Pulser energy calibrated: " << PulserEnergycalibrated[e] << endl;
                            cout << e + 1 << " Pulser energy calibrated by voltage: " << PulserEnergycalibratedbyVoltage << endl;
                            cout << e + 1 << " Pulser energy calibrated percentage difference: " << Pulsercalibrateddifference << endl;
                        }
                    }
                }
                
                //___________________________________________________________________________________//

                crystals_errors << endl;
            
            }
            else
            {
            	cout << "Number of peaks found does not correspond to the expected number!" << endl;
            }
            
            delete ss;
            
            cout << endl;    
	    }    
    }
    
    crystals_errors << endl << endl;
    
    //______________histograms_______________//
   
    //create histograms
    TH2F* OffsetGammaVsCrystal = new TH2F("gamma_offset_vs_CrystalID", "OffsetGamma vs Crystal ID;Crystal ID;Offset", fNumCrystals, 0, fNumCrystals, 1200, -100, 100);
    TH2F* SlopeGammaVsCrystal = new TH2F("gamma_slope_vs_CrystalID", "SlopeGamma vs Crystal ID;Crystal ID;Slope", fNumCrystals, 0, fNumCrystals, 1000, 0, 10);
    TH2F* OffsetProtonVsCrystal = new TH2F("proton_offset_vs_CrystalID", "OffsetProton vs Crystal ID; Crystal ID; Offset", fNumCrystals, 0, fNumCrystals, 2000, -150, 50);
    TH2F* SlopeProtonVsCrystal = new TH2F("proton_slope_vs_CrystalID", "SlopeProton vs Crystal ID; Crystal ID; Slope", fNumCrystals, 0, fNumCrystals, 1000, 0, 60);
    TH2F* Chi2ProtonVsCrystal = new TH2F("proton_chi2_vs_CrystalID", "Chi2Proton vs Crystal ID; Crystal ID; Chi2", fNumCrystals, 0, fNumCrystals, 1000, 0, 1);
    TH2F* PValueProtonVsCrystal = new TH2F("proton_pvalue_vs_CrystalID", "PValueProton vs Crystal ID; Crystal ID; PValue", fNumCrystals, 0, fNumCrystals, 1100, 0, 1.1);
    
    vector<TH2F*> gamma_SigmaVsCrystal(Num_all_peaks_gamma);
    vector<TH2F*> proton_SigmaVsCrystal(Num_all_peaks_proton);
    
    for (int j = 0; j < fNumPeaks; ++j) {
        string histName = "gamma_source_Peak_" + to_string(j + 1) + "_Sigma_vs_CrystalID";
        string histTitle = "gamma source Peak " + to_string(j + 1) + " Sigma vs Crystal ID; Crystal ID; Sigma";

        gamma_SigmaVsCrystal[j] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 1000, 0, 100);
    }
    
    for (int j = 0; j < fNumVoltages_gamma; ++j) {
        string histName = "gamma_pulser_Peak_" + to_string(j + 1) + "_Sigma_vs_CrystalID";
        string histTitle = "gamma pulser Peak " + to_string(j + 1) + " Sigma vs Crystal ID; Crystal ID; Sigma";

        gamma_SigmaVsCrystal[j+fNumPeaks] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 1000, 0, 100);
    }
    
    for (int j = 0; j < fNumPeaks; ++j) {
        string histName = "proton_source_Peak_" + to_string(j + 1) + "_Sigma_vs_CrystalID";
        string histTitle = "proton source Peak " + to_string(j + 1) + " Sigma vs Crystal ID; Crystal ID; Sigma";

        proton_SigmaVsCrystal[j] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 1000, 0, 30);
    }
    
    for (int j = 0; j < fNumVoltages_proton; ++j) {
        string histName = "proton_pulser_Peak_" + to_string(j + 1) + "_Sigma_vs_CrystalID";
        string histTitle = "proton pulser Peak " + to_string(j + 1) + " Sigma vs Crystal ID; Crystal ID; Sigma";

        proton_SigmaVsCrystal[j+fNumPeaks] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 1000, 0, 30);
    }
   
   
    //fill histograms
    Int_t Crystalrange_low = xValues.front() - 5;
    Int_t Crystalrange_high = xValues.back() + 5;

    for (int i = 0; i < xValues_gamma.size(); ++i) 
    {
        OffsetGammaVsCrystal->Fill(xValues_gamma[i], yOffsetValues[i]);
        SlopeGammaVsCrystal->Fill(xValues_gamma[i], ySlopeValues[i]);
    }
    
    for (int i = 0; i < xValues_proton.size(); ++i)   
    { 
        OffsetProtonVsCrystal->Fill(xValues_proton[i], yOffsetPulserValues[i]);
        SlopeProtonVsCrystal->Fill(xValues_proton[i], ySlopePulserValues[i]);
        Chi2ProtonVsCrystal->Fill(xValues_proton[i], ychi2[i]);
        PValueProtonVsCrystal->Fill(xValues_proton[i], ypvalue[i]);
    }

    OffsetGammaVsCrystal->SetMarkerStyle(20); 
    OffsetGammaVsCrystal->SetMarkerSize(0.5);
    SlopeGammaVsCrystal->SetMarkerStyle(20);     
    SlopeGammaVsCrystal->SetMarkerSize(0.5);
    PeaksFoundVsCrystal->SetMarkerStyle(20);    
    PeaksFoundVsCrystal->SetMarkerSize(0.5);
    CalibratedEnergies_VoltageVsCrystal->SetMarkerStyle(20);    
    CalibratedEnergies_VoltageVsCrystal->SetMarkerSize(0.5);
    SourcePeaksFoundGammaVsCrystal->SetMarkerStyle(20);
    SourcePeaksFoundGammaVsCrystal->SetMarkerSize(0.5);    
    SourcePeaksFoundProtonVsCrystal->SetMarkerStyle(20);
    SourcePeaksFoundProtonVsCrystal->SetMarkerSize(0.5);     
    PeaksCalibratedGammaVsCrystal->SetMarkerStyle(20);
    PeaksCalibratedGammaVsCrystal->SetMarkerSize(0.5);
    PeaksCalibratedProtonVsCrystal->SetMarkerStyle(20);
    PeaksCalibratedProtonVsCrystal->SetMarkerSize(0.5);    
    OffsetProtonVsCrystal->SetMarkerStyle(20);
    OffsetProtonVsCrystal->SetMarkerSize(0.5);
    SlopeProtonVsCrystal->SetMarkerStyle(20);
    SlopeProtonVsCrystal->SetMarkerSize(0.5); 
    Chi2ProtonVsCrystal->SetMarkerStyle(20);
    Chi2ProtonVsCrystal->SetMarkerSize(0.5);
    PValueProtonVsCrystal->SetMarkerStyle(20);
    PValueProtonVsCrystal->SetMarkerSize(0.5);
    RangeFactorVsCrystal->SetMarkerStyle(20);
    RangeFactorVsCrystal->SetMarkerSize(0.5);
            
    OffsetGammaVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, fNumCrystals/2);
    SlopeGammaVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, fNumCrystals/2);
    PeaksFoundVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
    CalibratedEnergies_VoltageVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);
    SourcePeaksFoundGammaVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, fNumCrystals/2);
    SourcePeaksFoundProtonVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);
    PeaksCalibratedGammaVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, fNumCrystals/2);
    PeaksCalibratedProtonVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);
    SlopeProtonVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high); 
    OffsetProtonVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);
    Chi2ProtonVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);      
    PValueProtonVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);      
    RangeFactorVsCrystal->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);

    PeaksFoundVsCrystal->Write();

    for (int j = 0; j < Num_all_peaks_gamma; ++j) 
    {
        for (int i = 0; i < xValues_gamma.size(); ++i) 
        {
            gamma_SigmaVsCrystal[j]->Fill(xValues_gamma[i], gamma_ySigmaValues[j][i]);
        }
        
        gamma_SigmaVsCrystal[j]->SetMarkerStyle(20);
        gamma_SigmaVsCrystal[j]->SetMarkerSize(0.5);
        
        gamma_SigmaVsCrystal[j]->GetXaxis()->SetRangeUser(Crystalrange_low, fNumCrystals/2);
        
        gamma_SigmaVsCrystal[j]->Write();
        
    }
    
    SourcePeaksFoundGammaVsCrystal->Write();
    OffsetGammaVsCrystal->Write();
    SlopeGammaVsCrystal->Write();
    PeaksCalibratedGammaVsCrystal->Write();
    
    for (int j = 0; j < Num_all_peaks_proton; ++j) 
    {
        for (int i = 0; i < xValues_proton.size(); ++i) 
        {
            proton_SigmaVsCrystal[j]->Fill(xValues_proton[i], proton_ySigmaValues[j][i]);
        }

        proton_SigmaVsCrystal[j]->SetMarkerStyle(20);
        proton_SigmaVsCrystal[j]->SetMarkerSize(0.5);
        
        proton_SigmaVsCrystal[j]->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);
        
        proton_SigmaVsCrystal[j]->Write();
        
    }
    
    SourcePeaksFoundProtonVsCrystal->Write();
    OffsetProtonVsCrystal->Write();
    SlopeProtonVsCrystal->Write();
    PeaksCalibratedProtonVsCrystal->Write();
    Chi2ProtonVsCrystal->Write();
    PValueProtonVsCrystal->Write();
    CalibratedEnergies_VoltageVsCrystal->Write();
    
    for (int j = 0; j < (fNumVoltages_proton-fNumVoltages_gamma); ++j) 
    {
        CalibratedEnergies_differenceVsCrystal[j]->SetMarkerStyle(20);
        CalibratedEnergies_differenceVsCrystal[j]->SetMarkerSize(0.5);
        
        CalibratedEnergies_differenceVsCrystal[j]->GetXaxis()->SetRangeUser(fNumCrystals/2, Crystalrange_high);
        
        CalibratedEnergies_differenceVsCrystal[j]->Write();
        
    } 
    
    RangeFactorVsCrystal->Write();
    
    //_________________________________________//
    
    //set CalParameters
    fCal_Par->setChanged();
    
    FairParAsciiFileIo *parOut = new FairParAsciiFileIo();
	parOut->open(fCalName, "out");
	
	rtdb->setOutput(parOut);
    rtdb->saveOutput();
    
    spectrum_file_gamma->Close();
    spectrum_file_proton->Close();
    
    local_outputFile->Close();
    
    crystals_errors.close();
    

    
    return;
}


Double_t R3BCalifaMapped2CrystalCalPar::FindChisquare(Double_t* X, Double_t* Y, Double_t* eX, Int_t ndf, TF1* f)
{
    cout<<"FindChisquare() called";
    Double_t chi2 = 0.;
    for (Int_t i=0; i<ndf; i++)
    {
	  // gaussian weighted average
	  chi2 += TMath::Exp(-eX[i]*eX[i])*TMath::Sq(f->Eval(X[i])-Y[i])/TMath::Sq(Y[i]);
    }

    return chi2/ndf;
}

ClassImp(R3BCalifaMapped2CrystalCalPar)
