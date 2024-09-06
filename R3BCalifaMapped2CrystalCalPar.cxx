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
        if (fMap_Par->GetInUse(i + 1) == 1)
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
        if (fMap_Par->GetInUse(crystalId) == 1)
        {
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
	    }	
	
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
    
        std::cerr << "No gamma_spectrum.root file found. Please perform histogram creation first!" << std::endl;
        return;
    }
    else
    {
        cout << "use gamma_spectrum file: " << fSpectrumName_gamma << endl;
    }
    
    TFile* spectrum_file_proton = TFile::Open(fSpectrumName_proton, "READ");
    if (!spectrum_file_proton) 
    {
    
        std::cerr << "No proton_spectrum.root file found. Please perform histogram creation first!" << std::endl;
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
    Double_t Num_all_peaks_gamma = fNumVoltages_gamma+fNumPeaks;
    Double_t Num_all_peaks_proton = fNumVoltages_proton+fNumPeaks; 
    
    // create histograms
    TH1F* local_fh_Map_energy_crystal[fNumCrystals];

    TH2F* OffsetPulserVsCrystal = new TH2F("OffsetPulser_vs_CrystalID", "OffsetPulser vs Crystal ID; Crystal ID; Offset", fNumCrystals, 0, fNumCrystals, 30000, -15000, 15000);
    TH2F* SlopePulserVsCrystal = new TH2F("SlopePulser_vs_CrystalID", "SlopePulser vs Crystal ID; Crystal ID; Slope", fNumCrystals, 0, fNumCrystals, 500, 0, 50);
    TH2F* OffsetVsCrystal = new TH2F("Offset_vs_CrystalID", "Offset vs Crystal ID;Crystal ID;Offset", fNumCrystals, 0, fNumCrystals, 1200, -80, 40);
    TH2F* SlopeVsCrystal = new TH2F("Slope_vs_CrystalID", "Slope vs Crystal ID;Crystal ID;Slope", fNumCrystals, 0, fNumCrystals, 2000, 0, 20);
    TH2F* PeaksFoundVsCrystal = new TH2F("PeaksFound_vs_Crystal", "PeaksFound vs Crystal ID; Crystal ID; PeaksFound", fNumCrystals, 0, fNumCrystals, 30000, 0, 30000);
    TH2F* PeaksFoundCalibratedGammaVsCrystal = new TH2F("PeaksFoundCalibratedGamma_vs_Crystal", "PeaksFoundCalibratedGamma vs Crystal ID; Crystal ID; PeaksFoundCalibratedGamma", fNumCrystals, 0, fNumCrystals, 30000, 0, 30000);
    TH2F* PeaksFoundCalibratedProtonVsCrystal = new TH2F("PeaksFoundCalibratedProton_vs_Crystal", "PeaksFoundCalibratedProton vs Crystal ID; Crystal ID; PeaksFoundCalibratedProton", fNumCrystals, 0, fNumCrystals, 50000, 0, 500000);
    
    
    vector<TH2F*> gamma_SigmaVsCrystal(Num_all_peaks_gamma);
    vector<TH2F*> proton_SigmaVsCrystal(Num_all_peaks_gamma);
    
    for (int j = 0; j < Num_all_peaks_gamma; ++j) {
        string histName = "gamma_Peak_" + to_string(j + 1) + "_Sigma_vs_CrystalID";
        string histTitle = "gamma_Peak " + to_string(j + 1) + " Sigma vs Crystal ID; Crystal ID; Sigma";

        gamma_SigmaVsCrystal[j] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 1000, 0, 100);
    }
    
    for (int j = 0; j < Num_all_peaks_proton; ++j) {
        string histName = "proton_Peak_" + to_string(j + 1) + "_Sigma_vs_CrystalID";
        string histTitle = "proton_Peak " + to_string(j + 1) + " Sigma vs Crystal ID; Crystal ID; Sigma";

        proton_SigmaVsCrystal[j] = new TH2F(histName.c_str(), histTitle.c_str(), fNumCrystals, 0, fNumCrystals, 1000, 0, 100);
    }
    
    
    // Arrays to store the x and y values ​​for the graphs
    vector<double> xValues, xValues_gamma, xValues_proton, yOffsetValues, ySlopeValues, yOffsetPulserValues, ySlopePulserValues;

    vector<vector<double>> gamma_ySigmaValues(Num_all_peaks_gamma, vector<double>());
    vector<vector<double>> proton_ySigmaValues(Num_all_peaks_proton, vector<double>());
    
    
    Int_t numPars = 2; // Number of parameters=2 by default
    Int_t nfound = 0;
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
    
    	//fNumPeaks corresponds to the number of energy values ​​specified in the macro, fNumVoltages corresponds to the number of different pulser signals
    
    for (Int_t i=0; i<fNumCrystals; i++)
    {    

        Int_t Num_all_peaks = (i < fNumCrystals / 2) ? Num_all_peaks_gamma : Num_all_peaks_proton;

        TSpectrum* ss = new TSpectrum(Num_all_peaks);

        char histName[100];
        sprintf(histName, "fh_Map_energy_crystal_spectrum_%i", i + 1);

        if (i<fNumCrystals/2)
            local_fh_Map_energy_crystal[i] = (TH1F*)spectrum_file_gamma->Get(histName);
        else
            local_fh_Map_energy_crystal[i] = (TH1F*)spectrum_file_proton->Get(histName);
        

        if ((local_fh_Map_energy_crystal[i]) && (local_fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics))
        {
            nfound = ss->Search(local_fh_Map_energy_crystal[i], fSigma, "", fThreshold);	// "goff" to turn off drawing
            
            //safe histograms in outputFile    
            local_outputFile->cd();
            local_fh_Map_energy_crystal[i]->Write();

            fChannelPeaks = (Double_t*) ss->GetPositionX();
            Int_t idx[nfound];
            TMath::Sort(nfound, fChannelPeaks, idx, kFALSE);	//kFALSE: sort in ascending order (from lowest to highest) -> sorted indexes: idx (kTRUE: in decscending)

            // Calibrated Spectrum
            Double_t X[nfound + 1];	//array with size nfound+1 elements
            Double_t Y[nfound + 1];

            cout << "Crystal " << i + 1 << endl;  
                
            if (nfound == Num_all_peaks)	//Check that the number of source peaks found matches the expected number
            {          
                xValues.push_back(i + 1);
            
                crystals_errors << i + 1;
                
        	    //gauss fit each found peak and check sigma and events
                for (Int_t j=0; j<Num_all_peaks; j++)
                {
                    Double_t posX = fChannelPeaks[idx[j]];
                    X[j] = posX;
                        
                    PeaksFoundVsCrystal->Fill(i + 1, posX);
                    
                    if (j<fNumPeaks) 
                    {
                    	Y[j] = fEnergyPeaks->GetAt(j);	//fills the Y array with reference energy values from fEnergyPeaks
                    	
                    	cout << "Source Peak " << j + 1 << ", bin number: " << X[j] << endl;
	                }
	                else
	                {
	                    cout << "Pulser Peak " << j + 1 - fNumPeaks << ", bin number: " << X[j] << endl;
	                }
	                
	                // Error checks
                    TF1 * gaussfit;
                    
                    if (i<fNumCrystals/2) {
		                gaussfit = new TF1("gaussfit","gaus",posX-fSigma*15,posX+fSigma*15);
		            }
                    else {
                    	gaussfit = new TF1("gaussfit","gaus",posX-fSigma*1.5,posX+fSigma*1.5);
                    }
	                TH1F* h_copy = (TH1F*) local_fh_Map_energy_crystal[i]->Clone("h_copy");
	                local_fh_Map_energy_crystal[i]->Fit("gaussfit","RQ");
	                
	                Double_t peakHeight = gaussfit->GetParameter(0);
	                Double_t mean = gaussfit->GetParameter(1);
	                Double_t sigma = gaussfit->GetParameter(2);
	                
	                
	                if (i < fNumCrystals / 2)
	                    gamma_ySigmaValues[j].push_back(sigma);
                    else
                        proton_ySigmaValues[j].push_back(sigma);
                    
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
                    
                TF1* f1fit = nullptr;
                
                if(fNumParam)
                {
                    if (fNumParam == 1) 
                        f1fit = new TF1("f1fit", "[0]*x", fleft, fright);
                    if (fNumParam == 2)
                        f1fit = new TF1("f1fit", "[0]+[1]*x", fleft, fright);
                }
                else
                {
                    LOG(warn)
                        << "R3BCalifaMapped2CrystalCalPar:: No input number of fit parameters, therefore, by default "
                           "NumberParameters=2";
                    f1fit = new TF1("f1fit", "[0]+[1]*x", fleft, fright);
                }
                
                TGraph* graph = new TGraph(fNumPeaks, X, Y);
                graph->Fit("f1fit", "Q"); // Quiet mode (minimum printing)
                
                cout << "Fit X values: ";
                for (int a = 0; a < fNumPeaks; ++a) {
                    cout << X[a];
                    if (a < fNumPeaks - a) 
                    {
                        cout << ", ";
                    }
                }
                cout << endl;
                
                cout << "Fit Y values: ";
                for (int b = 0; b < fNumPeaks; ++b) {
                    cout << Y[b];
                    if (b < fNumPeaks - b) 
                    {
                        cout << ", ";
                    }
                }
                cout << endl;
                
                //get offset and slope
                Double_t CalParams[fNumParam];  //standard, [0]: offset, [1]: slope
                
                for (Int_t h = 0; h < fNumParam; h++)
                {
                    CalParams[h] = f1fit->GetParameter(h);
                }

                cout << "Offset: " << CalParams[0] << endl;                
                cout << "Slope: " << CalParams[1] << endl;

                // Fill the histograms
                yOffsetValues.push_back(CalParams[0]);
                ySlopeValues.push_back(CalParams[1]);

                //_____________________________pulser_calibration_____________________________________//

                

                //_________________pulser_gamma_________________________//
                
                if (i<fNumCrystals/2)
                {  
                    xValues_gamma.push_back(i + 1);
                    
                    Double_t PulserEnergycalibrated[fNumVoltages_gamma];
                    
                    Double_t PulserEnergyexpected[fNumVoltages_gamma];
                    Double_t pulser_sigma = fSigma/100;
                    
                    
                    for (Int_t v = 0; v < fNumVoltages_gamma; v++)
                    {
                        PulserEnergycalibrated[v] = CalParams[0] + X[fNumPeaks + v] * CalParams[1];
                        cout << v + 1 << " Pulser energy calibrated: " << PulserEnergycalibrated[v] << endl;
                        
                        PeaksFoundCalibratedGammaVsCrystal->Fill(i + 1, PulserEnergycalibrated[v]);
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
                    xValues_proton.push_back(i + 1);
                    
                    Double_t PulserEnergycalibrated[fNumVoltages_proton];
                    
                    PulserEnergycalibrated[0] = CalParams[0] + X[fNumPeaks] * CalParams[1];
                    
                    cout << "1 Pulser energy calibrated: " << PulserEnergycalibrated[0] << endl;
                    
                    PeaksFoundCalibratedProtonVsCrystal->Fill(i + 1, PulserEnergycalibrated[0]);
                    
                    for (Int_t e = 1; e < fNumVoltages_proton; e++)
                    {
                        PulserEnergycalibrated[e] = PulserEnergycalibrated[0] * V_proton[e]/V_proton[0];
                        cout << e + 1 << " Pulser energy calibrated: " << PulserEnergycalibrated[e] << endl;
                        
                        PeaksFoundCalibratedProtonVsCrystal->Fill(i + 1, PulserEnergycalibrated[e]);
                    }
                    
                    Double_t X_Pulser[fNumVoltages_proton];
                    
                    for (Int_t r = 0; r < fNumVoltages_proton; r++)
                    {
                        X_Pulser[r] = X[fNumPeaks + r];
                    }
                    
                    TF1* f2fit = nullptr;
                    f2fit = new TF1("f2fit", "[0]+[1]*x", fleft, fright);
                    
                    TGraph* graph2 = new TGraph(fNumVoltages_proton, X_Pulser, PulserEnergycalibrated);
                    graph2->Fit("f2fit", "Q"); // Quiet mode (minimum printing)
                    
                    Double_t PulserParams[fNumParam];
                    
                    for (Int_t f = 0; f < fNumParam; f++)
                    {
                        PulserParams[f] = f2fit->GetParameter(f);
                    }
                    
                    cout << "Pulser offset: " << PulserParams[0] << endl;
                    cout << "Pulser slope: " << PulserParams[1] << endl;
                       
                    yOffsetPulserValues.push_back(PulserParams[0]);
                    ySlopePulserValues.push_back(PulserParams[1]);
                    
                    //pass slope and offset from pulser calibration
                    for (Int_t p = 0; p < fNumParam; p++)
                    {
                        fCal_Par->SetCryCalParams(PulserParams[p], fNumParam*i+p);

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
    //set histogram boundaries
    Int_t Crystalrange_low = xValues.front() - 5;
    Int_t Crystalrange_high = xValues.back() + 5;

    for (int i = 0; i < xValues.size(); ++i) 
    {
        OffsetVsCrystal->Fill(xValues[i], yOffsetValues[i]);
        SlopeVsCrystal->Fill(xValues[i], ySlopeValues[i]);
        SlopePulserVsCrystal->Fill(xValues[i], ySlopeValues[i]);
        OffsetPulserVsCrystal->Fill(xValues[i], ySlopeValues[i]);
    }

    for (int j = 0; j < Num_all_peaks_gamma; ++j) 
    {
        for (int i = 0; i < xValues_gamma.size(); ++i) 
        {
            gamma_SigmaVsCrystal[j]->Fill(xValues_gamma[i], gamma_ySigmaValues[j][i]);
        }
        
        gamma_SigmaVsCrystal[j]->SetMarkerStyle(20);
        gamma_SigmaVsCrystal[j]->SetMarkerSize(0.5);
        
        gamma_SigmaVsCrystal[j]->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
        
        gamma_SigmaVsCrystal[j]->Write();
        
    }
    
    for (int j = 0; j < Num_all_peaks_proton; ++j) 
    {
        for (int i = 0; i < xValues_proton.size(); ++i) 
        {
            proton_SigmaVsCrystal[j]->Fill(xValues_proton[i], proton_ySigmaValues[j][i]);
        }

        proton_SigmaVsCrystal[j]->SetMarkerStyle(20);
        proton_SigmaVsCrystal[j]->SetMarkerSize(0.5);
        
        proton_SigmaVsCrystal[j]->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
        
        proton_SigmaVsCrystal[j]->Write();
        
    }

    OffsetVsCrystal->SetMarkerStyle(20); 
    OffsetVsCrystal->SetMarkerSize(0.5);
    SlopeVsCrystal->SetMarkerStyle(20);
    SlopeVsCrystal->SetMarkerSize(0.5);
    PeaksFoundVsCrystal->SetMarkerStyle(20);
    PeaksFoundVsCrystal->SetMarkerSize(0.5);
    PeaksFoundCalibratedGammaVsCrystal->SetMarkerStyle(20);
    PeaksFoundCalibratedGammaVsCrystal->SetMarkerSize(0.5);
    PeaksFoundCalibratedProtonVsCrystal->SetMarkerStyle(20);
    PeaksFoundCalibratedProtonVsCrystal->SetMarkerSize(0.5);    
    OffsetPulserVsCrystal->SetMarkerStyle(20);
    OffsetPulserVsCrystal->SetMarkerSize(0.5);
    SlopePulserVsCrystal->SetMarkerStyle(20);
    SlopePulserVsCrystal->SetMarkerSize(0.5); 
        
    OffsetVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
    SlopeVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
    PeaksFoundVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
    PeaksFoundCalibratedGammaVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
    PeaksFoundCalibratedProtonVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);
    SlopePulserVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high); 
    OffsetPulserVsCrystal->GetXaxis()->SetRangeUser(Crystalrange_low, Crystalrange_high);      
    
    OffsetVsCrystal->Write();
    SlopeVsCrystal->Write();
    PeaksFoundVsCrystal->Write();
    PeaksFoundCalibratedGammaVsCrystal->Write();
    PeaksFoundCalibratedProtonVsCrystal->Write();
    
    SlopePulserVsCrystal->Write();
    OffsetPulserVsCrystal->Write();    

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
