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
    , fSourceName("fitting")
    , fMaxSigma(50.0)
    , fMinPeakEvents(100)
    , fPulserNumber(3)
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
    if (fSourceName == "spectrum")	// fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe" || fSourceName == "152Eu" || fSourceName == "22Na_pulser" || fSourceName == "60Co_pulser" || fSourceName == "AmBe_pulser"  || fSourceName == "152Eu_pulser" || 
    {
        fh_Map_energy_crystal = new TH1F*[fNumCrystals];
    } else if (fSourceName == "fitting")
    {
        fh2_peak_cal = new TH2F*[fNumCrystals];
	fh2_residual_energy = new TH2F*[fNumCrystals];
	fh2_sig_crystal = new TH2F*[3]; // 3 sources
	fh_Map_fit = new TH1F*[3*fNumCrystals];
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
	    else if (fSourceName=="fitting")
	    {
	        fh2_peak_cal[i] = new TH2F(name2, name2, fbins, fleft, fright, 250, 0, 5000);
		fh2_residual_energy[i] = new TH2F(name6, name6, 500, 0, 5000, 100, -50, 50);
		fh_Map_fit[i] = new TH1F(name3,name3,fbins,fleft,fright);
		fh_Map_fit[fNumCrystals+i] = new TH1F(name4,name4,fbins,fleft,fright);
		fh_Map_fit[2*fNumCrystals+i] = new TH1F(name5,name5,fbins,fleft,fright);
	    }
        }

    if (fSourceName == "spectrum")
    {
        fh2_Map_crystal_gamma = new TH2F("fh2_Map_crystal_gamma","fh2_Map_crystal_gamma;crystal ID;Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh2_Map_crystal_proton = new TH2F("fh2_Map_crystal_proton","fh2_Map_crystal_proton;crystal ID;Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
        fh_peak_crystal_gamma = new TH2F("fh_peak_crystal_gamma","fh_peak_crystal_gamma;crystal ID;peak Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh_peak_crystal_proton = new TH2F("fh_peak_crystal_proton","fh_peak_crystal_proton;crystal ID;peak Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
        fh_sigma_crystal_gamma = new TH2F("fh_sigma_crystal_gamma","fh_sigma_crystal_gamma;crystal ID;sigma Map Energy",fNumCrystals/2,0,fNumCrystals/2,15*fSigma,0,15*fSigma);
        fh_sigma_crystal_proton = new TH2F("fh_sigma_crystal_proton","fh_sigma_crystal_proton;crystal ID;sigma Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,15*fSigma,0,15*fSigma);
    } 
    else if (fSourceName == "fitting")
    {
	fh2_resolution_crystalID  = new TH2F("fh2_resolution_crystalID",";crystal ID;resolution (percentage)",fNumCrystals,0,fNumCrystals,40,0,20);
        fh2_slope_crystalID = new TH2F("fh2_slope_crystalID",";crystal ID;slope",fNumCrystals,0,fNumCrystals,40,0,20);
	fh2_intercept_crystalID = new TH2F("fh2_intercept_crystalID",";crystal ID;intercept",fNumCrystals,0,fNumCrystals,50,-100,50);
	fh2_sig_crystal[0] = new TH2F("fh2_sig_crystal_22Na","22Na;crystal ID;sigma",fNumCrystals,0,fNumCrystals,200,0,100);
	fh2_sig_crystal[1] = new TH2F("fh2_sig_crystal_60Co","60Co;crystal ID;sigma",fNumCrystals,0,fNumCrystals,200,0,100);
	fh2_sig_crystal[2] = new TH2F("fh2_sig_crystal_AmBe","AmBe;crystal ID;sigma",fNumCrystals,0,fNumCrystals,200,0,100);
	fh2_chi2_crystal = new TH2F("fh2_chi2_crystal",";crystal ID;Log10 weighted sum of squares",fNumCrystals,0,fNumCrystals,100,-20,0);
	fh_numPeak = new TH1F("fh_numPeak","fh_numPeak",10,0,10);
    }

    return kSUCCESS;
}


InitStatus R3BCalifaMapped2CrystalCalPar::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

//fill Histogramms in spectrum.root when function spectrum is selected, otherwise use already existing spectrum.root file
void R3BCalifaMapped2CrystalCalPar::Exec(Option_t* opt)
{
    if (fSourceName != "spectrum") return;	//???????????????????????????????????????????????????????????????????????????
    
    // cout << "fill spectrum.root" << endl;

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
	    } else
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
	        } else
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
    
    if (fSourceName == "22Na" || fSourceName == "60Co" || fSourceName == "AmBe"  || fSourceName == "152Eu" )
    {
	SearchPeaks();
    }
    
    if (fSourceName == "22Na_pulser" || fSourceName == "60Co_pulser" || fSourceName == "AmBe_pulser"  || fSourceName == "152Eu_pulser" )
    {
	PulserCalibration();
	fCal_Par->printParams();
    }
    
    if (fSourceName == "fitting")
    {
	FitPeaks();
    	fCal_Par->printParams();
    }

    outrootfile->Write();
    outrootfile->Close();
}



//_____________Search_Peaks____________________//

void R3BCalifaMapped2CrystalCalPar::SearchPeaks() 
{
    for (Int_t i=0; i<fNumCrystals; i++)
    {
	Int_t nfound = 0;
	TSpectrum* ss = new TSpectrum(fNumPeaks);

	if ((fMap_Par->GetInUse(i+1) == 1) && (fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics))
	{
	    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "", fThreshold); // "goff" to turn off drawing
            fh_Map_energy_crystal[i]->Write();

            fChannelPeaks = (Double_t*) ss->GetPositionX();
	    Int_t idx[nfound];
            TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);

	    for (Int_t j=0; j<nfound; j++)
	    {
		// gaussian fit
		Double_t posX = fChannelPeaks[idx[nfound-j-1]];
    		TF1 * gaussfit;
    
    		if (i<fNumCrystals/2)
			gaussfit = new TF1("gaussfit","gaus",posX-fSigma*15,posX+fSigma*15);
    		else
        		gaussfit = new TF1("gaussfit","gaus",posX-fSigma*1.5,posX+fSigma*1.5);
		
		TH1F* h_copy = (TH1F*) fh_Map_energy_crystal[i]->Clone("h_copy");
		fh_Map_energy_crystal[i]->Fit("gaussfit","RQ");
		Double_t mean = gaussfit->GetParameter(1);
		Double_t sigma = gaussfit->GetParameter(2);
	
		Double_t pmX[1] = {(Double_t) i};
		Double_t pmY[1] = {mean};
		Double_t pmZ[1] = {sigma};
		
		TPolyMarker *pm1 = new TPolyMarker(1,pmX,pmY);
		pm1->SetMarkerStyle(23);
		pm1->SetMarkerColor(kRed);
		TPolyMarker *pm2 = new TPolyMarker(1,pmX,pmZ);
		pm2->SetMarkerStyle(23);
		pm2->SetMarkerColor(kBlack);

		if (i<fNumCrystals/2)
		{	
			fh_peak_crystal_gamma->Fill(i,mean);
			//fh_peak_crystal_gamma->GetListOfFunctions()->Add(pm1);
			//fh_peak_crystal_gamma->GetListOfFunctions()->Print();
			fh_sigma_crystal_gamma->Fill(i,sigma);
			fh2_Map_crystal_gamma->GetListOfFunctions()->Add(pm1);
			fh2_Map_crystal_gamma->GetListOfFunctions()->Print();
			fh_sigma_crystal_gamma->GetListOfFunctions()->Add(pm2);
			fh_sigma_crystal_gamma->GetListOfFunctions()->Print();
		}
    		else
		{
			fh_peak_crystal_proton->Fill(i,mean);
			//fh_peak_crystal_proton->GetListOfFunctions()->Add(pm1);
			//fh_peak_crystal_proton->GetListOfFunctions()->Print();
			fh_sigma_crystal_proton->Fill(i,sigma);
			fh2_Map_crystal_proton->GetListOfFunctions()->Add(pm1);
			fh2_Map_crystal_proton->GetListOfFunctions()->Print();
			fh_sigma_crystal_proton->GetListOfFunctions()->Add(pm2);
			fh_sigma_crystal_proton->GetListOfFunctions()->Print();
		}
	    }
	}

  	if (ss)
		delete ss;
    }

    fh2_Map_crystal_gamma->Write("colz");
    fh2_Map_crystal_proton->Write("colz");
    fh_peak_crystal_gamma->Write();
    fh_peak_crystal_proton->Write();
    fh_sigma_crystal_gamma->Write();
    fh_sigma_crystal_proton->Write();
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
    TFile* local_outrootfile = TFile::Open(fSpectrumName, "READ");
    if (!local_outrootfile) 
    {
        std::cerr << "No spectrum.root file found. Please perform histogram creation first!" << std::endl;
        return;
    }
    else
    {
        cout << "use spectrum file: " << fSpectrumName << endl;
    }
  
    //open calibrated.root file
    TFile* local_outputFile = TFile::Open(foutputName, "UPDATE");
    if (!local_outputFile) 
    {
        std::cerr << "No calibrated.root file " << foutputName << " found. Please perform histogram creation first!" << std::endl;
        return;
    }
    else
    {
        cout << "use calibrated file: " << foutputName << endl;
    }
  
    //___________________________________________________calibration____________________________________________________________________//
  
    TH1F* local_fh_Map_energy_crystal[fNumCrystals];

    Int_t numPars = 2; // Number of parameters=2 by default
    Int_t nfound = 0;
    if (fNumParam)  // if num of parameters is explicitly specified then get this value
    {
        numPars = fNumParam;
    }
    
    fCal_Par->SetNumCrystals(fNumCrystals);
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(numPars * fNumCrystals);
    Int_t fright, fleft;
    TSpectrum* ss = new TSpectrum(fNumPeaks + fPulserNumber);	//fNumPeaks corresponds to the number of energy values ​​specified in the macro, fPulserNumber corresponds to the number of different pulser signals
    
    for (Int_t i=0; i<fNumCrystals; i++)
    {
    
        char histName[100];
        sprintf(histName, "fh_Map_energy_crystal_spectrum_%i", i + 1);

        local_fh_Map_energy_crystal[i] = (TH1F*)local_outrootfile->Get(histName);


        if (local_fh_Map_energy_crystal[i]) 
        {
            if (local_fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics) 
            {
        
	            nfound = ss->Search(local_fh_Map_energy_crystal[i], fSigma, "", fThreshold);	// "goff" to turn off drawing
                
                //safe histrogramms in outputFile    
	            local_outputFile->cd();
                local_fh_Map_energy_crystal[i]->Write();

                fChannelPeaks = (Double_t*) ss->GetPositionX();
	            Int_t idx[nfound];
                TMath::Sort(nfound, fChannelPeaks, idx, kFALSE);	//kFALSE: sort in ascending order (from lowest to highest) -> sorted indexes: idx (kTRUE: in decscending)

                // Calibrated Spectrum
                Double_t X[nfound + 1];	//array with size nfound+1 elements
                Double_t Y[nfound + 1];

                cout << "Crystal " << i + 1 << endl;  
                    
                if ((nfound - fPulserNumber) == fNumPeaks)	//Check that the number of source peaks found matches the expected number
	            {    
	        	    //gauss fit each found peak and check sigma and events
		            for (Int_t j=0; j<nfound; j++)
		            {
			        Double_t posX = fChannelPeaks[idx[j]];
		                X[j] = posX;
		            
		            if (j<fNumPeaks) {
		            	Y[j] = fEnergyPeaks->GetAt(j);	//fills the Y array with reference energy values from fEnergyPeaks
			        }
			        
			        cout << "Peak " << j + 1 << ", uncalibrated found energy: " << X[j] << endl;
			        
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
			        
		            Double_t peakArea = peakHeight * sigma * sqrt(2 * M_PI);

			        if (sigma > fMaxSigma) // Check if the standard deviation is too large
		                LOG(error)  << "Peak " << j + 1 << ": Sigma too large (" << sigma << ")!";
		            
		            if (peakArea < fMinPeakEvents) // Check if the number of events at the peak is too small
		                LOG(error)  << "Peak " << j + 1 << ": Too few events at peak (" << peakArea << ")!";
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
	                
		            cout << "Fit X, Y Werte: " << X[0] << ", " << X[1] << ", " << Y[0] << ", " << Y[1] << endl;

		            //pass slope and offset
	                for (Int_t h = 0; h < numPars; h++)
	                {
	                    fCal_Par->SetCryCalParams(f1fit->GetParameter(h), numPars*i+h); //1-base
	                
	                    Double_t parameterValue = f1fit->GetParameter(h);
                        fCal_Par->SetCryCalParams(parameterValue, numPars*i+h);
                        cout << "Parameter " << numPars*i+h << " set to " << parameterValue << endl;

	                }
  
                }
                else
                {
                	cout << "Number of peaks found does not correspond to the expected number!!!!" << endl;
                }
	            
	        }     
	    }
    }
    
    delete ss;
    
    fCal_Par->setChanged();
    
    FairParAsciiFileIo *parOut = new FairParAsciiFileIo();
	parOut->open(fCalName, "out");
	
	rtdb->setOutput(parOut);
    rtdb->saveOutput();
     
    
    local_outrootfile->Close();
    local_outputFile->Close();
    

    
    return;
}



//_____________Fit_Peaks____________________//

void R3BCalifaMapped2CrystalCalPar::FitPeaks()
{

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
