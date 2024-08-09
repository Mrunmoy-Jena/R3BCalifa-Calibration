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
    outrootfile = TFile::Open("/media/mrunmoy/MyDisk/data_calibration/calibrated/spectrum.root","UPDATE");
    if (!outrootfile) {
    	outrootfile = new TFile("/media/mrunmoy/MyDisk/data_calibration/calibrated/spectrum.root","RECREATE");
	outrootfile->cd();
	outroottree = new TTree("genT","General Tree");
    }

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
    if (fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe" || fSourceName=="152Eu")
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
	    if (fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe" || fSourceName == "152Eu")
	    {
                fh_Map_energy_crystal[i] = new TH1F(name1, name1, fbins, fleft, fright);
	    } else if (fSourceName=="fitting")
	    {
	        fh2_peak_cal[i] = new TH2F(name2, name2, fbins, fleft, fright, 250, 0, 5000);
		fh2_residual_energy[i] = new TH2F(name6, name6, 500, 0, 5000, 100, -50, 50);
		fh_Map_fit[i] = new TH1F(name3,name3,fbins,fleft,fright);
		fh_Map_fit[fNumCrystals+i] = new TH1F(name4,name4,fbins,fleft,fright);
		fh_Map_fit[2*fNumCrystals+i] = new TH1F(name5,name5,fbins,fleft,fright);
	    }
        }

    if (fSourceName=="60Co" || fSourceName=="22Na" || fSourceName=="AmBe" || fSourceName == "152Eu")
    {
        fh2_Map_crystal_gamma = new TH2F("fh2_Map_crystal_gamma","fh2_Map_crystal_gamma;crystal ID;Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh2_Map_crystal_proton = new TH2F("fh2_Map_crystal_proton","fh2_Map_crystal_proton;crystal ID;Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
        fh_peak_crystal_gamma = new TH2F("fh_peak_crystal_gamma","fh_peak_crystal_gamma;crystal ID;peak Map Energy",fNumCrystals/2,0,fNumCrystals/2,fMapHistos_bins,fMapHistos_left,fMapHistos_right);
        fh_peak_crystal_proton = new TH2F("fh_peak_crystal_proton","fh_peak_crystal_proton;crystal ID;peak Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,fMapHistos_binsp,fMapHistos_leftp,fMapHistos_rightp);
        fh_sigma_crystal_gamma = new TH2F("fh_sigma_crystal_gamma","fh_sigma_crystal_gamma;crystal ID;sigma Map Energy",fNumCrystals/2,0,fNumCrystals/2,15*fSigma,0,15*fSigma);
        fh_sigma_crystal_proton = new TH2F("fh_sigma_crystal_proton","fh_sigma_crystal_proton;crystal ID;sigma Map Energy",fNumCrystals/2,fNumCrystals/2,fNumCrystals,15*fSigma,0,15*fSigma);
    } else if (fSourceName == "fitting")
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

void R3BCalifaMapped2CrystalCalPar::Exec(Option_t* opt)
{
    cout<<"Exec() called";
    if (fSourceName == "fitting") return;

    Int_t nHits = fCalifaMappedDataCA->GetEntries();
    if (!nHits)
        return;

    R3BCalifaMappedData** MapHit = new R3BCalifaMappedData*[nHits];
    Int_t crystalId = 0;

    for (Int_t i = 0; i < nHits; i++)
    {
        MapHit[i] = dynamic_cast<R3BCalifaMappedData*>(fCalifaMappedDataCA->At(i));
        crystalId = MapHit[i]->GetCrystalId();
        // Fill histograms
        if (fMap_Par->GetInUse(crystalId) == 1) {
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

void R3BCalifaMapped2CrystalCalPar::Reset() {cout<<"Reset() called";}

void R3BCalifaMapped2CrystalCalPar::FinishEvent() {cout<<"FinishTask() called";}

void R3BCalifaMapped2CrystalCalPar::FinishTask()
{
    cout<<"finishTask() called";
    if (fSourceName == "22Na" || fSourceName == "60Co" || fSourceName == "AmBe"  || fSourceName == "152Eu" )
    {
	SearchPeaks();
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

// find and record
// mapping level peaks
void R3BCalifaMapped2CrystalCalPar::SearchPeaks()
{
    cout<<"SearchPeaks() called";
    Int_t numPars = 2; // Number of parameters=2 by default
    Int_t nfound = 0;
    if (fNumParam)  // if num of parameters is explicitly specified then get this value
    {
        numPars = fNumParam;
    }
    
    fCal_Par->SetNumCrystals(fNumCrystals);  //set the cal parameters
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(numPars * fNumCrystals); //total num of parameters I guess
    Int_t fright, fleft;
    TSpectrum* ss = new TSpectrum(fNumPeaks);
    
    for (Int_t i=0; i<fNumCrystals; i++)
    {
	//Int_t nfound = 0;
	//TSpectrum* ss = new TSpectrum(fNumPeaks);

	if ((fMap_Par->GetInUse(i+1) == 1) && (fh_Map_energy_crystal[i]->GetEntries() > fMinStadistics))
	{
	    nfound = ss->Search(fh_Map_energy_crystal[i], fSigma, "", fThreshold); // "goff" to turn off drawing
            fh_Map_energy_crystal[i]->Write();

            fChannelPeaks = (Double_t*) ss->GetPositionX();
	    Int_t idx[nfound];
            TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);
            
            // Calibrated Spectrum
            Double_t X[nfound + 1];
            Double_t Y[nfound + 1];
                
	    for (Int_t j=0; j<nfound; j++)
	    {
		// gaussian fit
		//Double_t posX = fChannelPeaks[idx[nfound-j-1]];
                //TF1 * gaussfit;
                X[j] = fChannelPeaks[idx[nfound - j - 1]];
                Y[j] = fEnergyPeaks->GetAt(nfound - j - 1);
                
                 /* if (i<fNumCrystals/2)
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
			fh_sigma_crystal_gamma->Fill(i,sigma);
			fh2_Map_crystal_gamma->GetListOfFunctions()->Add(pm1);
			fh2_Map_crystal_gamma->GetListOfFunctions()->Print();
			fh_peak_crystal_gamma->GetListOfFunctions()->Add(pm1);
			fh_peak_crystal_gamma->GetListOfFunctions()->Print();
			fh_sigma_crystal_gamma->GetListOfFunctions()->Add(pm2);
			fh_sigma_crystal_gamma->GetListOfFunctions()->Print();
		}
                else
		{
			fh_peak_crystal_proton->Fill(i,mean);
			fh_sigma_crystal_proton->Fill(i,sigma);
			fh2_Map_crystal_proton->GetListOfFunctions()->Add(pm1);
			fh2_Map_crystal_proton->GetListOfFunctions()->Print();
			fh_peak_crystal_proton->GetListOfFunctions()->Add(pm1);
			fh_peak_crystal_proton->GetListOfFunctions()->Print();
			fh_sigma_crystal_proton->GetListOfFunctions()->Add(pm2);
			fh_sigma_crystal_proton->GetListOfFunctions()->Print();
		} */
	    }
	    X[nfound]=0.;
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
                    {
                        f1fit = new TF1("f1fit", "[0]*x", fleft, fright);
                    }
                    if (fNumParam == 2)
                    {
                        f1fit = new TF1("f1fit", "[0]+[1]*x", fleft, fright);
                    }
                }
                else
                {
                    LOG(warn)
                        << "R3BCalifaMapped2CrystalCalPar:: No input number of fit parameters, therefore, by default "
                           "NumberParameters=2";
                    f1fit = new TF1("f1fit", "[0]+[1]*x", fleft, fright);
                }
                
                TGraph* graph = new TGraph(fNumPeaks + 1, X, Y);
                graph->Fit("f1fit", "Q"); // Quiet mode (minimum printing)

                for (Int_t h = 0; h < numPars; h++)
                {
                    fCal_Par->SetCryCalParams(f1fit->GetParameter(h), numPars*i+h); //1-base
                }
                
	}
    }
    
    delete ss;
    fCal_Par->setChanged();
    //fh2_Map_crystal_gamma->Write("colz");
    //fh2_Map_crystal_proton->Write("colz");
    //fh_peak_crystal_gamma->Write();
    //fh_peak_crystal_proton->Write();
    //fh_sigma_crystal_gamma->Write();
    //fh_sigma_crystal_gamma->Write();
    //fh_sigma_crystal_proton->Write();
    return;
}

void R3BCalifaMapped2CrystalCalPar::FitPeaks()
{
    cout<<"FitPeaks() called";
    Double_t chi2 = 0.;
    Int_t deleted_peaks = 0;
    fCal_Par->SetNumCrystals(fNumCrystals);
    fCal_Par->SetNumParametersFit(fNumParam);
    fCal_Par->GetCryCalParams()->Set(fNumParam * fNumCrystals);
    Double_t X[fNumCrystals][fMaxPeaks];
    Double_t Y[fNumCrystals][fMaxPeaks];
    Double_t eX[fNumCrystals][fMaxPeaks];
    Double_t eY[fNumCrystals][fMaxPeaks];
    // inference
    TH1F* AmBeSpectra[fNumCrystals];

    for (Int_t i=0; i<fNumCrystals; i++)
    {
	for (Int_t j=0; j<fMaxPeaks; j++)
	{
	    X[i][j] = -1;
	    Y[i][j] = -1;
	    eX[i][j] = -1;
	    eY[i][j] = -1;
	}
    }

    for (auto k : *outrootfile->GetListOfKeys())
    {
	//cout<<"Inside the k loop, k="<<k<<endl;
	TKey *key = static_cast<TKey*>(k);
	TClass *cl = gROOT->GetClass(key->GetClassName());
	std::string title = key->GetName();
	std::string begintitle = title.substr(0,13);
        
        if (begintitle != "fh_Map_energy") continue;
	if (!(cl->InheritsFrom("TH1"))) continue;
        cout<<"Still going?"<<endl;    //Doesn't output, meaning that the above two statements make the flow skip to the next iteration, i.e. nothing after these statements in this for loop runs.. 
	std::string sourceName = title.substr(22,4);
	std::string crystalName = title.substr(27,4);
	Int_t cryId = std::atoi(crystalName.c_str());
	
        //cout<<"Source:"<<sourceName<<endl;
        //cout<<"Crystal name:"<<crystalName<<endl;
        //cout<<"Crystal Id:"<<cryId<<endl;
        //cout<<"Title:"<<title<<endl;   
        //cout<<"Begintitle:"<<begintitle<<endl;
        
	TH1F* h = key->ReadObject<TH1F>();
	fEnergyPeaks->Reset();
	Double_t histMin = h->GetXaxis()->GetXmin();
	Double_t histMax = h->GetXaxis()->GetXmax();
        //cout<<"histMin"<<histMin<<endl;
        //cout<<"histMax"<<histMax<<endl;
        
	Int_t histIdx = 0; // index for fit
        if (sourceName == "22Na")
	{
	  fNumPeaks = 2;
	  fEnergyPeaks->Set(2);
	  fEnergyPeaks->AddAt(1274.5,0);
	  fEnergyPeaks->AddAt(511.0,1);
	  histIdx = cryId-1+fNumCrystals;
	} 
	else if (sourceName == "60Co")
	{
	  fNumPeaks = 2;
	  fEnergyPeaks->Set(2);
	  fEnergyPeaks->AddAt(1332.5,0);
	  fEnergyPeaks->AddAt(1173.2,1);
	  histIdx = cryId-1+fNumCrystals;
	} else if (sourceName == "AmBe")
	{
	  fNumPeaks = 3;
	  fEnergyPeaks->Set(3);
	  fEnergyPeaks->AddAt(4438.0,0);
	  fEnergyPeaks->AddAt(3927.0,1);
	  fEnergyPeaks->AddAt(3416.0,2);
	  histIdx = cryId-1+2*fNumCrystals;

	  AmBeSpectra[cryId-1] = (TH1F*) key->ReadObject<TH1F>();
	} else if (sourceName == "152Eu")
	{
	  fNumPeaks = 1;
	  fEnergyPeaks->Set(1);
	  fEnergyPeaks->AddAt(1408,0);
	  histIdx = cryId-1+3*fNumCrystals;
	}

	// copy the histograms
	fh_Map_fit[histIdx]->SetBins(h->GetNbinsX(),histMin,histMax);
	for (Int_t bin=0; bin<h->GetNbinsX(); bin++)
	{
	    fh_Map_fit[histIdx]->SetBinContent(bin,h->GetBinContent(bin));
	}

	if ((fMap_Par->GetInUse(cryId) == 1) && (h->GetEntries() > fMinStadistics) && fNumPeaks==2)
	{
	    TSpectrum* ss = new TSpectrum(fNumPeaks);
	    TH1F *h2 = (TH1F*)h->Clone("h2");
	    Int_t nfound = ss->Search(h2,fSigma,"goff",fThreshold);
	    TH1F* h_bkg = (TH1F*)ss->Background(h2);
	    fChannelPeaks = (Double_t*) ss->GetPositionX();
	    Int_t idx[nfound];
	    TMath::Sort(nfound, fChannelPeaks, idx, kTRUE);

	    if (nfound != fNumPeaks) continue;

	    Double_t tempX[fNumPeaks];
	    Double_t tempY[fNumPeaks];
	    Double_t tempS[fNumPeaks];

	    for (Int_t tempi=0; tempi<fNumPeaks; tempi++)
	    {
		tempX[tempi] = -1.0;
		tempY[tempi] = -1.0;
		tempS[tempi] = -1.0;
	    }

	    TF1 *fExpo = new TF1("fExpo","expo",h_bkg->GetXaxis()->GetXmin(),h_bkg->GetXaxis()->GetXmax());
	    h_bkg->Fit("fExpo","RQ");

	    TH1F *h3 = (TH1F*) h->Clone("h3");
	    TH1F *h_peak = (TH1F*) h->Clone("h_peak");
	    h_peak->Add(h_bkg,-1);

	    for (Int_t j=0; j<fNumPeaks; j++)
	    {
		    //if (sourceName=="22Na" && j==0) continue;

		    Double_t xPos = fChannelPeaks[idx[fNumPeaks-j-1]];
		    Double_t ratio = TMath::Sqrt(fEnergyPeaks->At(fNumPeaks-j-1)/fGausBaseEnergy);
		    Double_t gRange = 0.;
		    if ((cryId-1)<fNumCrystals/2)
		    {
			    gRange = fGausRange;
		    } else
		    {
			    gRange = fGausRangeP;
		    }
		    Double_t gausMin = TMath::Max(xPos-gRange*ratio,histMin);
		    Double_t gausMax = TMath::Min(xPos+gRange*ratio,histMax);

		    TF1* fGaus1 = new TF1("fGaus1","gaus",gausMin,gausMax);
		    h_peak->Fit("fGaus1","RQ+");

		    TF1* fFit = new TF1("fFit","expo(0)+gaus(2)",gausMin,gausMax);
		    fFit->SetParameter(0,fExpo->GetParameter(0));
		    fFit->SetParameter(1,fExpo->GetParameter(1));
		    fFit->SetParameter(2,fGaus1->GetParameter(0));
		    fFit->SetParameter(3,fGaus1->GetParameter(1));
		    fFit->SetParameter(4,fGaus1->GetParameter(2));
		    h3->Fit("fFit","RQ+");

		    TF1* fExpBkg = new TF1("fExpBkg","expo",gausMin,gausMax);
		    fExpBkg->SetParameter(0,fFit->GetParameter(0));
		    fExpBkg->SetParameter(1,fFit->GetParameter(1));

		    fGaus1->SetLineColor(kBlack);
		    fExpBkg->SetLineColor(kBlack);

		    if (sourceName=="22Na")
		    {
			    fh_Map_fit[cryId-1]->GetListOfFunctions()->Add(fExpBkg);
			    fh_Map_fit[cryId-1]->GetListOfFunctions()->Add(fGaus1);
			    fh_Map_fit[cryId-1]->GetListOfFunctions()->Add(fFit);
		    } 
		     else if (sourceName=="60Co")
		    {
			    fh_Map_fit[cryId-1+fNumCrystals]->GetListOfFunctions()->Add(fExpBkg);
			    fh_Map_fit[cryId-1+fNumCrystals]->GetListOfFunctions()->Add(fGaus1);
			    fh_Map_fit[cryId-1+fNumCrystals]->GetListOfFunctions()->Add(fFit);

			     if (j==1)
			     {
				   // fill resolution plot
				   Double_t pmX[1] = {(Double_t) cryId};
				   Double_t pmY[1] = {100.*fFit->GetParameter(4)/fFit->GetParameter(3)};
				   TPolyMarker *pm = new TPolyMarker(1,pmX,pmY);
				   pm->SetMarkerStyle(23);
				   pm->SetMarkerColor(kBlack);
				   fh2_resolution_crystalID->GetListOfFunctions()->Add(pm);
			     }
		    }

		    tempX[j] = fFit->GetParameter(3);
		    tempY[j] = fEnergyPeaks->At(fNumPeaks-j-1);
		    tempS[j] = fFit->GetParError(3);
	    }

	    int start = 0;
	    while (X[cryId-1][start]>0.) start++;

	    Int_t curIdx = 0;
	    for (Int_t j=0; j<fNumPeaks; j++)
	    {
		if (tempX[j]>0.)
		{
		    X[cryId-1][start+curIdx] = tempX[j];
		    Y[cryId-1][start+curIdx] = tempY[j];
		    eX[cryId-1][start+curIdx] = tempS[j];
		    eY[cryId-1][start+curIdx] = 0.;
		    curIdx++;
		}
	    }
	}
    } 

    // Fit with graph
    Int_t fleft, fright;
    for (Int_t i=0; i<fNumCrystals; i++)
    {
	if (fMap_Par->GetInUse(i+1) == 1)
	{
	    if (i<fMap_Par->GetNumCrystals()/2)
	    {
		fright = fMapHistos_right;
		fleft = fMapHistos_left;
	    } else
	    {
		fright = fMapHistos_rightp;
		fleft = fMapHistos_leftp;
	    }

	    if (fNumParam == 2)
	    {
		f1 = new TF1("f1","[0]+[1]*x",fleft,fright);
	    }

	    Int_t numPeak = 0;
	    while (X[i][numPeak]>0.)
	    {
		    numPeak++;
	    }

	    if (numPeak >= 2)
	    {
		// fit for function
		auto graph = new TGraphErrors(numPeak,X[i],Y[i],eX[i],eY[i]);
		Int_t crystalID = i+1;
		graph->Fit("f1","WQN");
		graph->Fit("f1","MQ+");

		chi2 = TMath::Log10(FindChisquare(X[i],Y[i],eX[i],numPeak,f1));

		if (chi2 > fChi2Threshold)
		{
		    // eliminate calibration points that are off
		    Int_t curPeak = 0;
		    while (curPeak<fMaxPeaks)
		    {
			if (X[i][curPeak]>0.)
			{
			    Double_t Chi2Point = TMath::Exp(-eX[i][curPeak]*eX[i][curPeak])*TMath::Sq(Y[i][curPeak]-f1->Eval(X[i][curPeak]))/TMath::Sq(Y[i][curPeak]);

			    if (chi2 > fChi2Threshold)
			    {
				// eliminate the point
				for (Int_t ipeak=curPeak; ipeak<fMaxPeaks-1; ipeak++)
				{
				    X[i][ipeak] = X[i][ipeak+1];
				    Y[i][ipeak] = Y[i][ipeak+1];
				    eX[i][ipeak] = eX[i][ipeak+1];
				    eY[i][ipeak] = eY[i][ipeak+1];
				}
				X[i][fMaxPeaks-1] = 0.;
				Y[i][fMaxPeaks-1] = 0.;
				eX[i][fMaxPeaks-1] = 0.;
				eY[i][fMaxPeaks-1] = 0.;
				curPeak--;
				numPeak--;
			    }
			}
			curPeak++;
		    }

		    if (numPeak>=2)
		    {
		    	auto modGraph = new TGraphErrors(numPeak,X[i],Y[i],eX[i],eY[i]);
		    	modGraph->Fit("f1","WQN");
		    	modGraph->Fit("f1","MQ+");

			chi2 = TMath::Log10(FindChisquare(X[i],Y[i],eX[i],numPeak,f1));
		    } else
		    {
			    f1->SetParameter(0,0.);
			    f1->SetParameter(1,0.);
			    numPeak = 0;
		    }
		}

		if (chi2 < fChi2Threshold && AmBeSpectra[i] && numPeak>=2)
		{
		  Double_t intercept = f1->GetParameter(0); // y position of f1 at 0
		  Double_t slope = f1->GetParameter(1);

	          Double_t AmBePeaks[3] = {3416.0, 3927.0, 4438.0};

		  Double_t gRange = 0.;
		  if (i<fNumCrystals/2)
		  {
	              gRange = fGausRange;
		  } else
		  {
		      gRange = fGausRangeP;
		  }

	          TH1F* h4 = (TH1F*) AmBeSpectra[i]->Clone("h4");
		  Double_t fitMin = h4->GetXaxis()->GetXmin();
		  Double_t fitMax = h4->GetXaxis()->GetXmax();
		  TF1 *fExpo = new TF1("fExpo","expo",fitMin,fitMax);
		  h4->Fit("fExpo","RQ");

		  TH1F* h_bkg = new TH1F("h_bkg","h_bkg",h4->GetNbinsX(),fitMin,fitMax);
		  for (Int_t bin=0; bin<h4->GetNbinsX(); bin++)
		  {
			h_bkg->SetBinContent(bin,fExpo->Eval((bin+0.5)*h4->GetXaxis()->GetBinWidth(0)+fitMin));
		  }
		  TH1F* h_peak = (TH1F*) AmBeSpectra[i]->Clone("h_peak");
		  h_peak->Add(h_bkg,-1);

		  for (Int_t t=0; t<3; t++)
		  {
	              Double_t expPeak = (AmBePeaks[t]-intercept)/slope;
		      Double_t gaus2min = TMath::Max(expPeak-gRange*TMath::Sqrt(AmBePeaks[t]/fGausBaseEnergy),h4->GetXaxis()->GetXmin());
		      Double_t gaus2max = TMath::Min(expPeak+gRange*TMath::Sqrt(AmBePeaks[t]/fGausBaseEnergy),h4->GetXaxis()->GetXmax());
		      TF1 *fGaus2 = new TF1("fGaus2","gaus",gaus2min,gaus2max);
		      fGaus2->SetParLimits(1,gaus2min,gaus2max);
		      h_peak->Fit("fGaus2","RQ+");

		      Double_t height = fGaus2->GetParameter(0);
		      Double_t mean = fGaus2->GetParameter(1);
		      Double_t sigma = fGaus2->GetParameter(2);

		      TF1 *fFit = new TF1("fFit","expo(0)+gaus(2)",gaus2min,gaus2max);
		      fFit->SetParameter(0,fExpo->GetParameter(0));
		      fFit->SetParameter(1,fExpo->GetParameter(1));
		      fFit->SetParameter(2,height);
		      fFit->SetParameter(3,mean);
		      fFit->SetParameter(4,sigma);
		      fFit->SetParLimits(1,2*fExpo->GetParameter(1),0.);
		      fFit->SetParLimits(2,0.,20*height);
		      fFit->SetParLimits(3,gaus2min,gaus2max);
		      fFit->SetParLimits(4,0.5*sigma,1.5*sigma);

		      h4->Fit("fFit","RQ+");
		      TF1 *fExpBkg = new TF1("fExpBkg","expo",gaus2min,gaus2max);
		      fExpBkg->SetParameter(0,fFit->GetParameter(0));
		      fExpBkg->SetParameter(1,fFit->GetParameter(1));

		      fExpBkg->SetLineColor(kBlack);
		      fGaus2->SetLineColor(kBlack);

		      fh_Map_fit[i+2*fNumCrystals]->GetListOfFunctions()->Add(fExpBkg);
		      fh_Map_fit[i+2*fNumCrystals]->GetListOfFunctions()->Add(fGaus2);
		      fh_Map_fit[i+2*fNumCrystals]->GetListOfFunctions()->Add(fFit);

		      Double_t xPos = fFit->GetParameter(3);
		      Double_t sigY = fFit->GetParameter(2);
		      Double_t bkgY = fFit->Eval(xPos) - sigY;
		      Double_t significance = sigY/TMath::Sqrt(bkgY);
		      if (significance>fSigLowThreshold && significance<fSigHighThreshold && fFit->GetParameter(4)>fMinWidth)
		      {
			      X[i][numPeak] = fFit->GetParameter(3);
			      Y[i][numPeak] = AmBePeaks[t];
			      eX[i][numPeak] = fFit->GetParError(3);
			      eY[i][numPeak] = 0.;
			      numPeak++;
		      }
		  }

		  fh_Map_fit[i+2*fNumCrystals]->Write();

		  auto graph1 = new TGraphErrors(numPeak,X[i],Y[i],eX[i],eY[i]);
		  graph1->Fit("f1","WQN");
		  graph1->Fit("f1","MQ+");

		  chi2 = TMath::Log10(FindChisquare(X[i],Y[i],eX[i],numPeak,f1));
		}

		if (chi2>fChi2Threshold)
		{
			f1->SetParameter(0,0.);
			f1->SetParameter(1,0.);
			deleted_peaks++;
		}

		if (i < fNumCrystals/2 && f1->GetParameter(1)>0.) // reset problematic crystal
		{
		    if (f1->GetParameter(1)<fMinSlope or f1->GetParameter(1)>fMaxSlope)
		    {
			f1->SetParameter(0,0.);
			f1->SetParameter(1,0.);
			deleted_peaks++;
		    }
		} else if (i >= fNumCrystals/2 && f1->GetParameter(1)>0.)
		{
		    if (f1->GetParameter(1)<fMinSlopeP or f1->GetParameter(1)>fMaxSlopeP)
		    {
			f1->SetParameter(0,0.);
			f1->SetParameter(1,0.);
			deleted_peaks++;
		    }
		}

		// write histograms
		TPolyMarker *pm1 = new TPolyMarker(numPeak,X[i],Y[i]);
		pm1->SetMarkerStyle(23);
		pm1->SetMarkerColor(kBlack);
		fh2_peak_cal[i]->GetListOfFunctions()->Add(pm1);
		for (Int_t peakIdx=0; peakIdx<numPeak; peakIdx++)
		{
		    TLine *line = new TLine(X[i][peakIdx]-eX[i][peakIdx],Y[i][peakIdx],X[i][peakIdx]+eX[i][peakIdx],Y[i][peakIdx]);
		    fh2_peak_cal[i]->GetListOfFunctions()->Add(line);
		}
		fh2_peak_cal[i]->GetListOfFunctions()->Add(f1);
		fh2_peak_cal[i]->SetStats(0);
		fh2_peak_cal[i]->Write();

		// Fill slope histogram
		Double_t pm2X[1] = {(Double_t) crystalID};
		Double_t pm2Y[1] = {f1->GetParameter(1)};
		TPolyMarker *pm2 = new TPolyMarker(1,pm2X,pm2Y);
		pm2->SetMarkerStyle(23);
		pm2->SetMarkerColor(kBlack);
		fh2_slope_crystalID->GetListOfFunctions()->Add(pm2);

		Double_t pm3X[1] = {(Double_t) crystalID};
		Double_t pm3Y[1] = {chi2};
		TPolyMarker *pm3 = new TPolyMarker(1,pm3X,pm3Y);
		pm3->SetMarkerStyle(23);
		pm3->SetMarkerColor(kBlack);
		fh2_chi2_crystal->GetListOfFunctions()->Add(pm3);

		Double_t pm4X[1] = {(Double_t) crystalID};
		Double_t pm4Y[1] = {f1->GetParameter(0)};
		TPolyMarker *pm4 = new TPolyMarker(1,pm4X,pm4Y);
		pm4->SetMarkerStyle(23);
		pm4->SetMarkerColor(kBlack);
		fh2_intercept_crystalID->GetListOfFunctions()->Add(pm4);

		for (Int_t peakIdx=0; peakIdx<numPeak; peakIdx++)
		{
		    Double_t pm5X[1] = {Y[i][peakIdx]};
		    Double_t pm5Y[1] = {f1->Eval(X[i][peakIdx])-Y[i][peakIdx]};
		    TPolyMarker *pm5 = new TPolyMarker(1,pm5X,pm5Y);
		    pm5->SetMarkerStyle(23);
		    pm5->SetMarkerColor(kBlack);
		    fh2_residual_energy[i]->GetListOfFunctions()->Add(pm5);
		}
		TLine *zeroLine = new TLine(0,0,5000,0);
		zeroLine->SetLineColor(kRed);
		fh2_residual_energy[i]->GetListOfFunctions()->Add(zeroLine);
		fh2_residual_energy[i]->SetStats(0);
		fh2_residual_energy[i]->Write();

		for (Int_t h=0; h<fNumParam; h++)
		{
		    fCal_Par->SetCryCalParams(f1->GetParameter(h), fNumParam*i+h);
		}
	    }

            fh_numPeak->Fill(numPeak);
	}
    }

    std::cout << "deleted peaks: " << deleted_peaks << std::endl;

    fh2_slope_crystalID->Write();
    fh2_intercept_crystalID->Write();
    fh2_sig_crystal[0]->Write();
    fh2_sig_crystal[1]->Write();
    fh2_sig_crystal[2]->Write();
    fh2_chi2_crystal->Write();
    fh2_resolution_crystalID->Write();
    fh_numPeak->Write();

    fCal_Par->setChanged();
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
