#include <iostream>
#include <fstream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include <math.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TChain.h"
#include <string>
#include <vector>



void analysis1000(){
//ofstream myfile("./pu_c_data/real_data_co60.txt");
//ofstream energyfile("./pu_c_data/clustered_energies.txt");
ofstream myfile("./pu_c_data/real_data_puc.txt");
ofstream energyfile("./pu_c_data/clustered_energies.txt");
char fname[200];
char hist_name[500];
int good_counts;
TH1F* h1_gamma_spec;
sprintf(hist_name, "Gamma spec");
h1_gamma_spec = new TH1F(hist_name,hist_name,1000,0,10);
h1_gamma_spec->GetXaxis()->SetTitle("Uncalibrated Energy, mapped level");
h1_gamma_spec->GetYaxis()->SetTitle("Counts");
h1_gamma_spec->GetXaxis()->CenterTitle(true);
h1_gamma_spec->GetYaxis()->CenterTitle(true);
h1_gamma_spec->GetYaxis()->SetLabelSize(0.045);
h1_gamma_spec->GetYaxis()->SetTitleSize(0.045);

TH1F* h1_califa_crystal_mapped_mult;
sprintf(hist_name, "CALIFA MAPPED MULTIPLICITY");
h1_califa_crystal_mapped_mult = new TH1F(hist_name,hist_name,100,0,100);
h1_califa_crystal_mapped_mult->GetXaxis()->SetTitle("CALIFA MULTIPLICITY, mapped level");
h1_califa_crystal_mapped_mult->GetYaxis()->SetTitle("Counts");
h1_califa_crystal_mapped_mult->GetXaxis()->CenterTitle(true);
h1_califa_crystal_mapped_mult->GetYaxis()->CenterTitle(true);
h1_califa_crystal_mapped_mult->GetYaxis()->SetLabelSize(0.045);
h1_califa_crystal_mapped_mult->GetYaxis()->SetTitleSize(0.045);
  
TChain* chain = new TChain("evt");

//sprintf(fname, "/home/e12exp/tobias_jenegger_python/new_cluster_025_main0146_0001_stitched.root");
//sprintf(fname, "/home/e12exp/tobias_jenegger_python/pu_c_data/no_cluster_main0244_0157_stitched.root");
//sprintf(fname, "/home/e12exp/tobias_jenegger_python/pu_c_data/new_cluster_025_main0244_0157_stitched.root");
sprintf(fname, "/home/e12exp/tobias_jenegger_python/no_cluster_main0146_0001_stitched.root");

	chain->Reset();
	chain->Add(fname);

      TClonesArray*  MappedCA = new TClonesArray("R3BCalifaMappedData",5);
      R3BCalifaMappedData** mappedCA;
      TBranch *branchMappedCA = chain->GetBranch("CalifaMappedData");
      branchMappedCA->SetAddress(&MappedCA);
      
      TClonesArray*  CalibCA = new TClonesArray("R3BCalifaCrystalCalData",5);
      R3BCalifaCrystalCalData** calibratedCA;
      TBranch *branchcalibratedCA = chain->GetBranch("CalifaCrystalCalData");
      branchcalibratedCA->SetAddress(&CalibCA);
      
      TClonesArray*  HitCA = new TClonesArray("R3BCalifaClusterData",5);
      R3BCalifaClusterData** hitCA;
      TBranch *branchhitCA = chain->GetBranch("CalifaClusterData");
      branchhitCA->SetAddress(&HitCA);

      Long64_t MappedCAPerEvent = 0;
  
  	Long64_t nevent = chain->GetEntries();

	
std::vector<int> R4PA5 = {3391,3392,3387,3388,3390,3389,3386,3385,3519,3520,3515,3516,3518,3517,3514,3513};

std::vector<int> R3PA9 = {3647,3648,3643,3644,3646,3645,3642,3641,3775,3776,3771,3772,3774,3773,3770,3769}; 

	cout << "Events = " << nevent << endl;
int numbers_to_analyse = 0;
for(Long64_t i = 0; i < nevent; i++  ){
//for(Long64_t i = 0; i < 10000; i++  ){
    if(i%100000==0 && i>0)cout<<"Processing event " <<i<< " => " << (static_cast<double>(i)/static_cast<double>(nevent))*100. << " %" << endl;
	MappedCA->Clear();
	chain->GetEvent(i);
	MappedCAPerEvent = MappedCA->GetEntriesFast();
	int64_t evtnr = i;
	//cout << "MappedCAPerEvent:\t"  << MappedCA->GetEntriesFast() <<   "and calib mult per event:\t" << CalibCA->GetEntriesFast() << "with eventnr. \t" << evtnr << endl;
	
       if(MappedCAPerEvent>0){
	  h1_califa_crystal_mapped_mult->Fill(MappedCAPerEvent);
	  mappedCA = new R3BCalifaMappedData*[MappedCAPerEvent];
	  for(int p = 0; p < MappedCAPerEvent; p++){
			mappedCA[p] =(R3BCalifaMappedData*)MappedCA->At(p);
			UShort_t CrystalId = (mappedCA[p]->GetCrystalId());
			Int_t Energy = (mappedCA[p]->GetEnergy());
			//h1_gamma_spec->Fill(double(Energy)/1000.);	
			//what I have to do:
			//make a vector of vectors having as entries: entrynumber / energy [MeV] / polar angle [degrees] /azimuthal angle [degrees] / time [ns]
			//then sort the vectors of vectors according to the time vec[0].time < vec[1].time
			//then substract from all vectors of vectors the time vev[0].time
			//write it as above in comma separated file
			//ah, and you should be in the cal level

			
            Int_t Nf = (mappedCA[p]->GetNf());
            Int_t Ns = (mappedCA[p]->GetNs());
	    //cout << "this is crystalid:\t" << CrystalId << endl;
	 
	  }
 	
      }
      //write the hit  part:
        vector< vector<std::variant<double,ULong64_t>> > event_vec;
	if (HitCA->GetEntriesFast()){
		numbers_to_analyse +=1;
		if (numbers_to_analyse > 100000) break;
		hitCA = new R3BCalifaClusterData*[(HitCA->GetEntriesFast())];
		for (int i = 0; i < HitCA->GetEntriesFast(); i++){
			std::vector<std::variant<double,ULong64_t>> hit_vec;
			hitCA[i] = (R3BCalifaClusterData*)HitCA->At(i);
			Double_t energy = (hitCA[i]->GetEnergy())/1000.;
			hit_vec.push_back(energy);
			Double_t theta = ((hitCA[i]->GetTheta())/3.14159265359)*180.;
			hit_vec.push_back(theta);
			Double_t phi = ((hitCA[i]->GetPhi())/3.14159265359)*180.;
			hit_vec.push_back(phi);
			ULong64_t time = (hitCA[i]->GetTime());
			hit_vec.push_back(time);
			event_vec.push_back(hit_vec);
			h1_gamma_spec->Fill(energy);
			energyfile << energy << endl;

			if (energy > 1. && energy < 1.4) good_counts += 1;

			cout << "this should be the energy:\t" << std::get<double>(hit_vec[0]) << endl;
		}
		delete [] hitCA;
		std::sort(event_vec.begin(), event_vec.end(),[](std::vector<std::variant<double,ULong64_t>>& a, std::vector<std::variant<double,ULong64_t>>& b) {return a[3] < b[3];});
		if (event_vec.size() ){
			cout << "begin of vector sorting...." << endl;	
			ULong64_t low_val = std::get<ULong64_t>(event_vec[0][3]);
			for (int i = 0; i < event_vec.size() ; i++){
				cout << "energy:\t" << std::get<double>(event_vec[i][0]) << "theta:\t" << std::get<double>(event_vec[i][1]) << "phi:\t" << std::get<double>(event_vec[i][2]) << "time:\t" << std::get<ULong64_t>(event_vec[i][3]) << endl;
				ULong64_t now_val = std::get<ULong64_t>(event_vec[i][3]);
				event_vec[i][3]= now_val - low_val + 100;
				}
			for (int i = 0; i < event_vec.size() ; i++){

				cout << "energy:\t" << std::get<double>(event_vec[i][0]) << "theta:\t" << std::get<double>(event_vec[i][1]) << "phi:\t" << std::get<double>(event_vec[i][2]) << "time:\t" << std::get<ULong64_t>(event_vec[i][3]) << endl;
				//make a vector of vectors having as entries: entrynumber / energy [MeV] / polar angle [degrees] /azimuthal angle [degrees] / time [ns]
				myfile << evtnr << "," << std::get<double>(event_vec[i][0]) << "," << std::get<double>(event_vec[i][1]) << "," << std::get<double>(event_vec[i][2]) << "," << std::get<ULong64_t>(event_vec[i][3]) << endl;		
				}

			cout << "END of vector sorting...." << endl;
		}
	}

      //if (MappedCAPerEvent>2) cout << "this is event with many entries:\t" << evtnr << endl;
   }
h1_gamma_spec->Draw();
cout << "Entries in Histogram:\t" << h1_gamma_spec->GetEntries() << endl;
cout << "good_counts:\t" << good_counts << endl;
cout << "percentage number:\t" << double(good_counts)/double(h1_gamma_spec->GetEntries())*100 << "\%" << endl;
myfile.close();
energyfile.close();
//h1_califa_crystal_mapped_mult->Draw();
} 
