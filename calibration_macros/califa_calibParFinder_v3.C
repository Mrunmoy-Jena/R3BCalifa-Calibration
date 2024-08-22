#include<boost/filesystem.hpp>
#include"TSystem.h"

typedef struct EXT_STR_h101_t {
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_CALIFA_t califa;
} EXT_STR_h101;


bool confirmExecution(const string& action) {
    string input;
    cout << "You have selected \"" << action << "\". This may overwrite existing data. Do you want to continue? (yes/no): ";
    cin >> input;

    // Convert input to lowercase for case-insensitive comparison
    for (auto& c : input) {
        c = tolower(c);
    }

    if (input == "yes" || input == "y") {
        return true;
    } else {
        cout << "Execution cancelled by user." << endl;
        return false;
    }
}


void califa_calibParFinder_v3()
{
  TStopwatch timer;
  timer.Start();

  const Int_t nev = -1; // number of events to read, -1 - until CTRL+C
  Int_t fRunId = 213;
  TString cRunId = Form("%04d", fRunId);
  
  
  //_____________________________choose function____________________________________________//
  
  TString filename, outputFileName, SpectrumFileName, outputCalFile, califamapfilename;
  TString ynchoice, sourcename;
  cout<<"Do you have a spectrum.root file ? (file with energy spectra histograms for each crystal)"<<endl;
  cout<<"Enter yes/no:";
  cin>>ynchoice;
  
  Bool_t flag1; //flag=true for a valid choice, false for incorrect choice
  do{
  	if(ynchoice=="y"||ynchoice=="yes"){	//if spectrum file exists
  	flag1=true;
  	cout<<"Enter name of spectrum.root file (with full path):";
  	Bool_t flag2;
	  	do{
	  	cin>>SpectrumFileName;
	  	ifstream file(SpectrumFileName);
	  	string specfile = string(SpectrumFileName);
	  		if(file && specfile.substr(specfile.length() - 5)==".root"){
	  		flag2=true;
	  		}
	  		else{
	  		flag2=false;
	  		cout<<"File does not exist ! Please re-enter the correct filename:";
	  		}
	  	}while(!flag2);
	sourcename="22Na_pulser";  	
  	}
  	
  	else if(ynchoice=="n"||ynchoice=="no"){	//if spectrum file doesn't exist
  	flag1=true; 
  	cout<<"Enter name of unpacked input .root file (with full path):";
  	Bool_t flaginp;
  		do{
  		cin>>filename;
  		ifstream fileinp(filename);
  		string strfile = string(filename);
  			if(fileinp && strfile.substr(strfile.length() - 5)==".root"){	//input .root file exists
  			flaginp=true;
  			}
  			else{			//input .root doesn't exist
  			flaginp=false;
  			cout<<"File does not exist ! Please re-enter the correct filename:";
  			}
  		}while(!flaginp);
  		
  	cout<<"Enter name for output spectrum .root file (with full path):";
  	Bool_t flagspect;
  		do{
    		cin>>SpectrumFileName;
    		string strspect = string(SpectrumFileName);  
    		size_t found = strspect.find_last_of("/");
    		string path = strspect.substr(0, found);    
	    		if(strspect.substr(strspect.length() - 5)!=".root"){
	    		flagspect=false;
	    		cout<<"Output spectrum file must end with '.root' !"<<endl;
	    		cout<<"please re-enter correct filename:";
	    		}
	    		else if(!gSystem->cd(path.c_str())){
	    		flagspect=false;
	    		cout<<"Specified directory doesn't exist !"<<endl;
	    		cout<<"Please re-enter correct filepath:";
	    		}
			else{
			flagspect=true;
			}
    		}while(!flagspect);	
    	sourcename = "spectrum";	//source is set to spectrum, will now create spectrum.root file
  	}
  	else{
	cout<<"Invalid choice !";
	cout<<"\n"<<"Please re-enter your choice:";
	cin>>ynchoice;
	}
  }while(!flag1);
  
  //cout<<"Sourcename:"<<sourcename<<endl;
  
  //Mapping parameters for CALIFA
  Bool_t flag3;
  cout<<"Choose input mapping parameters file (with full path):";
  do{
  	cin>>califamapfilename;
	ifstream file(califamapfilename);
	string calfile = string(califamapfilename);
		if(file && calfile.substr(calfile.length() - 4)==".par"){
	  	flag3=true;
	  	}
	  	else{
	  	flag3=false;
	  	cout<<"File does not exist ! Please re-enter the correct filename:";
	  	}
   }while(!flag3);
   califamapfilename.ReplaceAll("//", "/");  

  //run HistoFill first if you don't already have the specific spectrum.root file for your data!!!!!!
  
  //HistoFill: Fill Histogramms, create spectrum.root (spectrum)
  //FitPeaks: ... (fitting)
  //PulserCalibration: ... (22Na_pulser, 60Co_pulser, AmBe_pulser, 152Eu_pulser)
  

  // confirm execution for "fitting" and "spectrum"   --> Deactivated by MJ for now
  /*
  if (sourcename == "spectrum" || sourcename == "fitting") 
  {
     
     if (!confirmExecution(string(sourcename.Data()))) 
     {
         return;
     }
     
  }
  */
  
  //________________________Create source using ucesb for input___________________________//
  /*
  TString name_short = "analog_gamma_07_07_all";      //board_gamma_07_07_all  analog_gamma_07_07_all  analog_gamma_45mV_07_07_comb

  TString filename = TString("/home/e12exp/ssd/data_calibration/unpacked/7_7_24/") + name_short + ".root";
  TString outputFileName = TString("/home/e12exp/ssd/data_calibration/calibrated/") + name_short + "_calibrated.root";
  cout << "after in/out files" << endl;

  //output spectrum file with filled histograms
  TString SpectrumFileName = TString("/home/e12exp/ssd/data_calibration/spectra/") + name_short + "_spectrum.root";
  cout << "after spectrum files" << endl;
  
  // CALIFA output file with the parameters calibrated in keV
  TString outputCalFile = TString("/home/e12exp/ssd/data_calibration/CalPar/") + name_short + "_calibrated.par";
  cout << "output califa cal file" << endl;
  
  //Parameters for CALIFA
  TString califamapfilename = "/home/e12exp/ssd/R3BRoot/califa/calibration/calibration_macros/califamapping_v1.par";
  califamapfilename.ReplaceAll("//", "/");
  */
  
  /*
  //Input file
  cout<<"Choose input file name (with full path):";  
  cin>>filename;
  cout<<"\n"<<"Input filename entered"<<endl;
  
  //Mapping parameters for CALIFA
  cout<<"Choose input mapping parameters file (with full path):";
  cin>>califamapfilename;
  califamapfilename.ReplaceAll("//", "/");  
  cout<<"\n"<<"Mapping .par filename entered"<<endl;
  */
  
  //________________________________________HistoFill, FitPeaks________________________________________________//

  if (sourcename == "spectrum" /*|| sourcename == "fitting" */) 
  {         
          
          // Create run  --------------------------------------------------------------
	  FairRunAna *run = new FairRunAna();
	  //cout << "before eventheader" << endl;
	  
	  // Set up R3BHeader  --------------------------------------------------------
	  R3BEventHeader *EvntHeader = new R3BEventHeader();
	  run->SetEventHeader(EvntHeader);
	  //run->SetSink(new FairRootFileSink(outputFileName));   //commented out by MJ (need not write to calibrated.root)
	  cout << "\n"<<"after eventheader" << endl;

	  // Runtime data base ------------------------------------
	  FairRuntimeDb *rtdb = run->GetRuntimeDb();
	  R3BFileSource *source = new R3BFileSource(filename);
	  // source->AddFile("/igfae_lustre/NUCL/s509/rootfiles/s509_0213_002x_map_20230429_175930.root");
	  run->SetSource(source);

	  // CALIFA mapping
	  FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo(); // Ascii
	  parIo1->open(califamapfilename, "in");
	  rtdb->setFirstInput(parIo1);
	  rtdb->addRun(fRunId);
	  rtdb->initContainers(fRunId);
	  rtdb->print();
          
	  R3BEventHeaderPropagator *RunIdTask = new R3BEventHeaderPropagator();
	  run->AddTask(RunIdTask);
                        
	  // R3BCalifaMapped2CrystalCalPar ----
	  cout << "before setting calib peak parameters..." << endl;
	  R3BCalifaMapped2CrystalCalPar *CalPar = new R3BCalifaMapped2CrystalCalPar();
	  cout << "after creating calpar pointer" << endl;
	  
	  CalPar->SetSpectrumName(SpectrumFileName);
	  CalPar->SetSourceName(sourcename);
	  
	  CalPar->SetMinStadistics(10);
	  CalPar->SetNumParameterFit(2); // OPTIONAL by default 2
	  
	  // Gamma range
	  CalPar->SetCalRange_left(300); // Na: 700 Co: 650 AmBe: 2200 Eu: 950 fitting: 0
	  CalPar->SetCalRange_right(30000); // Na: 1700 Co: 1700 AmBe: 5000 Eu: 1700 fitting: 5000
	  CalPar->SetCalRange_bins(2700); // Na: 140 Co: 105 AmBe: 200 Eu: 75 fitting: 250
	  cout << "before particle range" << endl;
	  
	  // particle range
	  CalPar->SetCalRangeP_left(30); // Na: 80 Co: 70 AmBe: 210 Eu: 90 fit: 0
	  CalPar->SetCalRangeP_right(30000); // Na: 120 Co: 130 AmBe: 360 Eu: 150 fit: 400
	  CalPar->SetCalRangeP_bins(2700); // Na: 40 Co: 60 AmBe: 140 Eu: 60 fit:200
	  
	  CalPar->SetSigma(2.0); //3.0
	  CalPar->SetThreshold(0.05);
	  CalPar->SetMaxPeaks(20);
	  CalPar->SetChi2Threshold(-3);
	  CalPar->SetSigLowThreshold(3.0);
	  CalPar->SetSigHighThreshold(500.);
	  CalPar->SetMinWidth(1.0);
	  CalPar->SetGausRange(30.0);
	  CalPar->SetGausRangeP(5.0);
	  CalPar->SetGausBaseEnergy(511.0);
	  CalPar->SetMinSlope(1.0);
	  CalPar->SetMaxSlope(1.6);
	  CalPar->SetMinSlopeP(10.0);
	  CalPar->SetMaxSlopeP(16.0);	  
	  CalPar->SetDebugMode(0);

	  cout << "before adding the task" << endl;
	  run->AddTask(CalPar);
	  cout << "before initialization" << endl;
	  
	  // Initialize -------------------------------------------
	  run->Init();
	  //FairLogger::GetLogger()->SetLogScreenLevel("WARNING");
	  //FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
	  FairLogger::GetLogger()->SetLogScreenLevel("INFO");
	  
	  //output Ascii file with the Calibration Parameters         //commented out by MJ (no need to write output .par file here)
	  //FairParAsciiFileIo *parOut = new FairParAsciiFileIo();
	  //parOut->open(outputCalFile, "out");
	  //rtdb->setOutput(parOut);

	  //cout << "after setting parOut" << endl;
	  
      if (nev > -1) {
          run->Run(nev);
      }
      else
      {
          run->Run();
      }
      
      //rtdb->saveOutput();
      
  TString nextchoice; 
  Bool_t flagok;   
  cout<<"Completed writing to spectrum file:"<<SpectrumFileName;
  cout<<"Do you want to continue with the calibration procedure? : (yes/no)";
  do{
	  cin>>nextchoice;
	  if(nextchoice=="yes" || "y"){
	  flagok=true;
	  sourcename = "22Na_pulser";
	  }    
	  else if(nextchoice=="no" || "n"){
	  flagok=true;
	  cout<<"You have chosen not to continue with calibration.";
	  }
	  else{
	  flagok=false;
	  cout<<"Invalid choice entered ! Please re-enter your choice:";
	  }
  }while(!flagok);
  
  }


  //______________________________________PulserCalibration__________________________________________________//


  if (sourcename == "22Na_pulser" /*|| sourcename == "60Co_pulser" || sourcename == "AmBe_pulser"  || sourcename == "152Eu_pulser" */ )
  {
                
    // CALIFA output file with the parameters calibrated in keV
    Bool_t flag4;
    cout<<"\n"<<"Enter name for output parameters file (with full path, with extension as .par):";
    do{
    	cin>>outputCalFile;
    	string stroutCalFile = string(outputCalFile);  
    	size_t found = stroutCalFile.find_last_of("/");
    	string path = stroutCalFile.substr(0, found);    
    		if(stroutCalFile.substr(stroutCalFile.length() - 4)!=".par"){
    		flag4=false;
    		cout<<"Output parameter file must end with '.par' !"<<endl;
    		cout<<"please re-enter correct filename:";
    		}
    		else if(!gSystem->cd(path.c_str())){
    		flag4=false;
    		cout<<"Specified directory doesn't exist !"<<endl;
    		cout<<"Please re-enter correct filepath:";
    		}
		else{
		flag4=true;
		}
    }while(!flag4);
    
    //Calibrated .root file
    Bool_t flag5;
    cout<<"\n"<<"Enter name for output calibrated root file (with full path, with extension as .root):";
    do{
    	cin>>outputFileName;
    	string stroutrootfile = string(outputFileName);  
    	size_t found = stroutrootfile.find_last_of("/");
    	string path = stroutrootfile.substr(0, found);    
    		if(stroutrootfile.substr(stroutrootfile.length() - 5)!=".root"){
    		flag5=false;
    		cout<<"Output root file must end with '.root' !"<<endl;
    		cout<<"please re-enter correct filename:";
    		}
    		else if(!gSystem->cd(path.c_str())){
    		flag5=false;
    		cout<<"Specified directory doesn't exist !"<<endl;
    		cout<<"Please re-enter correct filepath:";
    		}
		else{
		flag5=true;
		}
    }while(!flag5);
    
    
    // R3BCalifaMapped2CrystalCalPar ----
    cout << "before setting calib peak parameters..." << endl;
    R3BCalifaMapped2CrystalCalPar *CalPar = new R3BCalifaMapped2CrystalCalPar();
    cout << "after creating calpar pointer" << endl;   
  
    TArrayF* EnergythePeaks = new TArrayF();
    
    // e1 must be always the lowest energy
    Float_t e1 = 511.0;         // 511.0 1173.2 3927
    Float_t e2 = 1274.5;        // 1274.5 1332.5 4438 1408.0
    EnergythePeaks->Set(2);     //2 or more Peaks
    EnergythePeaks->AddAt(e1, 0);
    EnergythePeaks->AddAt(e2, 1);

    CalPar->SetEnergyPeaks(EnergythePeaks);
    cout <<" after setting energy peaks" << endl;
   
    CalPar->SetSpectrumName(SpectrumFileName);
    CalPar->SetoutputName(outputFileName);
    CalPar->SetCalName(outputCalFile);
    CalPar->SetRunId(fRunId);
    CalPar->Setcalifamapfilename(califamapfilename);
    CalPar->SetSourceName(sourcename);
    
    CalPar->SetDebugMode(0);
    CalPar->SetMinStadistics(10);
    CalPar->SetNumParameterFit(2);
    
    // Gamma range
    CalPar->SetCalRange_left(300);
    CalPar->SetCalRange_right(30000);
    CalPar->SetCalRange_bins(2700);
    
    // particle range
    CalPar->SetCalRangeP_left(30);
    CalPar->SetCalRangeP_right(30000);
    CalPar->SetCalRangeP_bins(2700);
    
    CalPar->SetSigma(2.0);
    CalPar->SetThreshold(0.05);
    CalPar->SetMaxPeaks(20);
    CalPar->SetChi2Threshold(-3);
    CalPar->SetSigLowThreshold(3.0);
    CalPar->SetSigHighThreshold(500.);
    CalPar->SetMinWidth(1.0);
    CalPar->SetGausRange(30.0);
    CalPar->SetGausRangeP(5.0);
    CalPar->SetGausBaseEnergy(511.0);
    CalPar->SetMinSlope(1.0);
    CalPar->SetMaxSlope(1.6);
    CalPar->SetMinSlopeP(10.0);
    CalPar->SetMaxSlopeP(16.0);
    CalPar->SetMaxSigma(50.0);
    CalPar->SetMinPeakEvents(100.0);
    
    //number of pulser signals, the energies must be higher than the source energies!!
    CalPar->SetPulserNumber(3);
    
    //call PulserCalibration function
    cout<<"Calling pulser calibration function";
    CalPar->PulserCalibration();
    cout << "Output file is " << outputFileName << endl;
  }



  //________________________________________________________________________________________________________//


  /*
  else
  {
    cout << "You have selected a non-existent function: '" << sourcename << "', Change sourcename!" << endl;
  }
  */

  // Finish -----------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime() / 60.;
  Double_t ctime = timer.CpuTime() / 60.;
  cout << endl << endl;
  cout << "Macro finished succesfully." << endl;
  cout << "Real time " << rtime << " min, CPU time " << ctime << " min"
            << endl
            << endl;
  gApplication->Terminate();
}
