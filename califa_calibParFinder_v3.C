typedef struct EXT_STR_h101_t {
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_CALIFA_t califa;
} EXT_STR_h101;

void califa_calibParFinder_v3() {
  TStopwatch timer;
  timer.Start();

  const Int_t nev = -1; // number of events to read, -1 - until CTRL+C
  Int_t fRunId = 213;
  TString cRunId = Form("%04d", fRunId);
  // Create source using ucesb for input ------------------
  TString filename = "/home/e12exp/data_calibration/Juli_data/7_7/rootfiles/all_analog_gamma_07_07.root";
  TString outputFileName = "/home/e12exp/data_calibration/Juli_data/7_7/rootfiles/all_analog_gamma_07_07_calibrated.root";
  cout << "after in/out files" << endl;

  // Parameters for CALIFA
  TString dir = "/home/e12exp/data_calibration/param_files"; // gSystem->Getenv("VMCWORKDIR");
  TString califamapdir = dir;
  TString califamapfilename = dir +"/califamapping_v1.par";

  califamapfilename.ReplaceAll("//", "/");

  // CALIFA output file with the parameters calibrated in keV
  TString outputCalFile = "/home/e12exp/data_calibration/Juli_data/7_7/CalPar/all_analog_gamma_07_07_calibrated.par";
  cout << "output califa cal file" << endl;

  // Create source using ucesb for input ------------------
  // EXT_STR_h101 ucesb_struct;

  // Create run  --------------------------------------------------------------
  FairRunAna *run = new FairRunAna();
  cout << "before eventheader" << endl;
  
  // Set up R3BHeader  --------------------------------------------------------
  R3BEventHeader *EvntHeader = new R3BEventHeader();
  run->SetEventHeader(EvntHeader);
  run->SetSink(new FairRootFileSink(outputFileName));
  cout << "after eventheader" << endl;

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

  // TFile *file = TFile::Open(filename);
  // TTree *tree = (TTree *)file->Get("evt");
  // TClonesArray *califaData = new TClonesArray("R3BCalifaMappedData");
  // TBranch *branchCalifaData = tree->GetBranch("CalifaMappedData");
  // branchCalifaData->SetAddress(&califaData);
  // Add analysis task ------------------------------------
  
  R3BEventHeaderPropagator *RunIdTask = new R3BEventHeaderPropagator();
  run->AddTask(RunIdTask);
  

   // ---- calibration ---- //
   
   // R3BCalifaMapped2CrystalCalPar ----
   cout << "before setting calib peak parameters..." << endl;
   R3BCalifaMapped2CrystalCalPar *CalPar = new R3BCalifaMapped2CrystalCalPar();
   cout << "after creating calpar pointer" << endl;
  
   TArrayF *EnergythePeaks = new TArrayF();
   // e1 must be always the lowest energy
   Float_t e1 = 511.0; // 511.0 1173.2 3927
   Float_t e2 = 1274.5; // 1274.5 1332.5 4438 1408.0
   EnergythePeaks->Set(2); //2 or more Peaks
   EnergythePeaks->AddAt(e1, 0);
   EnergythePeaks->AddAt(e2, 1);
  
  
  //_____________________________choose function____________________________________________//
  
  //Search Peak: find peaks and draw histograms (22Na, 60Co, AmBe, 152Eu)
  //Fit Peak: ... (fitting)
  //Pulser Calibration: ... (22Na_pulser, 60Co_pulser, AmBe_pulser, 152Eu_pulser)
  
  CalPar->SetSourceName("22Na_pulser");
  
  //________________________________________________________________________________________//
  
  
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
  CalPar->SetMaxSigma(50.0);
  CalPar->SetMinPeakEvents(100.0);
  
  //number of pulser signals, the energies must be higher than the source energies!!
  CalPar->SetPulserNumber(3);
  
  
  cout <<" before setting energy peaks" << endl;
  CalPar->SetEnergyPeaks(EnergythePeaks);
  cout <<" after setting energy peaks" << endl;
  
  CalPar->SetDebugMode(0);
  cout << "after debug mode" << endl;
  
  cout << "before adding the task" << endl;
  run->AddTask(CalPar);
  cout << "before initialization" << endl;
  
  // Initialize -------------------------------------------
  run->Init();
  //FairLogger::GetLogger()->SetLogScreenLevel("WARNING");
  //FairLogger::GetLogger()->SetLogScreenLevel("DEBUG");
  FairLogger::GetLogger()->SetLogScreenLevel("INFO");

  // Choose Root or Ascii file
  // 1-Root file with the Calibration Parameters
  /*Bool_t kParameterMerged = kTRUE;
  FairParRootFileIo *parOut = new FairParRootFileIo(kParameterMerged);
  parOut->open(outputCalFile);
  rtdb->setOutput(parOut);*/
  
  // 2-Ascii file with the Calibration Parameters // the one used
  FairParAsciiFileIo *parOut = new FairParAsciiFileIo();
  parOut->open(outputCalFile, "out");
  rtdb->setOutput(parOut);

  cout << "after setting parOut" << endl;

  // Run --------------------------------------------------
  if (nev > -1)
    run->Run(nev);
  else
    run->Run();
  // Save parameters (if needed) --------------------------
  rtdb->saveOutput();

  // Finish -----------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime() / 60.;
  Double_t ctime = timer.CpuTime() / 60.;
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << "Output file is " << outputFileName << std::endl;
  std::cout << "Real time " << rtime << " min, CPU time " << ctime << " min"
            << std::endl
            << std::endl;
  gApplication->Terminate();
}
