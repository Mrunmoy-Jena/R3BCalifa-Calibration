typedef struct EXT_STR_h101_t {
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_CALIFA_t califa;
} EXT_STR_h101;

#include <filesystem>


//_______________________________________confirmation_______________________________________________//
bool confirmExecution() 
{
    string input;
    cout << "A spectrum file already exists. Do you want to overwrite these? (yes/no),  ";
    cin >> input;

    // Convert input to lowercase for case-insensitive comparison
    for (auto& c : input) {
        c = tolower(c);
    }

    if (input == "yes" || input == "y") {
        return true;
    } else {
        cout << "The already existing spectrum is used." << endl;
        return false;
    }
}


//________________________________________HistoFill________________________________________________//
void createSpectrum() 
{
    // Create run
    FairRunAna *run = new FairRunAna();

    // Set up R3BHeader 
    R3BEventHeader *EvntHeader = new R3BEventHeader();
    run->SetEventHeader(EvntHeader);
    run->SetSink(new FairRootFileSink(outputFileName));

    // Runtime data base 
    FairRuntimeDb *rtdb = run->GetRuntimeDb();
    R3BFileSource *source = new R3BFileSource(filename);
    run->SetSource(source);

    // CALIFA mapping
    FairParAsciiFileIo *parIo1 = new FairParAsciiFileIo();
    parIo1->open(califamapfilename, "in");
    rtdb->setFirstInput(parIo1);
    rtdb->addRun(fRunId);
    rtdb->initContainers(fRunId);
    rtdb->print();

    R3BEventHeaderPropagator *RunIdTask = new R3BEventHeaderPropagator();
    run->AddTask(RunIdTask);
                    
    // R3BCalifaMapped2CrystalCalPar
    R3BCalifaMapped2CrystalCalPar *CalPar = new R3BCalifaMapped2CrystalCalPar();
    cout << "after creating calpar pointer" << endl;

    CalPar->SetSpectrumName(SpectrumFileName);
    CalPar->SetSourceName(sourcename);

    CalPar->SetMinStadistics(200);
    CalPar->SetNumParameterFit(2);

    // Gamma range
    CalPar->SetCalRange_left(0);
    CalPar->SetCalRange_right(30000);
    CalPar->SetCalRange_bins(3000);

    // particle range
    CalPar->SetCalRangeP_left(0);
    CalPar->SetCalRangeP_right(30000);
    CalPar->SetCalRangeP_bins(30000);

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
    CalPar->SetDebugMode(0);

    cout << "before adding the task" << endl;
    run->AddTask(CalPar);
    cout << "before initialization" << endl;

    // Initialize -------------------------------------------
    run->Init();
    FairLogger::GetLogger()->SetLogScreenLevel("INFO");   // INFO, WARNING, DEBUG

    //output Ascii file with the Calibration Parameters
    FairParAsciiFileIo *parOut = new FairParAsciiFileIo();
    parOut->open(outputCalFile, "out");
    rtdb->setOutput(parOut);

    cout << "after setting parOut" << endl;

    if (nev > -1) {
      run->Run(nev);
    }
    else
    {
      run->Run();
    }

    rtdb->saveOutput();

    cout << endl << endl;
    cout << "Spectrum created succesfully." << endl;
    cout << "spectrum file is " << SpectrumFileName << endl;

}


//______________________________________PulserCalibration__________________________________________________//
void pulserCalibration() 
{     

    // R3BCalifaMapped2CrystalCalPar
    R3BCalifaMapped2CrystalCalPar *CalPar = new R3BCalifaMapped2CrystalCalPar();
    cout << "after creating calpar pointer" << endl;   
 
    TArrayF* EnergythePeaks = new TArrayF(sizeof(source_energies) / sizeof(source_energies[0]));
    TArrayF* PulserVoltages_gamma = new TArrayF(sizeof(voltages_gamma_values) / sizeof(voltages_gamma_values[0]));
    TArrayF* PulserVoltages_proton = new TArrayF(sizeof(voltages_proton_values) / sizeof(voltages_proton_values[0]));

    for (int i = 0; i < EnergythePeaks->GetSize(); ++i) {
        EnergythePeaks->AddAt(source_energies[i], i);
    }

    for (int i = 0; i < PulserVoltages_gamma->GetSize(); ++i) {
        PulserVoltages_gamma->AddAt(voltages_gamma_values[i], i);
    }

    for (int i = 0; i < PulserVoltages_proton->GetSize(); ++i) {
        PulserVoltages_proton->AddAt(voltages_proton_values[i], i);
    }

    CalPar->SetEnergyPeaks(EnergythePeaks);
    cout << "after setting energy peaks" << endl;
  
    CalPar->SetPulserVoltages_gamma(PulserVoltages_gamma);
    cout << "after setting pulser voltages_gamma" << endl;
    
    CalPar->SetPulserVoltages_proton(PulserVoltages_proton);
    cout << "after setting pulser voltages_proton" << endl;    
    
    //______________more_parameters__________________//
    CalPar->SetSpectrumName_gamma(SpectrumFileName_gamma);
    CalPar->SetSpectrumName_proton(SpectrumFileName_proton);
    CalPar->SetoutputName(outputFileName);
    CalPar->SetCalName(outputCalFile);
    CalPar->SetRunId(fRunId);
    CalPar->Setcalifamapfilename(califamapfilename);
    CalPar->SetSourceName(sourcename);
    CalPar->SetPeakErrors(PeakErrors);

    CalPar->SetDebugMode(0);
    CalPar->SetMinStadistics(200);
    CalPar->SetNumParameterFit(2);
    
    // Gamma range
    CalPar->SetCalRange_left(0);
    CalPar->SetCalRange_right(30000);
    CalPar->SetCalRange_bins(3000);
    
    // particle range
    CalPar->SetCalRangeP_left(0);
    CalPar->SetCalRangeP_right(30000);
    CalPar->SetCalRangeP_bins(30000);
    
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
    
    //call PulserCalibration function
    CalPar->PulserCalibration();


    cout << endl << endl;
    cout << "Calibration finished succesfully." << endl;
    cout << "Output file is " << outputFileName << endl;
    cout << ".par file is " << outputCalFile << endl;

}


//_______________________________________califa_calibParFinder_v3_____________________________________________________//
void califa_calibParFinder_v3()
{
    TStopwatch timer;
    timer.Start();

    califamapfilename.ReplaceAll("//", "/");

    if (sourcename == "spectrum") 
    {
        if ((!std::filesystem::exists(SpectrumFileName.Data())) || confirmExecution())
        {
            createSpectrum(); 
        }
    }
    else if (sourcename == "pulser")
    {
        pulserCalibration(); 
    }
    else
    {
        cout << "You have selected a non-existent detector range: '" << sourcename << "', Change range!" << endl;
    }
    
    
    // Finish -----------------------------------------------
    timer.Stop();
    Double_t rtime = timer.RealTime() / 60.;
    Double_t ctime = timer.CpuTime() / 60.;
    cout << "Real time " << rtime << " min, CPU time " << ctime << " min"
         << endl
         << endl;

}
