/*paCh2 Additional info:
 *
 *
 * One needs to set up the 2024 experiments, the unpacker is:
 *
 * at $UCESB_DIR/../upexps/202402_s091_s118
 *
 * Before executing the macro, set up the paths to lmd files or to stream server,
 * upexps and califa mapping file "califamapfilename". Also set up the name of file
 * "califacalfilename" that contains the energy calibration parameters of CALIFA.
 *
 * @since January 18th, 2024
 * */

typedef struct EXT_STR_h101_t
{
    EXT_STR_h101_unpack_t unpack;
    EXT_STR_h101_TPAT_t tpat;
    EXT_STR_h101_CALIFA_t califa;
    EXT_STR_h101_WRMASTER_t wrm;
    EXT_STR_h101_WRCALIFA_t wrcalifa;
} EXT_STR_h101;

void califa_offline_unpacker(int nev = -1)
{
    TStopwatch timer;

    FairLogger::GetLogger()->SetLogScreenLevel("info");
    FairLogger::GetLogger()->SetColoredLog(true);

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    std::ostringstream oss;
    oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
    TString time_suffix = oss.str();

    const int expId = 118; // select experiment: 118 

    // Create input -----------------------------------------
    // TString inputFileName="main0070_all.lmd";
    //TString inputFileName="--stream=lxlanddaq01:9000";
    TString inputFile;   //use wildcard
    TString in_dir = "/home/e12exp/ssd/data_calibration/lmd/"; //the directory of your lmd file
    TString inputFileName = in_dir+"analog_proton_760mV*.lmd";  //number of input lmd files (with wildcard *)
    // Output file ------------------------------------------
    TString out_dir = "/home/e12exp/ssd/data_calibration/unpacked/"; //directory where you want to store the unpacked file
    TString outputFileName = out_dir+"analog_proton_760mV_all.root";

    bool Cal_level = true;          // set true if there exists a file with the calibration parameters
    bool NOTstoremappeddata = false; // if true, don't store mapped data in the root file
    bool NOTstorecaldata = false;    // if true, don't store cal data in the root file
    bool NOTstorehitdata = false;    // if true, don't store hit data in the root file

    // Online server configuration --------------------------
    int refresh = 10; // Refresh rate for online histograms
    int port = 8867;

    // Calibration files ------------------------------------
    //TString dir = "/media/mrunmoy/MyDisk/calibration_macros"; // gSystem->Getenv("VMCWORKDIR");
    // CALIFA detector
    // Parameters for CALIFA
    TString currentDir = gSystem->Getenv("PWD");
    TString califamapfilename = currentDir+"/califamapping_v1.par";    
    califamapfilename.ReplaceAll("//", "/");

    TString califacalfilename = currentDir+"/califacalpar_v3.par";
    califacalfilename.ReplaceAll("//", "/");

    // UCESB configuration ----------------------------------
    TString ntuple_options = "RAW";
    TString ucesb_dir = getenv("UCESB_DIR");
    TString ucesb_path = ucesb_dir + "/../upexps/202402_s091_s118/202402_s091_s118 --allow-errors --input-buffer=200Mi";
    //TString ucesb_path = ucesb_dir + "/../upexps/202205_s522/202205_s522 --allow-errors --input-buffer=200Mi";
    //TString ucesb_path = ucesb_dir + "/../upexps/202205_s509/202205_s509 --allow-errors --input-buffer=180Mi";
    ucesb_path.ReplaceAll("//", "/");

    // Create online run ------------------------------------
    auto *run = new FairRunOnline();
    auto *EvntHeader = new R3BEventHeader();
    EvntHeader->SetExpId(expId);
    run->SetEventHeader(EvntHeader);
    run->SetRunId(1);
    run->SetSink(new FairRootFileSink(outputFileName));
    run->ActivateHttpServer(refresh, port);

    // Create source using ucesb for input ------------------
    EXT_STR_h101 ucesb_struct;
    R3BUcesbSource* source =
	    new R3BUcesbSource(inputFileName, ntuple_options, ucesb_path, &ucesb_struct, sizeof(ucesb_struct));
    source->SetMaxEvents(nev);

    // Definition of reader ---------------------------------
    // Add readers ------------------------------------------
    source->AddReader(new R3BUnpackReader((EXT_STR_h101_unpack*)&ucesb_struct.unpack, offsetof(EXT_STR_h101, unpack)));
    
    auto unpacktrl = new R3BTrloiiTpatReader(&ucesb_struct.tpat, offsetof(EXT_STR_h101, tpat));
    // unpacktrl->SetTrigger(1);
    // unpacktrl->SetTpatRange(1,12);
    source->AddReader(unpacktrl);

    source->AddReader(new R3BWhiterabbitMasterReader((EXT_STR_h101_WRMASTER *)&ucesb_struct.wrm,
                                                          offsetof(EXT_STR_h101, wrm), 0x1000));

    auto unpackcalifa = new R3BCalifaFebexReader((EXT_STR_h101_CALIFA*)&ucesb_struct.califa,
                                                                   offsetof(EXT_STR_h101, califa));
    unpackcalifa->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackcalifa);

    auto unpackWRC = new R3BWhiterabbitCalifaReader((EXT_STR_h101_WRCALIFA*)&ucesb_struct.wrcalifa,
			                            offsetof(EXT_STR_h101, wrcalifa), 0x1a00, 0x1b00);
    unpackWRC->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackWRC);

    run->SetSource(source);

    // Runtime data base ------------------------------------
    auto rtdb = run->GetRuntimeDb();



    // Load parameters --------------------------------------
    if (!Cal_level)
    {
	    // CALIFA mapping
	    auto parIo1 = new FairParAsciiFileIo(); // Ascii file
	    parIo1->open(califamapfilename, "in");
	    rtdb->setFirstInput(parIo1);
	    rtdb->print();
    } 
    else 
    {
	    // CALIFA mapping and calibration in keV (both parameters should be here)
	    auto parIo2 = new FairParAsciiFileIo(); // Ascii file
	    parIo2->open(califacalfilename, "in");
	    rtdb->setFirstInput(parIo2);
	    rtdb->print();
    }

    // Create analysis task ---------------------------------
    if (Cal_level)
    {
	    // R3BCalifaMapped2CrystalCal ---
	    auto Map2Cal = new R3BCalifaMapped2CrystalCal();
	    Map2Cal->SetOnline(NOTstorecaldata);
	    run->AddTask(Map2Cal);

	    //TJ commented out, not needed for now
//	    // R3BCalifaCrystalCal2Cluster ---
//	    auto Cal2Clus = new R3BCalifaCrystalCal2Cluster();
//	    Cal2Clus->SetCrystalThreshold(100.);//100keV
//	    Cal2Clus->SetProtonClusterThreshold(20000.);//12MeV
//	    Cal2Clus->SetGammaClusterThreshold(2000); //4MeV
//	    Cal2Clus->SetRoundWindow(0.3);
//	    Cal2Clus->SelectGeometryVersion(2024);
//            Cal2Clus->SetRandomization(kTRUE);
//            Cal2Clus->SetRandomizationFile("/u/land/r3broot/202402_s091_s118/R3BParams_s091_s118/califa/angular_distributions.root");
//	    Cal2Clus->SetOnline(NOTstorehitdata);
//	    run->AddTask(Cal2Clus);
    }

    	//TJ commented out, not needed for now
//    // R3BOnlineSpectra -------------------------------------
//    auto CalifaOnline = new R3BCalifaOnlineSpectra();
//    //CalifaOnline->SetRange_max(200000);//10MeV original *****
//    CalifaOnline->SetRange_max(360000);//10MeV original *****
//    //CalifaOnline->SetRange_max(10000);//10MeV
//    CalifaOnline->SetRange_bins(10000);
//    CalifaOnline->SetBinChannelFebex(500);
//    CalifaOnline->SetMaxBinFebex(30000);
//    CalifaOnline->SetTotHist(true);
//    //CalifaOnline->SetTpat(1,10);shold(0.0)

//    run->AddTask(CalifaOnline);

                                                                                                                                 
    //TJ commented out, not needed for now
//    auto genonline = new R3BGeneralOnlineSpectra();                                                                                
//    run->AddTask(genonline); 


    // Initialize -------------------------------------------
    timer.Start();
    run->Init();

    // Informations about portnumber and main data stream.
    cout<<"\n\n"<<endl;
    cout<<"Data stream is: "<<inputFileName<<endl;
    cout<<"Portnumber for califa online is: "<<port<<endl;
    cout<<"\n\n"<<endl;
    // Run --------------------------------------------------
    cout << "Just before run ..." << endl;
    run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);
	
    cout << "after execution of run->Run.." << endl;
    // -----   Finish   -------------------------------------
    cout << endl << endl;
    // Extract the maximal used memory an add is as Dart measurement
    // This line is filtered by CTest and the value send to CDash
    FairSystemInfo sysInfo;
    Float_t maxMemory = sysInfo.GetMaxMemory();
    cout << "MaxMemory: ";
    cout << maxMemory << endl;

    timer.Stop();
    Double_t rtime = timer.RealTime()/60.;
    Double_t ctime = timer.CpuTime()/60.;
    Float_t cpuUsage = ctime / rtime;
    cout << "CPU used: " << cpuUsage << endl;

    cout << endl;
    std::cout << "Output file is " << outputFileName << std::endl;
    cout << "Real time " << rtime << " min, CPU time " << ctime << " min" << endl << endl;
    cout << "Macro finished successfully." << endl;
    gApplication->Terminate();
}
