#include <TString.h>

TString sourcename;
TString filename;
TString SpectrumFileName;


//___________________________________________________________________change file names and select parameters/functions___________________________________________________________________//

TString short_filename_gamma = "gamma_pulser_comb";       //   gamma_onlysource_comb     gamma_pulser_comb 
TString folder_name_gamma = "/home/e12exp/data_calibration/September_data/11_9/";     // September_data/11_9

TString short_filename_proton = "proton_pulser_comb";       //  proton_onlysource_comb    proton_pulser_comb
TString folder_name_proton = "/home/e12exp/data_calibration/September_data/11_9/";     // September_data/11_9

// input root file with collected data
TString filename_gamma = folder_name_gamma + "rootfiles/" + short_filename_gamma + ".root";
TString filename_proton = folder_name_proton + "rootfiles/" + short_filename_proton + ".root";

// output spectrum file with filled histograms
TString SpectrumFileName_gamma = folder_name_gamma + "spectra/" + short_filename_gamma + "_spectrum.root";
TString SpectrumFileName_proton = folder_name_proton + "spectra/" + short_filename_proton + "_spectrum.root";

// output file with calibrated histograms
TString outputFileName = folder_name_proton + "rootfiles/" + short_filename_proton + "_calibrated.root";

// CALIFA output file with the parameters calibrated in keV
TString outputCalFile = folder_name_proton + "CalPar/" + short_filename_proton + ".par";

// Parameters for CALIFA
TString califamapfilename = "/home/e12exp/data_calibration/param_files/califamapping_v1.par";

// all errors from calibration: peak number, ...
TString PeakErrors = "/home/e12exp/data_calibration/September_data/peak_errors.txt";

// source_energies in keV
Float_t source_energies[] = {511.0, 1274.5};

// Pulser voltages in mV
Float_t voltages_gamma_values[] = {45, 270, 495};  
Float_t voltages_proton_values[] = {79, 7800};

const Int_t nev = -1; // number of events to read, -1 - until CTRL+C
Int_t fRunId = 213;
TString cRunId = Form("%04d", fRunId);


//___________________________________________________________________run califa_calibParFinder_v3_____________________________________________________________________________________//

#include "califa_calibParFinder_v3.C"

void pulser_calibration_v1() 
{
    
    sourcename = "spectrum";
    
    
    // create gamma_spectrum.root
    filename = filename_gamma;
    SpectrumFileName = SpectrumFileName_gamma;    
    califa_calibParFinder_v3();
    
    
    /*
    // create proton_spectrum.root
    filename = filename_proton;
    SpectrumFileName = SpectrumFileName_proton;
    califa_calibParFinder_v3();
    */
    
    /*
    // calibration 
    sourcename = "pulser";
    califa_calibParFinder_v3();
    */
    
    gApplication->Terminate();
}

