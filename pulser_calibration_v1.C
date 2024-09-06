#include <TString.h>

TString sourcename;
TString filename;
TString SpectrumFileName;


//___________________________________________________________________change file names and select parameters/functions___________________________________________________________________//

TString short_filename_gamma = "analog_gamma_07_07_all";       // analog_gamma_07_07_all   
TString folder_name_gamma = "/home/e12exp/data_calibration/Juli_data/7_7/";     // 7_7

TString short_filename_proton = "analog_proton_10x_wo780mV";       // analog_proton_10x_05_07_all   analog_proton_79mV_05_07_comb   analog_proton_780mV_05_07_comb    analog_proton_7_8V_05_07_comb   analog_proton_10x_wo780mV
TString folder_name_proton = "/home/e12exp/data_calibration/Juli_data/6_7/";     // 6_7

// input root file with collected data
TString filename_gamma = folder_name_gamma + "rootfiles/" + short_filename_gamma + ".root";
TString filename_proton = folder_name_proton + "rootfiles/" + short_filename_proton + ".root";

// output spectrum file with filled histograms
TString SpectrumFileName_gamma = folder_name_gamma + "spectra/" + short_filename_gamma + "_spectrum.root";
TString SpectrumFileName_proton = folder_name_proton + "spectra/" + short_filename_proton + "_spectrum.root";

// output file with calibrated histograms
TString outputFileName = folder_name_gamma + "rootfiles/" + short_filename_gamma + "_calibrated.root";

// CALIFA output file with the parameters calibrated in keV
TString outputCalFile = folder_name_gamma + "CalPar/" + short_filename_gamma + ".par";

// Parameters for CALIFA
TString califamapfilename = "/home/e12exp/data_calibration/param_files/califamapping_v1.par";

// all errors from calibration: peak number, ...
TString PeakErrors = "/home/e12exp/data_calibration/Juli_data/peak_errors.txt";

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

    // create proton_spectrum.root
    filename = filename_proton;
    SpectrumFileName = SpectrumFileName_proton;
    califa_calibParFinder_v3();
    

    // calibration 
    sourcename = "pulser";
    califa_calibParFinder_v3();
    
    gApplication->Terminate();
}


