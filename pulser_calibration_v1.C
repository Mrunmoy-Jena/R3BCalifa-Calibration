//____________________________________Instructions__________________________________________//
/*

Order:
bash_pulser_calibration.sh (bash script) -> pulser_calibration_v1.C (parameters, ...) 
-> califa_calibParFinder_v3.C (calibration settings) -> R3BCalifaMapped2CrystalCalPar.cxx (calibration macro) -> ...

You need:
-measurement_gamma_range.lmd files with source and pulser events 
-measurement_proton_range.lmd files with source and pulser (same as in gamma range + more with higher energy) events 

for both ranges each, combine all files and unpack it to .root file 

1) if needed change fit parameters, ... in califa_calibParFinder_v3.C 
 
2) change file and folder names, source energies and pulser voltages in pulser_calibration_v1.C

3) configurate r3b: .../build/r3broot_build/config.sh

4) build r3b: cmake --build .../build/r3broot_build

5) make bash script executable: chmod +x bash_pulser_calibration.sh
 
6) run bash_pulser_calibration.sh: ./bash_pulser_calibration.sh

*/
//________________________________________________________________________________________//

#include <TString.h>

TString sourcename;
TString filename;
TString SpectrumFileName;


//_________________change file names and select parameters/functions_____________________//

TString short_filename_gamma = "gamma_pulser_comb";       //   gamma_onlysource_comb     gamma_pulser_comb 
TString folder_name_gamma = "/home/e12exp/data_calibration/September_data/11_9/";

TString short_filename_proton = "proton_pulser_comb";       //  proton_onlysource_comb    proton_pulser_comb
TString folder_name_proton = "/home/e12exp/data_calibration/September_data/11_9/";

TString short_filename = "pulser_comb";
TString folder_name = "/home/e12exp/data_calibration/September_data/11_9/";


// input root file with collected data
TString filename_gamma = folder_name_gamma + "rootfiles/" + short_filename_gamma + ".root";
TString filename_proton = folder_name_proton + "rootfiles/" + short_filename_proton + ".root";

// output spectrum file with filled histograms
TString SpectrumFileName_gamma = folder_name_gamma + "spectra/" + short_filename_gamma + "_spectrum.root";
TString SpectrumFileName_proton = folder_name_proton + "spectra/" + short_filename_proton + "_spectrum.root";

// output file with calibrated histograms
TString outputFileName = folder_name + "rootfiles/" + short_filename + "_calibrated.root";

// CALIFA output file with the parameters calibrated in keV
TString outputCalFile = folder_name + "CalPar/" + short_filename + ".par";

// all errors from calibration: peak number, ...
TString PeakErrors = folder_name + "peak_errors.txt";

// Parameters for CALIFA
TString califamapfilename = "/home/e12exp/data_calibration/param_files/califamapping_v1.par";

// CrystalID <-> angles
TString anglesfilename = "/home/e12exp/crystalID_angles.txt";

// source_energies in keV
Double_t source_energies[] = {1274.5};

// offset: if 0, then offset as a variable
Double_t offset_calibration = -20;

// Pulser voltages in mV
Double_t voltages_gamma_values[] = {45, 252, 425};  
Double_t voltages_proton_values[] = {45, 252, 425, 1300, 3010, 5620};

const Int_t nev = -1; // number of events to read, -1 - until CTRL+C
Int_t fRunId = 213;
TString cRunId = Form("%04d", fRunId);


//________________________run califa_calibParFinder_v3________________________________//

#include "califa_calibParFinder_v3.C"

void pulser_calibration_v1(TString function) 
{
    
    if (function == "spectrum_gamma")
    {
        // create gamma_spectrum.root
        sourcename = "spectrum";
        filename = filename_gamma;
        SpectrumFileName = SpectrumFileName_gamma;    
        califa_calibParFinder_v3();   
    }    
    else if (function == "spectrum_proton")
    {
        // create proton_spectrum.root
        sourcename = "spectrum";
        filename = filename_proton;
        SpectrumFileName = SpectrumFileName_proton;
        califa_calibParFinder_v3();
    }
    else if (function == "pulser_calibration")
    {
        // calibration 
        sourcename = "pulser";
        califa_calibParFinder_v3();
    }
    
    gApplication->Terminate();
}


