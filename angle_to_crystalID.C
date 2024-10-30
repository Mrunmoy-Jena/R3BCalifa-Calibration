#include <iostream>
#include <fstream>
#include <map>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

//_________________________how_to__________________________//
/*
root
.x angle_to_crystalID.C

->crystalID_angles.txt
*/
//________________________________________________________//

// Event structure for ROOT tree
struct EventData {
    UShort_t fCrystalId;
    Double_t fTheta;
    Double_t fPhi;
};

// Function to determine the region based on Theta in degrees
std::string getRegionName(double thetaDeg) {
    if (thetaDeg >= 43.2 && thetaDeg <= 140.3) {
        return "Barrel";
    } else if (thetaDeg >= 19.3 && thetaDeg < 43.2) {
        return "Iphos";
    } else if (thetaDeg >= 7.2 && thetaDeg < 19.3) {
        return "CEPA";
    } else {
        return "Unknown";  // If Theta is outside known regions
    }
}

void angle_to_crystalID() {
    // Open file and tree
    TFile *file = TFile::Open("/home/e12exp/data_calibration/September_data/proton_3V_comb_small_Cal2Clus.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }
    
    TTree *tree = (TTree*)file->Get("evt");
    if (!tree) {
        std::cerr << "Error: Tree 'evt' not found!" << std::endl;
        file->Close();
        return;
    }
    
    // Structures to store data
    EventData califaCrystalCalData;
    EventData califaClusterData;
    
    // Set branches
    tree->SetBranchAddress("CalifaCrystalCalData.fCrystalId", &califaCrystalCalData.fCrystalId);
    tree->SetBranchAddress("CalifaClusterData.fTheta", &califaClusterData.fTheta);
    tree->SetBranchAddress("CalifaClusterData.fPhi", &califaClusterData.fPhi);
    
    // Maximum number of CrystalIDs and events
    const int maxCrystalId = 5088;
    const Long64_t maxEvents = 9000000; // Process up to 9,000,000 events
    
    // Map to store found CrystalIDs and their Theta, Phi, and region
    std::map<int, std::tuple<double, double, std::string>> crystalAngles;
    
    // Iterate over all events with progress display
    Long64_t nEntries = tree->GetEntries();
    Long64_t nProcessed = std::min(nEntries, maxEvents);

    for (Long64_t i = 0; i < nProcessed; ++i) {
        tree->GetEntry(i);
        
        // Progress display: show message every 1,000,000 events
        if (i % 1000000 == 0) {
            std::cout << "Processed events: " << i << " of " << nProcessed << std::endl;
        }

        // Check if this CrystalID has already been stored
        if (crystalAngles.find(califaCrystalCalData.fCrystalId) == crystalAngles.end()) {
            // Store CrystalID and associated angles (converted to degrees)
            double thetaDeg = califaClusterData.fTheta * 180.0 / M_PI;
            double phiDeg = califaClusterData.fPhi * 180.0 / M_PI;
            std::string region = getRegionName(thetaDeg);  // Determine region
            crystalAngles[califaCrystalCalData.fCrystalId] = std::make_tuple(thetaDeg, phiDeg, region);
        }
        
        // Break if all CrystalIDs have been found
        if (crystalAngles.size() >= maxCrystalId) {
            break;
        }
    }
    
    // Create and open the output file
    std::ofstream outFile("crystalID_angles.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file!" << std::endl;
        file->Close();
        return;
    }
    
    // Write header in the first line
    outFile << "CrystalID Theta(Deg) Phi(Deg) Area\n";
    
    // Write all CrystalIDs to the file
    for (int id = 1; id <= maxCrystalId; ++id) {
        if (crystalAngles.find(id) != crystalAngles.end()) {
            double theta = std::get<0>(crystalAngles[id]);
            double phi = std::get<1>(crystalAngles[id]);
            std::string region = std::get<2>(crystalAngles[id]);
            outFile << id << " " << theta << " " << phi << " " << region << "\n";
        } else {
            outFile << id << " 0 0 Unknown\n";  // If no angle is found, set dummy values (0)
        }
    }
    
    // Close the file
    outFile.close();
    
    std::cout << "The file 'crystalID_angles.txt' was successfully created." << std::endl;
    
    // Close ROOT file
    file->Close();
}

