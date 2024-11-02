
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"

using namespace std;

// Struktur für die CSV- und ROOT-Daten
struct CrystalData {
    int CrystalID = 0;
    int HALF = 0;
    int RING = 0;
    int PREAMP = 0;
    int CHANNEL = 0;
    int FEBEX_PC = 0;
    int FEBEX_SFP = 0;
    int FEBEX_CH = 0;
    int FEBEX_MOD = 0; 
    double AngleTheta = 0.0;
    double AnglePhi = 0.0;
    string Region = "Unknown";
};

// Hilfsfunktion für sichere Konvertierung
int safeStoi(const string &value, int defaultValue = 0) {
    try {
        return stoi(value);
    } catch (...) {
        return defaultValue;
    }
}

// Funktion zum Einlesen der CSV-Dateien und Zuordnen von `CRYSTALID` als `CrystalID`
void readCSV(const string &filename, map<int, CrystalData> &dataMap) {
    cout << "Lese CSV-Datei: " << filename << endl;
    ifstream file(filename);
    string line;

    if (!file.is_open()) {
        cerr << "Fehler beim Öffnen der Datei: " << filename << endl;
        return;
    }

    getline(file, line);  // Überspringe die Header-Zeile

    int linkedEntries = 0; // Zählt verknüpfte Einträge aus der CSV
    int outputCount = 0;   // Zählt die Anzahl der ausgegebenen Einträge für Terminalprüfung

    while (getline(file, line)) {
        // Überspringe leere Zeilen (z.B. Block-Trennungen)
        if (line.empty()) continue;

        stringstream ss(line);
        CrystalData entry;
        string value;

        // Lese die relevanten Spalten
        getline(ss, value, ','); entry.HALF = safeStoi(value);        // HALF
        getline(ss, value, ','); entry.RING = safeStoi(value);        // RING
        getline(ss, value, ','); entry.PREAMP = safeStoi(value);      // PREAMP
        getline(ss, value, ','); entry.CHANNEL = safeStoi(value);     // CHANNEL
        
        // Überspringe die nicht relevanten Spalten: `CRYSTAL TYPE`, `APD NUMBER`, `VOLTAGE`
        for (int i = 0; i < 3; ++i) getline(ss, value, ',');

        // Lese die restlichen relevanten Spalten
        getline(ss, value, ','); entry.FEBEX_PC = safeStoi(value);    // FEBEX PC
        getline(ss, value, ','); entry.FEBEX_SFP = safeStoi(value);   // FEBEX SFP
getline(ss, value, ','); entry.FEBEX_MOD = safeStoi(value);           // FEBEX_MOD
        getline(ss, value, ','); entry.FEBEX_CH = safeStoi(value);    // FEBEX CH
        getline(ss, value, ','); // Lese `CRYSTALID`
        entry.CrystalID = safeStoi(value);

        // Überspringe Zeilen, bei denen `CRYSTALID` leer oder null ist, ohne Warnung
        if (entry.CrystalID == 0) continue;

        // Wenn CrystalID bereits in der Map ist, fügen wir die CSV-Daten hinzu
        auto it = dataMap.find(entry.CrystalID);
        if (it != dataMap.end()) {
            it->second.HALF = entry.HALF;
            it->second.RING = entry.RING;
            it->second.PREAMP = entry.PREAMP;
            it->second.CHANNEL = entry.CHANNEL;
            it->second.FEBEX_PC = entry.FEBEX_PC;
            it->second.FEBEX_SFP = entry.FEBEX_SFP;
            it->second.FEBEX_MOD = entry.FEBEX_MOD;
            it->second.FEBEX_CH = entry.FEBEX_CH;
            linkedEntries++; // Ein erfolgreicher Link wird gezählt

            // Ausgabe der ersten 10 erfolgreichen Verknüpfungen zur Prüfung
            if (outputCount < 10) {
                cout << "CSV-Verknüpfung - CrystalID: " << entry.CrystalID
                     << ", HALF: " << entry.HALF << ", RING: " << entry.RING
                     << ", PREAMP: " << entry.PREAMP << ", CHANNEL: " << entry.CHANNEL
                     << ", FEBEX_PC: " << entry.FEBEX_PC << ", FEBEX_SFP: " << entry.FEBEX_SFP
                     << ", FEBEX_MOD: " << entry.FEBEX_MOD << ", FEBEX_CH: " << entry.FEBEX_CH 
                     << endl;
                outputCount++;
            }
        } else {
            cout << "Warnung: CrystalID " << entry.CrystalID << " in CSV-Datei, aber nicht in ROOT-Daten vorhanden." << endl;
        }
    }

    file.close();
    cout << "CSV-Datei " << filename << " erfolgreich eingelesen. Verknüpfte Einträge: " << linkedEntries << endl;
}

// Funktion zum Einlesen der Winkel aus den ROOT-Dateien
void readRootData(const string &filename, map<int, CrystalData> &dataMap, int minID, int maxID) {
    const Long64_t maxEvents = 5000000;
    cout << "Lese ROOT-Datei: " << filename << " für CrystalID-Bereich " << minID << " bis " << maxID << endl;
    TFile *file = TFile::Open(filename.c_str());
    if (!file || file->IsZombie()) {
        cerr << "Fehler beim Öffnen der ROOT-Datei: " << filename << endl;
        return;
    }

    TTree *tree = (TTree*)file->Get("evt");
    if (!tree) {
        cerr << "Fehler: TTree 'evt' in der Datei " << filename << " nicht gefunden!" << endl;
        file->Close();
        return;
    }

    UShort_t fCrystalId;
    Double_t fTheta, fPhi;

    tree->SetBranchAddress("CalifaCrystalCalData.fCrystalId", &fCrystalId);
    tree->SetBranchAddress("CalifaClusterData.fTheta", &fTheta);
    tree->SetBranchAddress("CalifaClusterData.fPhi", &fPhi);

    Long64_t nEntries = tree->GetEntries();
    Long64_t nProcessed = min(nEntries, maxEvents);
    cout << "Anzahl der Einträge in ROOT-Datei: " << nEntries << " (Maximal zu verarbeiten: " << nProcessed << ")" << endl;

    for (Long64_t i = 0; i < nProcessed; ++i) {
        tree->GetEntry(i);

        if (i % 1000000 == 0) {
            cout << "Verarbeitete Einträge: " << i << "/" << nProcessed << endl;
        }

        if (fCrystalId >= minID && fCrystalId <= maxID) {
            double thetaDeg = fTheta * 180.0 / M_PI;
            double phiDeg = fPhi * 180.0 / M_PI;
            string region;
            if (thetaDeg >= 43.2 && thetaDeg <= 140.3) {
                region = "Barrel";
            } else if (thetaDeg >= 19.3 && thetaDeg < 43.2) {
                region = "Iphos";
            } else if (thetaDeg >= 7.2 && thetaDeg < 19.3) {
                region = "CEPA";
            } else {
                region = "Unknown";
            }

            // Füge Winkelinformationen in die Map ein
            dataMap[fCrystalId] = {fCrystalId, 0, 0, 0, 0, 0, 0, 0, 0, thetaDeg, phiDeg, region};
        }
    }

    file->Close();
    cout << "ROOT-Datei " << filename << " erfolgreich eingelesen." << endl;
}



void all_to_crystalID() {
    int fNumCrystals = 5088; // Gesamtanzahl der CrystalIDs
    map<int, CrystalData> crystalDataMap;

    // Schritt 1: Einlesen der ROOT-Dateien (Winkel)
    readRootData("/home/e12exp/data_calibration/September_data/gamma_252mV_comb_small_Cal2Clus.root", crystalDataMap, 0, 2544);
    readRootData("/home/e12exp/data_calibration/September_data/proton_3V_comb_small_Cal2Clus.root", crystalDataMap, 2545, 5088);

    // Schritt 2: Einlesen der CSV-Dateien (weitere Daten)
    readCSV("/home/e12exp/califa_gamma.csv", crystalDataMap);
    readCSV("/home/e12exp/califa_proton.csv", crystalDataMap);

    // Schritt 3: Ausgabe der kombinierten Daten
    ofstream outFile("califa_specs.txt");
    if (!outFile.is_open()) {
        cerr << "Fehler beim Öffnen der Ausgabedatei!" << endl;
        return;
    }

    cout << "Schreibe Daten in die Ausgabedatei 'califa_specs.txt'..." << endl;
    outFile << "CrystalID Theta(Deg) Phi(Deg) Area HALF RING PREAMP CHANNEL FEBEX_PC FEBEX_SFP FEBEX_MOD FEBEX_CH\n";

    // Iteriere über alle CrystalIDs von 1 bis fNumCrystals und schreibe die Daten
    for (int i = 1; i <= fNumCrystals; ++i) {
        auto it = crystalDataMap.find(i);
        if (it != crystalDataMap.end()) {
            // Wenn CrystalID gefunden, schreibe die gespeicherten Werte
            const CrystalData &data = it->second;
            outFile << data.CrystalID << " "
                    << data.AngleTheta << " "
                    << data.AnglePhi << " "
                    << data.Region << " "
                    << data.HALF << " "
                    << data.RING << " "
                    << data.PREAMP << " "
                    << data.CHANNEL << " "
                    << data.FEBEX_PC << " "
                    << data.FEBEX_SFP << " "
                    << data.FEBEX_MOD << " " 
                    << data.FEBEX_CH << "\n";
        } else {
            // Wenn CrystalID nicht gefunden, schreibe 0er und "Unknown"
            outFile << i << " 0 0 Unknown 0 0 0 0 0 0 0 0\n";
        }
    }

    outFile.close();
    cout << "Die Datei 'califa_specs.txt' wurde erfolgreich erstellt." << endl;
}
