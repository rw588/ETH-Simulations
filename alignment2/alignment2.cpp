#include <TMath.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph2D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TVirtualFFT.h>
#include <iostream>
#include <TRandom3.h>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct Parameters {
    float gratingPeriod;
    float wavelength;
    int numberOfPeriods;
    int zSteps;
    int xySteps;
};

Parameters json();

float json()

int main(){
    Parameters params = json();
    float gratingPeriod;
    float wavelength;
    int numberOfPeriods;
    int zSteps;
    int xySteps;

    return 0;

}

Parameters json() {
    Parameters params;

    // Read JSON file
    std::ifstream ifs("settings.json");
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return params;
    }

    // Parse JSON
    nlohmann::json j;
    try {
        ifs >> j;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        return params;
    }

    // Accessing values
    float gratingPeriod = j["GRATING_PERIOD"];
    float wavelength = j["WAVELENGTH"];
    int numberOfPeriods = j["NUMBER_OF_PERIODS"];
    int zSteps = j["ZSTEPS"];
    int xySteps = j["XYSTEPS"];

    return params;
}

