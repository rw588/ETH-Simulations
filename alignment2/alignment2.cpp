#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

struct Parameters {
    float gratingPeriod;
    float wavelength;
    int numberOfPeriods;
    int zSteps;
    int xySteps;
};

Parameters readParametersFromJson();

int main() {
    Parameters params = readParametersFromJson();

    // Access parameters
    float gratingPeriod = params.gratingPeriod;
    float wavelength = params.wavelength;
    int numberOfPeriods = params.numberOfPeriods;
    int zSteps = params.zSteps;
    int xySteps = params.xySteps;

    // Print parameters
    std::cout << "Grating Period: " << params.gratingPeriod << std::endl;
    std::cout << "Wavelength: " << params.wavelength << std::endl;
    std::cout << "Number of Periods: " << params.numberOfPeriods << std::endl;
    std::cout << "Z Steps: " << params.zSteps << std::endl;
    std::cout << "XY Steps: " << params.xySteps << std::endl;

    return 0;
}

Parameters readParametersFromJson() {
    Parameters params;

    // Read JSON file
    std::ifstream ifs("settings.json");
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return params;
    }

    // Parse JSON
    json j;
    try {
        ifs >> j;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        return params;
    }

    // Accessing values
    try {
        params.gratingPeriod = j["GRATING_PERIOD"];
        params.wavelength = j["WAVELENGTH"];
        params.numberOfPeriods = j["NUMBER_OF_PERIODS"];
        params.zSteps = j["ZSTEPS"];
        params.xySteps = j["XYSTEPS"];
    } catch (const std::exception& e) {
        std::cerr << "Error accessing JSON values: " << e.what() << std::endl;
    }

    return params;
}
