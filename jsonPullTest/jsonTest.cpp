#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>

int main() {
    // Read JSON file
    std::ifstream ifs("settings.json");
    if (!ifs.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    // Parse JSON
    nlohmann::json j;
    try {
        ifs >> j;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        return 1;
    }

    // Accessing values
    float gratingPeriod = j["GRATING_PERIOD"];
    float wavelength = j["WAVELENGTH"];
    int numberOfPeriods = j["NUMBER_OF_PERIODS"];
    int zSteps = j["ZSTEPS"];
    int xySteps = j["XYSTEPS"];

    // Example usage
    std::cout << "Grating Period: " << gratingPeriod << std::endl;
    std::cout << "Wavelength: " << wavelength << std::endl;
    std::cout << "Number of Periods: " << numberOfPeriods << std::endl;
    std::cout << "Z Steps: " << zSteps << std::endl;
    std::cout << "XY Steps: " << xySteps << std::endl;

    return 0;
}