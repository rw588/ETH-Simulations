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
#include <omp.h>

// Define constants
const Int_t NUMBER_OF_PERIODS = 20;
const Int_t ZSTEPS = 3000;
const Int_t XYSTEPS = 300;
const Int_t RESOLUTION = ZSTEPS/NUMBER_OF_PERIODS;

// Forward declaration of functions
Double_t Pattern(TRandom3* randomGen, Double_t x, Double_t y, Double_t z, Double_t GRATING_PERIOD, Double_t TALBOT_LENGTH, Double_t Z0, Double_t L);
void plane_to_grid(TVectorD coeff, Double_t xlim, Double_t ylim, Int_t steps, TMatrixD& X, TMatrixD& Y, TMatrixD& Z);
void fft_shift(Double_t *fft_magnitude, Int_t window_size);
//TRandom3 randomGen;
void prepare_fft_axis(TH2D *hist, double dt, int window_size);

int main() {
    // Define constants
    Double_t GRATING_PERIOD = 400;
    Double_t WAVELENGTH = 5;
    Double_t TALBOT_LENGTH = TMath::Power(GRATING_PERIOD, 2) / WAVELENGTH;
//position of the centre of the pattern
    Double_t Z0 = TALBOT_LENGTH * NUMBER_OF_PERIODS / 2;
//width of the pattern (Good field region)
    Double_t L = TALBOT_LENGTH * 1000;

    std::cout << TALBOT_LENGTH << " ";

    // Create 2D grid for x and zOthers
    TVectorD x(XYSTEPS), y(XYSTEPS), z(ZSTEPS);
    for (Int_t i = 0; i < XYSTEPS; i++) {
        x[i] = (i-XYSTEPS/2) * (NUMBER_OF_PERIODS * GRATING_PERIOD) / XYSTEPS;
        y[i] = x[i];
    }
    for (Int_t i = 0; i < ZSTEPS; i++) {
        z[i] = i * (NUMBER_OF_PERIODS * TALBOT_LENGTH) / ZSTEPS;
    }

    // Define the third grating plane, starts upright
    TMatrixD third_grating(XYSTEPS, XYSTEPS);
    //third_grating = 1;

    // Tip and rotation angles
    Double_t tip_angle = TMath::Pi() / 12;
    Double_t rotation_angle_x = 0*TMath::Pi() / 12;
    Double_t rotation_angle_y = 0*TMath::Pi() / 12;

    // Define rotation matrices
    TMatrixD R_z(3, 3);
    R_z[0][0] = TMath::Cos(rotation_angle_x);
    R_z[0][1] = -TMath::Sin(rotation_angle_x);
    R_z[0][2] = 0;
    R_z[1][0] = TMath::Sin(rotation_angle_x);
    R_z[1][1] = TMath::Cos(rotation_angle_x);
    R_z[1][2] = 0;
    R_z[2][0] = 0;
    R_z[2][1] = 0;
    R_z[2][2] = 1;

    TMatrixD R_x(3, 3);
    R_x[0][0] = 1;
    R_x[0][1] = 0;
    R_x[0][2] = 0;
    R_x[1][0] = 0;
    R_x[1][1] = TMath::Cos(tip_angle);
    R_x[1][2] = -TMath::Sin(tip_angle);
    R_x[2][0] = 0;
    R_x[2][1] = TMath::Sin(tip_angle);
    R_x[2][2] = TMath::Cos(tip_angle);
    TMatrixD R_y(3, 3);
    R_y[0][0] = TMath::Cos(rotation_angle_y);
    R_y[0][1] = 0;
    R_y[0][2] = TMath::Sin(rotation_angle_y);
    R_y[1][0] = 0;
    R_y[1][1] = 1;
    R_y[1][2] = 0;
    R_y[2][0] = -TMath::Sin(rotation_angle_y);
    R_y[2][1] = 0;
    R_y[2][2] = TMath::Cos(rotation_angle_y);
    //total rotation matrix
    TMatrixD R = R_x * R_y * R_z;
    //mean pattern vector declaration
    TVectorD mean_pattern_values(z.GetNrows());

    //create plotter
    TGraph *graph = new TGraph(XYSTEPS);

    //it is required for each core to have its own random generator. If left at 24 this prevents problem between macOS and linux
#pragma omp parallel private(randomGen)
    //TRandom3 randomGen(omp_get_thread_num());
    TRandom3 randomGen (24);
#pragma omp for
    //generate the pattern as seen by the third grating
    for (Int_t i = 0; i < z.GetNrows(); i++) {
        //std::cout << i << " ";
        Double_t sum_pattern = 0.0;
        Int_t count = 0;
        //add for loop for j and k and then create a true mean
        for (Int_t j = 0; j < x.GetNrows(); j++) {
            for (Int_t k = 0; k < y.GetNrows(); k++) {
                //define the x,y,z point to be rotated and pattern generated
                TMatrixD point(3,1);
                point[0][0] = x[j];
                point[1][0] = y[k];
                point[2][0] = z[i];
                //rotate the point
                TMatrixD rotP = R * point;
                //extract the pattern
                Double_t pattern = Pattern(&randomGen, rotP[0][0], rotP[1][0], rotP[2][0], GRATING_PERIOD, TALBOT_LENGTH, Z0, L);
                sum_pattern += pattern;
                count++;
            }
            mean_pattern_values[i] = sum_pattern / count;

#pragma omp critical
            {
                graph->SetPoint(i, z[i], mean_pattern_values[i]);
            }
        }
        //Double_t pattern = Pattern(x[j], y[k], z[i], GRATING_PERIOD, TALBOT_LENGTH, Z0, L);
        //mean_pattern_values[i] = pattern; // Here I'm assuming that pattern is the mean pattern value
        //std::cout << pattern << " ";

    }

// Create a new TCanvas
    TCanvas *c1 = new TCanvas("c1", "Pattern Plot", 800, 600);

// Draw the graph
    graph->Draw("AL");

// Set titles
    graph->SetTitle("Pattern vs. Z; Z; Pattern");

// Update the canvas
    c1->Update();

    // Save the canvas to a .pdf file
    c1->SaveAs("meanPattern.png");

    std::cout << "mean pattern produced: ";

    /////////////////////////
    //fft
    // Define window size (2*Talbot length)
    Int_t window_size = 2 * RESOLUTION;

// Create a new TVirtualFFT object
    TVirtualFFT *fft = TVirtualFFT::FFT(1, &window_size, "R2C M K");

// Loop over z values with a step of window_size
    std::cout << "Rows: " << z.GetNrows() << " ";
    std::cout << "Window size: " << window_size << " ";

////////////////////////////////

// Create a 2D vector for FFT results
    std::vector<std::vector<Double_t>> fft_results(ZSTEPS - window_size, std::vector<Double_t>(window_size, 0));

#pragma omp parallel
#pragma omp for
    std::cout << "fftresults.size" << fft_results.size() << std::endl;
    for (Int_t i = 0; i < z.GetNrows() - window_size; i++) {
        std::cout << i << " / " << z.GetNrows() - window_size << std::endl;
        // Create an array for the window data
        Double_t window_data[window_size];

        // Fill the window data array
        for (Int_t k = 0; k < window_size; k++) {
            window_data[k] = mean_pattern_values[i + k];
        }

        // Set the data for the FFT
        fft->SetPoints(window_data);

        // Compute the FFT
        fft->Transform();

        // Get the magnitude of the FFT result
        Double_t fft_magnitude[window_size];
        for (Int_t k = 0; k < window_size; k++) {
            Double_t re, im;
            fft->GetPointComplex(k, re, im);
            fft_magnitude[k] = TMath::Sqrt(re * re + im * im);
        }

        //perform fft shift like in numpy
        fft_shift(fft_magnitude, window_size);

        // Store the FFT magnitudes
        for (Int_t k = 0; k < window_size; k++) {
            fft_results[i][k] = fft_magnitude[k];
        }
        std::cout << fft_results[i].size() << " " << std::endl;
    }
    std::cout << "got past the hard bit" << std::endl;


//now we have fft results
// Set the size of the array
    const Int_t N = ZSTEPS - window_size; // replace with your actual array size
    const Int_t M = window_size; // replace with your actual array size

    // Normalize each row to its maximum value
//    for (Int_t i = 0; i < N; ++i) {
//        // Find the maximum value in the row
//        Double_t max_value = fft_results[i][0];
//        for (Int_t j = 0; j < M; ++j) {
//            if (fft_results[i][j] > max_value) {
//                max_value = fft_results[i][j];
//            }
//        }
//
//        // Normalize the row
//        for (Int_t j = 0; j < M; ++j) {
//            fft_results[i][j] /= max_value;
//        }
//    }



    // Initialize the histogram
    TH2D *hist = new TH2D("hist", "FFT Results", ZSTEPS, 0, ZSTEPS, window_size, 0, window_size);

    // Fill the histogram with the values from fft_results
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            //hist->SetBinContent(i+1, j+1, fft_results[i][j]);
            hist->SetBinContent(i+1, j+1, mean_pattern_values[i+j]);
        }
    }

    // Create a canvas to draw on
    TCanvas *canvas = new TCanvas("canvas", "FFT Results", 800, 600);

    // Draw the histogram
    hist->Draw("COLZ");

    // Save the canvas as a PNG file
    canvas->SaveAs("preFft.png");



    // Create a 2D histogram
    TH2D *h2 = new TH2D("h2", "FFT Results", N, 0, N, M, 0, M);

    // Fill the histogram from the array
    for (Int_t i = 0; i < N; ++i) {
        for (Int_t j = 0; j < M; ++j) {
            h2->Fill(i, j, fft_results[i][j]);
        }
    }
    Double_t dt = (NUMBER_OF_PERIODS * TALBOT_LENGTH) / ZSTEPS;

    prepare_fft_axis(h2, dt, window_size);

    //Double_t bin_min = 450;
    //Double_t  bin_max = 550;


    TCanvas *c3 = new TCanvas("c3", "FFT Results", 800, 600);

// Get total number of bins in the histogram
    Int_t total_bins = h2->GetNbinsY();

// Calculate the midpoint bin
    Int_t midpoint_bin = total_bins / 2;

// Define a window around the midpoint bin
    Int_t window = 25; // Change this value as needed
    Int_t bin_min = midpoint_bin - window;
    Int_t bin_max = midpoint_bin + window;

// Ensure bin_min and bin_max are valid
    if (bin_min < 1) bin_min = 1;
    if (bin_max > total_bins) bin_max = total_bins;

    h2->GetYaxis()->SetRange(bin_min, bin_max);
    h2->Draw("COLZ");
    c3->SaveAs("fft_results.png");


    std::cout << "simulation is finished " << std::endl;
    return 0;


}

Double_t Pattern(TRandom3* randomGen1, Double_t x, Double_t y, Double_t z, Double_t GRATING_PERIOD, Double_t TALBOT_LENGTH, Double_t Z0, Double_t L) {
    Double_t noise = (randomGen1->Rndm() - 0.5) * 0;  // Random noise in range [-0.05, 0.05]
    return (TMath::Power(TMath::Cos(2 * TMath::Pi() / GRATING_PERIOD * x), 2) *
            TMath::Power(TMath::Cos(2 * TMath::Pi() / TALBOT_LENGTH * (z - Z0)), 2) *
            TMath::Exp(-(z - Z0) * (z - Z0) / (L * L))) + noise;
}

void plane_to_grid(TVectorD coeff, Double_t xlim, Double_t ylim, Int_t steps, TMatrixD& X, TMatrixD& Y, TMatrixD& Z) {
    Double_t A = coeff[0], B = coeff[1], C = coeff[2];
    TVectorD x = TVectorD(steps);
    TVectorD y = TVectorD(steps);
    for (Int_t i = 0; i < steps; i++) {
        x[i] = -xlim + 2*xlim/steps*i;
        y[i] = -ylim + 2*ylim/steps*i;
    }
    X.ResizeTo(steps, steps);
    Y.ResizeTo(steps, steps);
    Z.ResizeTo(steps, steps);
    for (Int_t i = 0; i < steps; i++) {
        for (Int_t j = 0; j < steps; j++) {
            X[i][j] = x[i];
            Y[i][j] = y[j];
            Z[i][j] = -(A*x[i] + B*y[j]) / C;
        }
    }
}

//this function does the same as a numpy fftshift it can be applied in the main or not depending on preference
//this func shifts a std::vector not a 2d vector so perform the shift before filling the 2d vec or TH2D
void fft_shift(Double_t *fft_magnitude, Int_t window_size) {
    Double_t shifted[window_size];
    Int_t half_window = window_size / 2;
    for (Int_t k = 0; k < window_size; ++k) {
        Int_t index_shifted = (k + half_window) % window_size;
        shifted[index_shifted] = fft_magnitude[k];
    }
    for (Int_t k = 0; k < window_size; ++k) {
        fft_magnitude[k] = shifted[k];
    }
}
//function to generate the correct axis for the frequency after the fft
void prepare_fft_axis(TH2D *hist, double dt, int window_size) {
    double df = 1.0 / (window_size * dt);
    double f_max = 1.0 / (2.0 * dt);
    std::cout << "dt: " << dt << std::endl;

    //adding the correct axis
    int label_step = window_size/10;
    for (int k = 1; k <= window_size; k+=label_step) {
        double freq = -f_max + (k - 1) * df;
        hist -> GetYaxis()->SetBinLabel(k, Form("%.4f", freq));
    }
}
