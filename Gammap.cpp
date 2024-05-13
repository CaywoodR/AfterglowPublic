#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

double Ej = 5.0E+50; // Initial energy in ergs
double G1 = 100.0;      // Initial Gamma
double t1 = 0.1;
double t2 = 0.1;

double rdeg(double x)
{
return M_PI/180*x;
}

double E(double x) // Energy curve takes in radians
{
    if (rdeg(x) < 0.1) {
        return Ej;
    } else {
        return Ej * pow(rdeg(x) / t2, -4.5);
    }
}

double G(double x) // Initial Gamma: takes in an angle in radians
{
    if (rdeg(x) < 0.1) {
        return G1;
    } else {
        return 1.0 + (G1 - 1.0) * pow(rdeg(x) / t1, -4.5);
    }
}

int main() {
    std::ofstream outputFile("gamma.txt"); // Open a file for writing

    if (!outputFile.is_open()) { // Check if file is opened successfully
        std::cerr << "Error opening the file." << std::endl;
        return 1;
    }

    for (double x = 0.0; x < 100; x++) {
        double energy = E(x);
        double gamma = G(x);

        outputFile << std::fixed << std::setprecision(6) << x << " " << energy << " " << gamma << std::endl;
    }

    outputFile.close(); // Close the file
    return 0;
}

