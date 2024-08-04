/******************************************************************************
*                                                                             *
*            vHLLE : a 3D viscous hydrodynamic code                           *
*            by Iurii Karpenko                                                *
*  contact:  yu.karpenko@gmail.com                                            *
*  For the detailed description please refer to:                              *
*  Comput. Phys. Commun. 185 (2014), 3016   arXiv:1312.4160                   *
*                                                                             *
*  This code can be freely used and redistributed, provided that this         *
*  copyright appear in all the copies. If you decide to make modifications    *
*  to the code, please contact the authors, especially if you plan to publish *
* the results obtained with such modified code. Any publication of results    *
* obtained using this code must include the reference to                      *
* arXiv:1312.4160 [nucl-th] or the published version of it.                   *
*                                                                             *
*******************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm> // for std::sort
#include <TGraph.h>
#include <TSpline.h>
#include <TCanvas.h>

// Structure to hold x and y values
struct DataPoint {
    double x;
    double y;
};

// Function to parse CSV file and extract data points
std::vector<DataPoint> parseCSV(const std::string& filename) {
    std::vector<DataPoint> dataPoints;
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return dataPoints;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string xStr, yStr;

        if (std::getline(iss, xStr, ',') && std::getline(iss, yStr, ',')) {
            double x = std::stod(xStr);
            double y = std::stod(yStr);
            dataPoints.push_back({x, y});
        }
    }

    file.close();
    return dataPoints;
}

void plotTGraph(TGraph* graph) {
    // Create a canvas for the plot
    TCanvas* canvas = new TCanvas("canvas", "Graph Plot", 800, 600);

    // Set the graph attributes
    graph->SetLineColor(kBlue);
    graph->SetLineWidth(2);

    // Draw the graph on the canvas
    graph->Draw("AL");

    // Update the canvas
    canvas->Update();

    // Wait for user input to exit the program
    canvas->WaitPrimitive();
}

TSpline3* performInterpolation(const std::string& filename) {
    std::vector<DataPoint> dataPoints = parseCSV(filename);
    if (dataPoints.empty()) {
        std::cerr << "No data points found." << std::endl;
        return nullptr;
    }

    // Sort data points based on x-values in case the data file is not sorted
    std::sort(dataPoints.begin(), dataPoints.end(), [](const DataPoint& a, const DataPoint& b) {
        return a.x < b.x;
    });

    size_t numDataPoints = dataPoints.size();
    TGraph graph(numDataPoints);
    for (size_t i = 0; i < numDataPoints; ++i) {
        graph.SetPoint(i, dataPoints[i].x, dataPoints[i].y);
    }
    TSpline3* spline = new TSpline3("spline", &graph);

    return spline;
}

// Function definition for obtaining the TSpline3 object
TSpline3* getInterpolationSpline(const std::string& filename) {
    // Call the function in the interpolation file to perform the interpolation
    // and return the TSpline3 object
    return performInterpolation(filename);
}

void saveTableToFile(const TSpline3* spline, const std::string& filename) {
    std::ofstream file(filename);

    // Generate a table of values between x=0.1 and x=20
    double x = 0.1;
    while (x <= 20.0) {
        double y = spline->Eval(x);
        file << x << "\t" << y << std::endl;
        x += 0.05;
    }

    file.close();
}