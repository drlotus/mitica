#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include <string>
#include "TSpline.h"

struct DataPoint {
    double x;
    double y;
};

std::vector<DataPoint> parseCSV(const std::string& filename);
TSpline3* performInterpolation(const std::string& filename);
TSpline3* getInterpolationSpline(const std::string& filename);
void plotTGraph(TGraph* graph);
void saveTableToFile(const TSpline3* spline, const std::string& filename);

#endif // INTERPOLATION_H