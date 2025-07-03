/*
HDF5Functions.h
functions for outputting to HDF5
*/

//circular build protection
#ifndef HDF5Functions_H
#define HDF5Functions_H

//libraries
#include <vector>
using std::vector;
#include <string>
using std::string;

//class includes
#include "Grid.h"

int dumpGridToHDF5File(string outputFolder, Grid* myGrid, int timestep, float time, int rank);
int dumpToHDF5File(string outputFolder, Grid* myGrid, int timestep, float time, int rank);

int hdf5outputVector(vector<vector<vector<vector<double>>>> data4D, string outputFilename, string Datasetname, int Dimx, int Dimy, int Dimz, float time);

int hdf5outputScalar(vector<vector<vector<double>>> data3D, string outputFilename, string Datasetname, int Dimx, int Dimy, int Dimz, float time);

int hdf5input(vector<vector<vector<double>>>& data3D, string inputFilename, string Datasetname, int Dimx, int Dimy, int Dimz);

#endif
