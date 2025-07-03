/*
HDF5Functions.cpp
functions for outputting to HDF5
Updates by Drew Murray
2/15/2017
*/

//libraries
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
#include <string>
using std::string;
#include <H5Cpp.h> 
using H5::DataType;
using H5::Exception;
using H5::H5File;
using H5::DataSpace;
using H5::DataSet;
using H5::FileIException;
using H5::DataSetIException;
using H5::DataSpaceIException;
using H5::DataTypeIException;
using H5::PredType;
#include <cmath>

//class includes
#include "HDF5Functions.h"
#include "Grid.h"
#include "Block.h"

//Write the physical coordinates and constants to file (runs only at the beginning of the simulation)
int dumpGridToHDF5File(string outputFolder, Grid* myGrid, int timestep, float time, int rank)
{
	int NUM_BLOCKS[3];
	int NUM_ELEMENTS[3];
	int ORDER=myGrid->getOrder();
	Block* readBlock;

	for (int itr=0;itr<3;itr++)
	{
		NUM_BLOCKS[itr]=myGrid->getNumBlocks()[itr];
		NUM_ELEMENTS[itr]=myGrid->getNumElements()[itr];
	}
	int size=NUM_BLOCKS[0]*NUM_BLOCKS[1]*NUM_BLOCKS[2];
	
	int dimx = NUM_ELEMENTS[0] * (ORDER + 1);
	int dimy = NUM_ELEMENTS[1] * (ORDER + 1);
	int dimz = NUM_ELEMENTS[2] * (ORDER + 1);
	if (size==1)
	{
		dimx = NUM_BLOCKS[0] * NUM_ELEMENTS[0] * (ORDER + 1);
		dimy = NUM_BLOCKS[1] * NUM_ELEMENTS[1] * (ORDER + 1);
		dimz = NUM_BLOCKS[2] * NUM_ELEMENTS[2] * (ORDER + 1);
	}

	vector<vector<vector<vector<double>>>> outGrid(dimx, vector<vector<vector<double>>>(dimy, vector<vector<double>>(dimz, vector<double>(3))));
	vector<vector<vector<double>>> outGridC0(dimx, vector<vector<double>>(dimy, vector<double>(dimz)));
	vector<vector<vector<double>>> outGridRho0(dimx, vector<vector<double>>(dimy, vector<double>(dimz)));
	vector<vector<vector<double>>> outGridBOverA(dimx, vector<vector<double>>(dimy, vector<double>(dimz)));
	vector<vector<vector<double>>> outGridDiff(dimx, vector<vector<double>>(dimy, vector<double>(dimz)));
	//
	for (int ii = 0; ii < dimx; ii++){
		for (int jj = 0; jj < dimy; jj++){
			for (int kk = 0; kk < dimz; kk++){
				int blk[3];
				int ele[3];
				int nod[3];
				if (size==1)
				{
					blk[0] = ii / NUM_ELEMENTS[0] / (ORDER + 1);
					blk[1] = jj / NUM_ELEMENTS[1] / (ORDER + 1);
					blk[2] = kk / NUM_ELEMENTS[2] / (ORDER + 1);
					ele[0] = (ii / (ORDER + 1)) % NUM_ELEMENTS[0]+1;//+1 to adjust for fact that 0th index is halo
					ele[1] = (jj / (ORDER + 1)) % NUM_ELEMENTS[1]+1;
					ele[2] = (kk / (ORDER + 1)) % NUM_ELEMENTS[2]+1;
					readBlock = &(myGrid->mBlocks[blk[0]][blk[1]][blk[2]]);
				}
				else
				{
					readBlock = &(myGrid->mBlock);
					ele[0] = (ii / (ORDER + 1))+1;//+1 to adjust for fact that 0th index is halo
					ele[1] = (jj / (ORDER + 1))+1;
					ele[2] = (kk / (ORDER + 1))+1;
				}
				nod[0] = ii % (ORDER + 1);
				nod[1] = jj % (ORDER + 1);
				nod[2] = kk % (ORDER + 1);
				outGrid[ii][jj][kk][0] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputCoords(nod)[0];
				outGrid[ii][jj][kk][1] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputCoords(nod)[1];
				outGrid[ii][jj][kk][2] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputCoords(nod)[2];
				outGridC0[ii][jj][kk] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getC0();
				outGridRho0[ii][jj][kk] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getRho0();
				outGridBOverA[ii][jj][kk] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getBoverA();
				outGridDiff[ii][jj][kk] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getDiffusionConstant();
				
				
			}
		}
	}
	string outputFilename0 = outputFolder+"/raw/coords_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";
	string outputFilename1 = outputFolder+"/raw/ba_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";
	string outputFilename2 = outputFolder+"/raw/diff_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";
	string outputFilename3 = outputFolder+"/raw/rho0_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";
	string outputFilename4 = outputFolder+"/raw/c0_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";
	cout << " - Printing to "<<outputFilename0 << "\n";
	try
	{
		hdf5outputVector(outGrid, outputFilename0, "coords", dimx, dimy, dimz,time);
		hdf5outputScalar(outGridBOverA, outputFilename1, "ba", dimx, dimy, dimz,time);
		hdf5outputScalar(outGridDiff, outputFilename2, "diff", dimx, dimy, dimz,time);
		hdf5outputScalar(outGridRho0, outputFilename3, "rho0", dimx, dimy, dimz,time);
		hdf5outputScalar(outGridC0, outputFilename4, "c0", dimx, dimy, dimz,time);
	}

	// catch failure caused by the H5File operations
	catch (FileIException error)
	{
		cout<<"HDF5 FileIException error\n";
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		cout<<"HDF5 DataSetIException error\n";
		return -2;
	}
	// catch failure caused by the DataSpace operations
	catch (DataSpaceIException error)
	{
		cout<<"HDF5 DataSpaceIException error\n";
		return -3;
	}
	// catch failure caused by the DataSpace operations
	catch (DataTypeIException error)
	{
		cout<<"HDF5 DataTypeIException error\n";
		return -4;
	}
	return 0;
}

//write the state variables to file
int dumpToHDF5File(string outputFolder, Grid* myGrid, int timestep, float time, int rank)
{
	int NUM_BLOCKS[3];
	int NUM_ELEMENTS[3];
	int ORDER=myGrid->getOrder();
	Block* readBlock;

	for (int itr=0;itr<3;itr++)
	{
		NUM_BLOCKS[itr]=myGrid->getNumBlocks()[itr];
		NUM_ELEMENTS[itr]=myGrid->getNumElements()[itr];
	}
	int size=NUM_BLOCKS[0]*NUM_BLOCKS[1]*NUM_BLOCKS[2];
	
	//size>1
	int dimx = NUM_ELEMENTS[0] * (ORDER + 1);
	int dimy = NUM_ELEMENTS[1] * (ORDER + 1);
	int dimz = NUM_ELEMENTS[2] * (ORDER + 1);
	if (size==1)
	{
		dimx = NUM_BLOCKS[0] * NUM_ELEMENTS[0] * (ORDER + 1);
		dimy = NUM_BLOCKS[1] * NUM_ELEMENTS[1] * (ORDER + 1);
		dimz = NUM_BLOCKS[2] * NUM_ELEMENTS[2] * (ORDER + 1);
	}

	vector<vector<vector<double>>> outRho(dimx, vector<vector<double>>(dimy, vector<double>(dimz)));
	vector<vector<vector<vector<double>>>> outVelocity(dimx, vector<vector<vector<double>>>(dimy, vector<vector<double>>(dimz, vector<double>(3))));
	double maxDRho = 0.;
	for (int ii = 0; ii < dimx; ii++){
		for (int jj = 0; jj < dimy; jj++){
			for (int kk = 0; kk < dimz; kk++){
				int blk[3];
				int ele[3];
				int nod[3];
				if (size==1)
				{
					blk[0] = ii / NUM_ELEMENTS[0] / (ORDER + 1);
					blk[1] = jj / NUM_ELEMENTS[1] / (ORDER + 1);
					blk[2] = kk / NUM_ELEMENTS[2] / (ORDER + 1);
					ele[0] = (ii / (ORDER + 1)) % NUM_ELEMENTS[0]+1;//+1 to adjust for fact that 0th index is halo
					ele[1] = (jj / (ORDER + 1)) % NUM_ELEMENTS[1]+1;
					ele[2] = (kk / (ORDER + 1)) % NUM_ELEMENTS[2]+1;
					readBlock = &(myGrid->mBlocks[blk[0]][blk[1]][blk[2]]);
				}
				else
				{
					readBlock = &(myGrid->mBlock);
					ele[0] = (ii / (ORDER + 1))+1;//+1 to adjust for fact that 0th index is halo
					ele[1] = (jj / (ORDER + 1))+1;
					ele[2] = (kk / (ORDER + 1))+1;
				}
				nod[0] = ii % (ORDER + 1);
				nod[1] = jj % (ORDER + 1);
				nod[2] = kk % (ORDER + 1);

				outRho[ii][jj][kk] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputRho(nod);
				if(std::abs(outRho[ii][jj][kk])>maxDRho || std::isnan(outRho[ii][jj][kk])) maxDRho=std::abs(outRho[ii][jj][kk]);
				outVelocity[ii][jj][kk][0] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputVelocity(nod)[0];
				outVelocity[ii][jj][kk][1] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputVelocity(nod)[1];
				outVelocity[ii][jj][kk][2] = readBlock->mElements[ele[0]][ele[1]][ele[2]].getOutputVelocity(nod)[2];
				
			}
		}
	}

	string outputFilename1 = outputFolder+"/raw/dRho_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";
	string outputFilename2 = outputFolder+"/raw/vel_"+to_string((long long int)rank)+"_"+to_string((long long int)timestep) + ".h5";

	cout << " - Printing to "<<outputFilename1 << "; Max dRho="<<maxDRho<<"\n";

	try
	{
		hdf5outputScalar(outRho, outputFilename1, "rho", dimx, dimy, dimz,time);
		hdf5outputVector(outVelocity, outputFilename2, "velocity", dimx, dimy, dimz,time);
	}

	// catch failure caused by the H5File operations
	catch (FileIException error)
	{
		cout<<"HDF5 FileIException error\n";
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		cout<<"HDF5 DataSetIException error\n";
		return -2;
	}
	// catch failure caused by the DataSpace operations
	catch (DataSpaceIException error)
	{
		cout<<"HDF5 DataSpaceIException error\n";
		return -3;
	}
	// catch failure caused by the DataSpace operations
	catch (DataTypeIException error)
	{
		cout<<"HDF5 DataTypeIException error\n";
		return -4;
	}
	return 0;
	
}

int hdf5outputVector(vector<vector<vector<vector<double>>>> data4D, string outputFilename, string Datasetname, int Dimx, int Dimy, int Dimz, float time){
	try {
		Exception::dontPrint();
		// HDF5 output starts here
		H5File file(outputFilename, H5F_ACC_TRUNC);
		H5std_string DATASET_NAME(Datasetname);          //convert a string to an H5 string JFK 3/4/17

		// dataset dimensions
		hsize_t dimsf[4];
		dimsf[0] = 3;
		dimsf[1] = Dimx;
		dimsf[2] = Dimy;
		dimsf[3] = Dimz;
		
		DataSpace dataspace(4, dimsf);

		DataType datatype(H5::PredType::NATIVE_DOUBLE);
		DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);
		// write the data, don't switch order x,y,z => z,y,x
		double data[3][Dimx][Dimy][Dimz];
		for (int h = 0; h < 3; h++)
			for (int i = 0; i < Dimx; i++)
				for (int j = 0; j < Dimy; j++)
					for (int k = 0; k < Dimz; k++)
						data[h][i][j][k] = data4D[i][j][k][h];

		dataset.write(data, H5::PredType::NATIVE_DOUBLE);


		// close the dataset
		dataset.close();
		dataspace.close();
		
	    // Write out time. 
		hsize_t dimsf2[1];
		dimsf2[0] = 1;
		double data2[1];
		data2[0] = time;
		DataSpace dataspace2(1, dimsf2);
		H5std_string DATASET_NAME2("time");
		DataSet dataset2 = file.createDataSet(DATASET_NAME2, datatype, dataspace2);
		dataset2.write(data2, H5::PredType::NATIVE_DOUBLE);
		dataset2.close();
		dataspace2.close();	
		
		file.close();
	}

	// catch failure caused by the H5File operations
	catch (FileIException error)
	{
		cout<<"HDF5 FileIException error\n";
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		cout<<"HDF5 DataSetIException error\n";
		return -2;
	}
	// catch failure caused by the DataSpace operations
	catch (DataSpaceIException error)
	{
		cout<<"HDF5 DataSpaceIException error\n";
		return -3;
	}
	// catch failure caused by the DataSpace operations
	catch (DataTypeIException error)
	{
		cout<<"HDF5 DataTypeIException error\n";
		return -4;
	}
	return 0;
}
int hdf5outputScalar(vector<vector<vector<double>>> data3D, string outputFilename, string Datasetname, int Dimx, int Dimy, int Dimz, float time){
	try {
		Exception::dontPrint();
		// HDF5 output starts here
		H5File file(outputFilename, H5F_ACC_TRUNC);
		H5std_string DATASET_NAME(Datasetname);         //convert a string to an H5 string JFK 3/4/17
		// dataset dimensions
		hsize_t dimsf[3];
		dimsf[0] = Dimx;
		dimsf[1] = Dimy;
		dimsf[2] = Dimz;
		DataSpace dataspace(3, dimsf);

		DataType datatype(H5::PredType::NATIVE_DOUBLE);
		//DataSet dataset = file.createDataSet(Datasetname, datatype, dataspace);
		DataSet dataset = file.createDataSet(DATASET_NAME, datatype, dataspace);

		// write the data, don't switch order x,y,z => z,y,x
		double data[Dimx][Dimy][Dimz];
		for (int i = 0; i < Dimx; i++)
			for (int j = 0; j < Dimy; j++)
				for (int k = 0; k < Dimz; k++)
					data[i][j][k] = data3D[i][j][k];

		dataset.write(data, H5::PredType::NATIVE_DOUBLE);

		// close the dataset
		dataset.close();
		dataspace.close();
		
		
		
		// Write out time. 
		hsize_t dimsf2[1];
		dimsf2[0] = 1;
		double data2[1];
		data2[0] = time;
		DataSpace dataspace2(1, dimsf2);
		H5std_string DATASET_NAME2("time");
		DataSet dataset2 = file.createDataSet(DATASET_NAME2, datatype, dataspace2);
		dataset2.write(data2, H5::PredType::NATIVE_DOUBLE);
		dataset2.close();
		dataspace2.close();
		
				
		file.close();
	}

	// catch failure caused by the H5File operations
	catch (FileIException error)
	{
		cout<<"HDF5 FileIException error\n";
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		cout<<"HDF5 DataSetIException error\n";
		return -2;
	}
	// catch failure caused by the DataSpace operations
	catch (DataSpaceIException error)
	{
		cout<<"HDF5 DataSpaceIException error\n";
		return -3;
	}
	// catch failure caused by the DataSpace operations
	catch (DataTypeIException error)
	{
		cout<<"HDF5 DataTypeIException error\n";
		return -4;
	}

	return 0;
}
int hdf5input(vector<vector<vector<double>>>& data3D, string inputFilename, string Datasetname, int Dimx, int Dimy, int Dimz){
	int i, j, k;
	double data_out[Dimz][Dimy][Dimx]; /* output buffer */
	for (i = 0; i < Dimz; i++)
	{
		for (j = 0; j < Dimy; j++)
		{
			for (k = 0; k < Dimx; k++)
				data_out[i][j][k] = 0;
		}
	}
	try
	{
		Exception::dontPrint();
		H5File file(inputFilename.c_str(), H5F_ACC_RDONLY);
		DataSet dataset = file.openDataSet(Datasetname.c_str());

		H5T_class_t type_class = dataset.getTypeClass();
		cout << type_class << "\n";

		DataSpace dataspace = dataset.getSpace();
		hsize_t dimsm[3];              /* memory space dimensions */
		dimsm[0] = Dimz;
		dimsm[1] = Dimy;
		dimsm[2] = Dimx;
		hsize_t rank = dataspace.getSimpleExtentDims(dimsm, NULL);


		DataSpace memspace(rank, dimsm);

		dataset.read(data_out, PredType::NATIVE_DOUBLE, dataspace);
		//dataset.read(data_out, PredType::NATIVE_DOUBLE, memspace, dataspace);
		for (int i = 0; i < Dimz; i++)
			for (int j = 0; j < Dimy; j++)
				for (int k = 0; k < Dimx; k++)
					data3D[k][j][i] = data_out[i][j][k];
		file.close();

	}
	// catch failure caused by the H5File operations
	catch (FileIException error)
	{
		cout<<"HDF5 FileIException error\n";
		return -1;
	}
	// catch failure caused by the DataSet operations
	catch (DataSetIException error)
	{
		cout<<"HDF5 DataSetIException error\n";
		return -2;
	}
	// catch failure caused by the DataSpace operations
	catch (DataSpaceIException error)
	{
		cout<<"HDF5 DataSpaceIException error\n";
		return -3;
	}
	// catch failure caused by the DataSpace operations
	catch (DataTypeIException error)
	{
		cout<<"HDF5 DataTypeIException error\n";
		return -4;
	}
	return 0;  // successfully terminated

}

