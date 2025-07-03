/*
driver.cpp
main driver for 3D-DG
Drew Murray
11/29/16
*/


//circular build protection
#ifndef DRIVER
#define DRIVER

//libraries
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::cerr;
#include <fstream>
using std::ofstream;
using std::ifstream;
#include <sstream>
using std::istringstream;
#include<ctime> //TODO: CONVERSION: must be replaced with mpitime
#include <sys/stat.h>
#include <sys/types.h>
#include <mpi.h>
#include <unistd.h>

//class includes
#include "Element.h"
#include "Node.h"
#include "Basis.h"
#include "BCParam.h"
#include "Grid.h"
#include "HDF5Functions.h"

//Error codes:
const int ERRCODE_SUCESSFUL_COMPLETION=0;
const int ERRCODE_INVALID_ARGUMENT=1;
const int ERRCODE_INVALID_PARAMETER=2;
const int ERRCODE_CANNOT_READ=3;
const int ERRCODE_CANNOT_WRITE=4;

//simulation constants (THESE ARE ONLY DEFAULTS, they will be replaced by the read-in from the file)
unsigned int NUM_TIME_SUBSTEPS=1;//number of timesteps code runs between output
unsigned int NUM_OUTSTEPS=1;//number of times outputed
unsigned int FREQ_TERMINAL=1;//how often to print status to terminal
bool BOOL_FILE=true;//whether to dump data to file
unsigned int NUM_ELEMENTS[3] = {4,3,2};
unsigned int NUM_BLOCKS[3] = {5,4,6};
unsigned short int ORDER = 4;
float ELEMENT_WIDTH[3]={1.F,2.F,3.F};//TODO coord float or double?
float STARTING_CORNER[3]={0.25F,0.5F,0.75F};//TODO coord float or double?
float DT = .001F;
float C0 = 1500.0F;
float RHO0 = 1000.0F;
float B_OVER_A = 0.0F;
float DIFF_CONST = 0.0F;
float ATTENUATION = 1.F;
bool eval_advect = false;
bool eval_diff = false;
bool eval_sponge = false;

//printed when command line arguments are incorrect
const string NEEDED_ARGUMENTS = "	PARAMETER_FILENAME : path to the file defining simulation parameters\n	OUTPUT_FILENAME : path to the file that the data will be written to";

int main(int argc, char* argv[])
{
///////////////////////////////////////////////////////////////////////////////
//initialize
	MPI_Init (&argc, &argv);      // initializes MPI
	//Turn on error return (which makes it continue instead of crash with invalid messages, so that the error codes can be checked.
	//https://www-internal.lsc.phy.cam.ac.uk/nmm1/MPI/Notes/notes_06.pdf
	int error = MPI_Comm_set_errhandler (MPI_COMM_WORLD,MPI_ERRORS_RETURN);
	//void *value; int isSet; MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,&value,&isSet); cout<<"Largest tag is "<<*(int *)value<<"\n";

	int rank, size;
	MPI_Comm_rank (MPI_COMM_WORLD, &rank); // get current MPI-process ID. O, 1, ...
	MPI_Comm_size (MPI_COMM_WORLD, &size); // get the total number of processes
	MPI_Barrier(MPI_COMM_WORLD);
	
	cout << " - Thread "<<rank<<" created\n";
	usleep(200);
	MPI_Barrier(MPI_COMM_WORLD);

///////////////////////////////////////////////////////////////////////////////
//Parse command line arguments
	string parameterFilename, outputFolder, logFilename;
	if (argc != 3)
	{
		cerr<<"Not enough command line arguments.\n"<<NEEDED_ARGUMENTS;
		return ERRCODE_INVALID_ARGUMENT;//errorcode
	}
	else if(argc==3)
	{
		parameterFilename=argv[1];
		outputFolder=argv[2];
		logFilename = outputFolder + "/_log.txt";
	}
	if(rank==0)cout << " - Parameters from   : "<<parameterFilename<<"\n";
	if(rank==0)cout << " - Log will be written to : "<<logFilename<<"\n";
	if(rank==0)cout << " - Results in same folder as log\n";


///////////////////////////////////////////////////////////////////////////////
//I/O management
	if(rank==0)cout << "Setting up I/O\n";
	ifstream parameterFile;
	string line;


	if (rank == 0)cout << " - Checking Folder\n";
	mkdir(outputFolder.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	mkdir((outputFolder+"/raw").c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	if (rank == 0)cout << " - Preparing to write\n";
	ofstream logFile;
	logFile.open(logFilename);
	if (!logFile.is_open())
	{
		cerr << "Unable to open log file " << logFilename << "\n";
		MPI_Finalize();
		return ERRCODE_CANNOT_WRITE;
	}

	if(rank==0)cout << " - Reading parameters\n";
	parameterFile.open(parameterFilename);
	BCParam myBCParam = BCParamDefault;//VERY IMPORTANT to ensure no memory garbage gets in since not all vals are overwritten

	if (parameterFile.is_open())
  	{
		string parameterName;
		string parameterValues;
		unsigned int position;
		logFile << "Parameters used:\n{\n";
		while ( getline (parameterFile,line) )
		{
			logFile << "  "<< line << "\n";
			if(line[0]!='#'&&line.length()>=2)//allow commenting
			{
				position = (unsigned int)line.find(" = ");
				parameterName=line.substr(0,position);
				parameterValues=line.substr(position+3);
				istringstream iss(parameterValues);
				if(parameterName=="NUM_TIME_SUBSTEPS") iss>>NUM_TIME_SUBSTEPS;
				else if(parameterName=="NUM_OUTSTEPS") iss>>NUM_OUTSTEPS;
				else if(parameterName=="FREQ_TERMINAL") iss>>FREQ_TERMINAL;
				else if(parameterName=="BOOL_FILE") iss>>BOOL_FILE;
				else if(parameterName=="C0") iss>>C0;
				else if(parameterName=="RHO0") iss>>RHO0;
				else if(parameterName=="B_OVER_A") iss>>B_OVER_A;
				else if(parameterName=="DIFF_CONST") iss>>DIFF_CONST;
				else if(parameterName=="NUM_ELEMENTS") iss>>NUM_ELEMENTS[0]\
													>>NUM_ELEMENTS[1]>>NUM_ELEMENTS[2];
				else if(parameterName=="NUM_BLOCKS") iss>>NUM_BLOCKS[0]\
													>>NUM_BLOCKS[1]>>NUM_BLOCKS[2];
				else if(parameterName=="ORDER") iss>>ORDER;
				else if(parameterName=="ELEMENT_WIDTH") iss>>ELEMENT_WIDTH[0]\
													>>ELEMENT_WIDTH[1]>>ELEMENT_WIDTH[2];
				else if(parameterName=="STARTING_CORNER") iss>>STARTING_CORNER[0]\
													>>STARTING_CORNER[1]>>STARTING_CORNER[2];
				else if(parameterName=="DT") iss>>DT;
				else if(parameterName=="EVALUATE_DIFFUSION") iss>>eval_diff;
				else if(parameterName=="EVALUATE_ADVECTION") iss>>eval_advect;
				else if(parameterName=="EVALUATE_SPONGE") iss>>eval_sponge;
			//BC parameters
				else if(parameterName=="BC_IS_PERIODIC") iss>>myBCParam.IS_PERIODIC[0]>>myBCParam.IS_PERIODIC[1]>>myBCParam.IS_PERIODIC[2];
				else if(parameterName=="BC_SIGMA_MAX") iss >> myBCParam.sigmaMax;
				else if(parameterName=="BC_SPONGE_WIDTH") iss >> myBCParam.spongeWidth[0] >> myBCParam.spongeWidth[1] >> myBCParam.spongeWidth[2] >> myBCParam.spongeWidth[3] >> myBCParam.spongeWidth[4] >> myBCParam.spongeWidth[5];
				else if(parameterName=="IC_X_GAUSS") iss>>myBCParam.IC_X_GAUSS;
				else if(parameterName=="IC_Y_GAUSS") iss>>myBCParam.IC_Y_GAUSS;
				else if(parameterName=="IC_Z_GAUSS") iss>>myBCParam.IC_Z_GAUSS;
				else if(parameterName=="IC_DRHO_AMPLITUDE") iss>>myBCParam.IC_drho_amplitude;
				else if(parameterName=="IC_VX_AMPLITUDE") iss>>myBCParam.IC_vx_amplitude;
				else if(parameterName=="IC_VY_AMPLITUDE") iss>>myBCParam.IC_vy_amplitude;
				else if(parameterName=="IC_VZ_AMPLITUDE") iss>>myBCParam.IC_vz_amplitude;
				else if(parameterName=="IC_SIGMA") iss>>myBCParam.IC_sigma;
				else if(parameterName=="IC_CENTER_X") iss>>myBCParam.IC_center_x;
				else if(parameterName=="IC_CENTER_Y") iss>>myBCParam.IC_center_y;
				else if(parameterName=="IC_CENTER_Z") iss>>myBCParam.IC_center_z;
				else cerr<<"WARNING: IGNORING INVALID PARAMETER LINE: \""<<line<<"\"\n";
			}
		}
		parameterFile.close();
		logFile << "}\n";
		logFile.close();
  	}
  	else
	{
		cerr << "Unable to open parameter file " << parameterFilename <<"\n";
		MPI_Finalize();
		return ERRCODE_CANNOT_READ;
	}

	if(rank==0)cout<< " - Reading basis\n";
	Basis myBasis("hardcoded_coeff.txt","hardcoded_coeff2.txt");//precalculated constants from coefficient basis functions, for each order
	Basis* basisPtr = &myBasis;

	const unsigned int NUM_TIMESTEPS=NUM_OUTSTEPS*NUM_TIME_SUBSTEPS;


///////////////////////////////////////////////////////////////////////////////
//initialize elements

	if(rank==0)cout << "Creating mesh\n";
	if(rank==0)cout << " - Initializing grid using matrices from 'matrices.txt'\n";
	Grid myGrid("matrices.txt", &myBCParam, ORDER, NUM_BLOCKS, NUM_ELEMENTS, ELEMENT_WIDTH, STARTING_CORNER, basisPtr, RHO0, C0, B_OVER_A, DIFF_CONST,eval_advect,eval_diff,eval_sponge);
	Grid myGrid2(myGrid);//create copy
	Grid k1(myGrid);//create copy
	Grid k2(myGrid);//create copy
	//prevent time-consuming copying by switching pointers instead of memory
	Grid* oldGrid = &myGrid;
	Grid* newGrid = &myGrid2;
	Grid* swapGrid = &myGrid;
	if(rank==0)
	{
		cout << " - Logging grid\n";
	}
	if (BOOL_FILE) dumpGridToHDF5File(outputFolder, newGrid, 0, 0.0,rank);
	if(rank==0)
	{
		cout << " - Logging initial conditions\n";
	}
	if (BOOL_FILE) dumpToHDF5File(outputFolder, newGrid, 0,0.,rank);
	

///////////////////////////////////////////////////////////////////////////////
//run simulation
	float t=0;
	double beginTime = MPI_Wtime();
	double prevTime = MPI_Wtime();
	double endTime;
	

	//begin timestepping
	MPI_Barrier(MPI_COMM_WORLD);
	beginTime = MPI_Wtime();
	prevTime = MPI_Wtime();
	if(rank==0)cout << "Timestepping\n";
	
	// Note: Complexity of 3D DG using tensor product hexes is O (N_e*order^4).  Assuming Nx=Ny=Nz=N_e^(1/3)*order, the constant
	// C scales with order.  Also, for an explicit method, Nt is proportional to order^2.  
	for(unsigned int i = 0;i<NUM_TIMESTEPS;i+=NUM_TIME_SUBSTEPS) //0222
	{
		if(rank==0)
		{
			
			endTime=MPI_Wtime();
			if (i>0)cout << "    - Elapsed runtime:" << int(endTime - beginTime) << "s\n";
			if (i>0)cout << "    - Estimated remaining:" << int((endTime - prevTime) * (NUM_TIMESTEPS - i)/NUM_TIME_SUBSTEPS) << "s\n";
			prevTime=endTime;
		}
		for(unsigned int j = 0;j<NUM_TIME_SUBSTEPS;j+=1) //0222
		{
			//fancy progress printing
			if(rank==0 && (i*NUM_TIME_SUBSTEPS+j)%FREQ_TERMINAL==0)
			{
				cout << " - Timestep ";
				cout.width(7);//increase this to give more room for large timestep numbers
				cout<<i+j<<"/"<<NUM_TIMESTEPS<<" (";
				cout.width(3);
				cout<<100*(i+j)/NUM_TIMESTEPS<<"%) : t=";
				cout.width(5);
				cout<<t;
				cout<<"\n";
			}
			//prevent time-consuming copying by switching pointers
			swapGrid=oldGrid;
			oldGrid=newGrid;
			newGrid=swapGrid;
			//newGrid->forwardEuler(DT,oldGrid,logFile);
			newGrid->RK33(DT,oldGrid,&k1,&k2,logFile);
			t += DT;
		}
		if (BOOL_FILE) dumpToHDF5File(outputFolder, newGrid, i+NUM_TIME_SUBSTEPS,t,rank);

		
		
	}

	if(rank==0)cout << "Post-processing\n";
	logFile.close();
	endTime = MPI_Wtime();
	if(rank==0)cout << "Exiting\nTotal runtime is "<<(endTime-beginTime)<<"s\n";

	MPI_Finalize();
	return ERRCODE_SUCESSFUL_COMPLETION;
}

#endif
