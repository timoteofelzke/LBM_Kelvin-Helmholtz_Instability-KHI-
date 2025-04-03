/**
 *  @file   parameters.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of functions needed to read the input data from INI files.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

#include "parameters.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;

stParameters readDataINI( std::string filename)
{
	stParameters prm;

	using boost::property_tree::ptree;
	ptree pt;
	read_ini(filename.c_str(),pt);

//	prm.tau =  pt.get<double>("General.tau");
	prm.numberOfSteps = pt.get<int>("General.number_of_steps");
	prm.initialRho = pt.get("General.initial_density",1.f);
	prm.geoFiletype = pt.get<string>("Geometry.filetype");

	if ( (prm.geoFiletype == "raw") || (prm.geoFiletype == "nofile") )
	{
		prm.geoSizeX = pt.get<int>("Geometry.size_x");
		prm.geoSizeY = pt.get<int>("Geometry.size_y");
		prm.geoSizeZ = pt.get<int>("Geometry.size_z");
	}

        if ( (prm.geoFiletype == "raw") || (prm.geoFiletype == "vtk") ) prm.geoFilename = pt.get<string>("Geometry.filename");
	prm.accX = pt.get<double>("Acceleration.x",0.f);
	prm.accY = pt.get<double>("Acceleration.y",0.f);
	prm.accZ = pt.get<double>("Acceleration.z",0.f);

	prm.writeZero  = pt.get<bool>("Output.write_zero_step",false);
	prm.writeFinal = pt.get<bool>("Output.write_final_step",false);

	prm.numberOfFiles = pt.get<int>("Output.number_of_files");
	prm.fileVelocity = pt.get<bool>("Output.velocity");
	prm.fileDensity = pt.get<bool>("Output.density");

	prm.fileFormat = 0;

	if ( pt.get<bool>("Output.binary",true ) ) prm.fileFormat += 1;

        if ( pt.get<bool>("Output.xml",false) )
        {
             prm.fileFormat += 2;
             if (pt.get<bool>("Output.compress",false) ) prm.fileFormat = 4;
 	}

	prm.recordInterval =   prm.numberOfSteps + 1;
	if (prm.numberOfFiles != 0) prm.recordInterval =  prm.numberOfSteps / prm.numberOfFiles;
	prm.viscosity = 1.0 / 3.0 * ( prm.tau - 0.5 );

	prm.setInitial = pt.get<bool>("Initial_Condition.use_files",false);
	prm.iniRhoFilename = pt.get<string>("Initial_Condition.density_file","");
	prm.iniVelFilename = pt.get<string>("Initial_Condition.velocity_file","");

	return prm;
}

void reportParameters( stParameters prm)
{
    cout << endl;
    cout << left << setw(40) << "Geometry file:";
    cout << right << setw(30) << prm.geoFilename << endl;

    cout << left << setw(40) << "Number of steps:";
    cout << right << setw(30) << prm.numberOfSteps << endl;

    cout << left << setw(40) << "Number of recorded VTK files:";
    cout << right << setw(30) << prm.numberOfFiles << endl;

    cout << left << setw(40) << "Relaxation time:";
    cout << right << setw(30)<< prm.tau << endl;

    cout << left << setw(40) << "Initial density:";
    cout << right << setw(30) << prm.initialRho << endl;

    stringstream acceleration;
    acceleration << "[" << prm.accX << " " << prm.accY << " "   << prm.accZ << "]";
    cout << left << setw(40) << "Acceleration:";
    cout << setw(30) << right << acceleration.str() << endl;
	
	if (prm.setInitial)
	{
		cout << "Reading initial conditions from files : " << prm.iniRhoFilename << " and " << prm.iniVelFilename << endl;
	}
	
}
