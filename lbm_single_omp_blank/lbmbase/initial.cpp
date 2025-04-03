/**
 *  @file   initial.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of the functions to read and apply initial conditions
 */

#include <string>
#include <fstream>

#include "geometry.h"
#include "lattice/lattice.h"
#include "initial.h"
#include "vtkio/vtkio.h"

using namespace std;

double * setInitial(stGeometry& geo, double rho)
{
    double *iniN = new double[ geo.numberOfPoints *  (NUM_OF_VEL  )];
    for ( int n  = 0; n  < geo.numberOfPoints; n++ )
	{
        double *f = iniN + mapMemory( n );
        equilibriumDistribution( 0.0 , 0.0 , 0.0 , rho ,  f);
	}
    return iniN;
}

double * setInitialFromVtk(stGeometry& geo,string rhoFilename, string velFilename)
{
    double *iniN = new double[ geo.numberOfPoints * NUM_OF_VEL ]; // Allocate memory for the distribution function
	
	ifstream rhoVtkFile( rhoFilename.c_str() , ifstream::binary );  // Open vtkfile
	stVtkHeader rhoHeader = readVtkHeader( rhoVtkFile );											// Read header of the vtkfile
	
	ifstream velVtkFile;
	if (velFilename != "") velVtkFile.open( velFilename.c_str() , ifstream::binary );  // Open vtkfile
	
	stVtkHeader velHeader = readVtkHeader( velVtkFile );
	
	for ( int z = 0; z < geo.nz; z++ )
	{
		for ( int y = 0; y < geo.ny; y++ )
		{
			for ( int x = 0; x < geo.nx; x++ )
			{				
				double vx = 0, vy = 0, vz = 0;
				double rho = readVtkFloat(rhoVtkFile, rhoHeader) ;
				
				if (velFilename != "")
				{
					vx  = readVtkFloat(velVtkFile, velHeader) ;
					vy  = readVtkFloat(velVtkFile, velHeader) ;
					vz  = readVtkFloat(velVtkFile, velHeader) ;
				}
				
				int id = getIndex(geo,x,y,z);

				if (id != 0 )
				{
                    double *f = iniN + mapMemory( id - 1 );
                    equilibriumDistribution( vx , vy , vz , rho ,  f);
				}
			}
		}
	}	
	
    return iniN;
}
