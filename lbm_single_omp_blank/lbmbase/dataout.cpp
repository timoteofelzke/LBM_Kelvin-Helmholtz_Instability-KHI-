/**
 *  @file    dataout.cpp
 *  @author  Diogo Nardelli Siebert
 *  @brief   Implementation of functions to write output files (VTK) of the flow fields (density and velocity)
 */

#include <sstream>
#include <iostream>
#include <fstream>
#include "dataout.h"
#include "lattice/lattice.h"
#include "force/force.h"
#include "geometry.h"

using namespace std;

string intToString(int number)
{
	const int NUMBER_OF_DIGITS = 9;    // The number of digits used in the string to write the int (if the int is to small zeros will be added in the left)
	stringstream ss;    // Stringstream to convert the int to string
	ss << number;
	string converted = ss.str();    //  Convert the stringstream to string
	converted.insert(0,NUMBER_OF_DIGITS-converted.size(),'0');   // Add the necessary number of zeros
	return converted;
}

void vtkVelocity(stGeometry& geo, double *ini_R, double *ini_B, unsigned int step, int fileFormat)
{
    // Set the filename and open the output file
    string filename = "vel_" + intToString(step) + ( (fileFormat<2) ? ".vtk" : ".vti" ) ;

    stVtkHeader header = setVtkHeader("Velocidade",fileFormat, "VECTORS","velocidade","float", geo.nx ,geo.ny , geo.nz);

    ofstream fvel(filename.c_str() ,ios::binary);
    writeVtkHeader( fvel, header) ;

    for ( int z = 0; z < geo.nz; z++ )
    {
        for ( int y = 0; y < geo.ny; y++ )
        {
            for ( int x = 0; x < geo.nx; x++ )
            {
                unsigned int id = getIndex(geo,x,y,z);
                double vx, vy, vz, rho;
                vx = vy = vz = 0;

                if ( id )
                {
                    double *R = ini_R + ( id  - 1 ) * NUM_OF_VEL;
                    double *B = ini_B + ( id  - 1 ) * NUM_OF_VEL;
                    double f[NUM_OF_VEL];

                    for (int i = 0; i < NUM_OF_VEL ; i++) f[i] = R[i] + B[i];

                    calculateMacro( f, vx, vy, vz, rho );
                }

                writeVtkFloat( fvel, header, vx);
                writeVtkFloat( fvel, header, vy);
                writeVtkFloat( fvel, header, vz);

            }
        }
    }

    writeVtkFooter( fvel, header) ;
    fvel.close();
}


void vtkDensity( stGeometry& geo, double *ini_R, double* ini_B, unsigned int step, int fileFormat)
{
    // Set the filename and open the output file
    string filename = "rho_" + intToString(step) + ( (fileFormat<2) ? ".vtk" : ".vti" ) ;
    stVtkHeader header = setVtkHeader("Density",fileFormat,"SCALARS","concentration","float", geo.nx ,geo.ny , geo.nz);

    // Set the filename and open the output file
    ofstream frho( filename.c_str() );

    // Write the vtk file header (http://www.vtk.org/VTK/img/file-formats.pdf)
    writeVtkHeader( frho, header) ;

    //Write the mass density field, if the point is a solid a zero density is written.
    for ( int z = 0; z < geo.nz; z++ )
    {
        for ( int y = 0; y < geo.ny; y++ )
        {
            for ( int x = 0; x < geo.nx; x++ )
            {
                unsigned int id = getIndex(geo,x,y,z);

                float rho = 0;

                if ( id )
                {
                    double *R = ini_R + ( id  - 1 ) * NUM_OF_VEL;
                    double *B = ini_B + ( id  - 1 ) * NUM_OF_VEL;

                    float rhoR = calculateMass(R);
                    float rhoB = calculateMass(B);
                    rho = rhoR + rhoB;
                }
                writeVtkFloat( frho , header, rho) ;
            }
            //frho << endl;
        }

    }

    writeVtkFooter( frho, header) ;
    frho.close();
}

void vtkConcentration( stGeometry& geo, double *ini_R, double* ini_B, unsigned int step, int fileFormat)
{
    // Set the filename and open the output file
    string filename = "wr_" + intToString(step) + ( (fileFormat<2) ? ".vtk" : ".vti" ) ;
    stVtkHeader header = setVtkHeader("Concentration",fileFormat,"SCALARS","concentration","float", geo.nx ,geo.ny , geo.nz);

    // Set the filename and open the output file
    ofstream frho( filename.c_str() );

    // Write the vtk file header (http://www.vtk.org/VTK/img/file-formats.pdf)
    writeVtkHeader( frho, header) ;

    //Write the mass density field, if the point is a solid a zero density is written.
    for ( int z = 0; z < geo.nz; z++ )
    {
        for ( int y = 0; y < geo.ny; y++ )
        {
            for ( int x = 0; x < geo.nx; x++ )
            {
                unsigned int id = getIndex(geo,x,y,z);

                float wR = 0;

                if ( id )
                {
                    double *R = ini_R + ( id  - 1 ) * NUM_OF_VEL;
                    double *B = ini_B + ( id  - 1 ) * NUM_OF_VEL;

                    float rhoR = calculateMass(R);
                    float rhoB = calculateMass(B);
                    wR = rhoR/(rhoR + rhoB);
                }
                writeVtkFloat( frho , header, wR) ;
            }
            //frho << endl;
        }

    }

    writeVtkFooter( frho, header) ;
    frho.close();

}

void vtkVelocity(stGeometry& geo, double *ini_N, stForce& force, unsigned int step, int fileFormat)
{
	// Set the filename and open the output file
    string filename = "vel_" + intToString(step) + ( (fileFormat<2) ? ".vtk" : ".vti" ) ;
	stVtkHeader header;

    header = setVtkHeader("Velocidade",fileFormat,"VECTORS","velocidade","float", geo.nx ,geo.ny , geo.nz);
	
	ofstream fvel(filename.c_str() ,ios::binary);
    writeVtkHeader( fvel, header) ;

	for ( int z = 0; z < geo.nz; z++ )
	{
		for ( int y = 0; y < geo.ny; y++ )
		{
			for ( int x = 0; x < geo.nx; x++ )
			{
				unsigned int id = getIndex(geo,x,y,z);
				double vx, vy, vz, rho;
				vx = vy = vz = 0;

				if ( id )
				{
					double *f = ini_N + ( id  - 1 ) * NUM_OF_VEL;
                    double* F = force.constant ? force.constantForce : force.forceField + ( id  - 1 ) * NUM_OF_DIM;
                    calculateMacroOutput( f, vx, vy, vz, rho , F[0], F[1], F[2] ) ;
				}

				writeVtkFloat( fvel, header, vx);
				writeVtkFloat( fvel, header, vy);
				writeVtkFloat( fvel, header, vz);

			}
		}
	}

    writeVtkFooter( fvel, header) ;
	fvel.close();
}

void vtkDensity( stGeometry& geo, double *ini_N, unsigned int step, int fileFormat, string prefix )
{
    // Set the filename and open the output file
    string filename = prefix + "_" + intToString(step) + ( (fileFormat<2) ? ".vtk" : ".vti" ) ;
    stVtkHeader header;

    header = setVtkHeader("Densidade",fileFormat,"SCALARS","densidade","float", geo.nx ,geo.ny , geo.nz);

    // Set the filename and open the output file
    ofstream frho( filename.c_str() );

    // Write the vtk file header (http://www.vtk.org/VTK/img/file-formats.pdf)
    writeVtkHeader( frho, header) ;

    //Write the mass density field, if the point is a solid a zero density is written.
    for ( int z = 0; z < geo.nz; z++ )
    {
        for ( int y = 0; y < geo.ny; y++ )
        {
            for ( int x = 0; x < geo.nx; x++ )
            {
                unsigned int id = getIndex(geo,x,y,z);

                float rho = 0;

                if ( id )
                {
                    double *f = ini_N + ( id - 1 ) * NUM_OF_VEL;
                    rho = calculateMass(f);
                }
                writeVtkFloat( frho , header, rho);
            }
            //frho << endl;
        }

    }

    writeVtkFooter( frho, header) ;
    frho.close();
}

void saveRecovery(const char* filename, double *pDistribution, int numberOfPoints,unsigned int step)
{
	cout << endl << "Recording Backup: .... " << flush;
	ofstream output(filename);
	output.write( (char*) &step , sizeof(int));
	output.write( (char*) pDistribution , sizeof(double)* numberOfPoints * NUM_OF_VEL);
	output.close();
	cout << "Done!" << endl;
}

void loadRecovery(const char* filename, double *pDistribution, int numberOfPoints, unsigned int& step)
{
	cout << endl << "Reading Backup: .... " << flush;
	ifstream input(filename);
	input.read( (char*) &step , sizeof(int));
	input.read( (char*) pDistribution , sizeof(double)* numberOfPoints * NUM_OF_VEL);
	input.close();
	cout << "Done! " << endl;
}
