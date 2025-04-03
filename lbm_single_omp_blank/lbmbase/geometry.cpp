/**
 *  @file    geometry.cpp
 *  @author  Diogo Nardelli Siebert
 *  @brief   Implementation of the functions to read and manage geometry information.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "vtkio/vtkio.h"
#include "geometry.h"

using namespace std;

stGeometry readGeometry(std::string filename ,std::string filetype, int sizeX, int sizeY, int sizeZ)
{
        if ( filetype == "raw" )
        {
                return readGeometryRAW( filename, sizeX, sizeY, sizeZ);
        }
        else if ( filetype == "vtk" )
        {
                return readGeometryVTK( filename );
        }
	else if ( filetype == "nofile" )
	{
		return createFullGeometry( sizeX, sizeY, sizeZ);
	}
}

void reportGeometry(stGeometry& geo)
{
    stringstream dimensions;
    dimensions << "[" << geo.nx << " "   << geo.ny << " "   << geo.nz << "]";

    cout << endl;
    cout << left << setw(40) << "Dimensions: ";
    cout << setw(30) << right << dimensions.str() << endl;

    cout << left << setw(40) << "Fraction of fluid nodes: ";
    cout << setw(30) << right << static_cast<double>(geo.numberOfPoints)/ static_cast<double>(geo.nx*geo.ny*geo.nz) << endl << endl;
}

void updateGeometry(stGeometry& geo)
{
	for (int pos = 0; pos < geo.nx*geo.ny*geo.nz; pos++)
	{
		int id = geo.index[pos];
		if (id)
		{
			geo.position[id-1] = pos;
		}
	}
}

unsigned int getPosition(stGeometry &geo,int id)
{
    return geo.position[id-1];
}

void getPosition(stGeometry &geo,int id, int& x,int &y,int &z)
{
    unsigned int pos = getPosition(geo,id);
	z =   pos / (geo.nx*geo.ny);
	unsigned int mod = pos % (geo.nx*geo.ny);
	y = ( mod ) / geo.nx;
	x = ( mod ) % geo.nx;
}

unsigned int getIndex(stGeometry &geo,int x,int y, int z)
{
	int xNew = ( x + geo.nx ) % geo.nx;
	int yNew = ( y + geo.ny ) % geo.ny;
	int zNew = ( z + geo.nz ) % geo.nz;
	return geo.index[ xNew + yNew * geo.nx + zNew * geo.nx * geo.ny ];
}

stGeometry createFullGeometry(int nx, int ny, int nz)
{
	stGeometry geo;
        geo.nx = nx;
        geo.ny = ny;
        geo.nz = nz;

        geo.index = new unsigned int[geo.nx*geo.ny*geo.nz];
        geo.numberOfPoints = geo.nx*geo.ny*geo.nz;

	for (int pos = 0; pos < geo.numberOfPoints; pos++) geo.index[pos] = pos+1;

        geo.position = new unsigned int[geo.numberOfPoints];
        updateGeometry(geo);
	return geo;
}

stGeometry readGeometryRAW(string geoFilename,int nx,int ny,int nz)
{
	stGeometry geo;
	ifstream geoFile( geoFilename.c_str() );

	geo.nx = nx;
	geo.ny = ny;
	geo.nz = nz;

	geo.index = new unsigned int[geo.nx*geo.ny*geo.nz];
	geo.numberOfPoints = 0;

	// Reads from the file if the pixel is solid or fluid
	// if the pixel is a solid stores 0 in the pGeometry array
	// if the pixel is a fluid it stores the fluid number starting from one

	unsigned char* rawData = new unsigned char[geo.nx*geo.ny*geo.nz];

	geoFile.read( (char*) rawData,geo.nx*geo.ny*geo.nz);

	for ( int z = 0; z < geo.nz; z++ )
	{
		for ( int y = 0; y < geo.ny; y++ )
		{
			for ( int x = 0; x < geo.nx; x++ )
			{
				int pos = x + y * geo.nx + z * geo.ny * geo.nx;
				geo.index[pos] = 0;
				int siteType = rawData[pos];
				if (  siteType == FLUID )
				{
					geo.numberOfPoints ++;
					geo.index[pos] = geo.numberOfPoints;
				}
			}
		}
	}

	geo.position = new unsigned int[geo.numberOfPoints];
	updateGeometry(geo);

	geoFile.close();
	
	delete rawData;
	return geo;
}

stGeometry readGeometryVTK(string geoFilename)
{
	stGeometry geo;

	// Open file containing the geometry
	// Dump all the VTK header with the exception of geometry size
	ifstream geoFile( geoFilename.c_str() );

	stVtkHeader header = readVtkHeader( geoFile );	
	geo.nx = header.sizeX;
	geo.ny = header.sizeY;
	geo.nz = header.sizeZ;

	geo.index = new unsigned int[geo.nx*geo.ny*geo.nz];
	geo.numberOfPoints = 0;

	// Reads from the file if the pixel is solid or fluid
	// if the pixel is a solid stores 0 in the pGeometry array
	// if the pixel is a fluid it stores the fluid number starting from one

	for (int pos = 0; pos < geo.nx * geo.ny * geo.nz ; pos ++)
	{
		geo.index[pos] = 0;
		int siteType = 2 ;
		if   (header.variableType == "char")     siteType = readVtkChar(geoFile,header);
		else if (header.variableType == "int")   siteType = readVtkInt (geoFile,header);
		if (  siteType == FLUID )
		{
			geo.numberOfPoints ++;
			geo.index[pos] = geo.numberOfPoints;
		}		
	}

	geo.position = new unsigned int[geo.numberOfPoints];
	updateGeometry(geo);

	geoFile.close();
	return geo;
}
