/**
 *  @file   drag.cpp
 *  @author Damylle Cristina Xavier Donati
 *  @brief  Implementation of functions to compute the drag around an object
 *
 *  @image  html drag.png
 */

#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cmath>
#include "../lbmbase/boundary/zouhe.h"
#include "../lbmbase/boundary/neumann.h"
#include "../lbmbase/boundary/convective.h"
#include "../lbmbase/geometry.h"
#include "../parameters.h"
#include "../lbmbase/lattice/lattice.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "drag.h"

using namespace std;

namespace dragProblem
{

stProblem problemParameters( string filename )
{

	stProblem dg;

	using namespace boost;
	using namespace property_tree;

	ptree pt;

	read_ini(filename.c_str(),pt);

	try
	{
		dg.U         	   =  pt.get<double>("Drag.U");
		dg.computeInterval =  pt.get<int>("Drag.compute_interval",100);
		dg.filename  	   =  pt.get<string>("Drag.filename","drag.csv");
		dg.precision       =  pt.get<double>("Drag.precision",0);
	}
	catch(const ptree_error &e)
	{
              cerr << "Error in the INI File: ";
              cerr << e.what() << endl;
	      exit(3);
	}

	return dg;
}

void problemInitialize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
{

}

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& prm, stProblem& dg, stForce& force, int step )
{
    DragVelYBeg(iniN, geo, dg.U);
    DragVelYEnd(iniN, geo, dg.U);
    DragVelXBeg(iniN, geo, dg.U);
    DragVelZBeg(iniN, geo, dg.U);
    DragVelZEnd(iniN, geo, dg.U);
    static stBoundaryConvective outlet = defineConvectiveYZ(geo, geo.nx-1, -1);
    convective( outlet , iniN);
}

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& prm, stProblem& dragData, stForce& force, int step )
{
    string filename = dragData.filename;
    bool keepRunning = true;

    if (step == 1)
    {
        setupDragComputation(geo,dragData);
        ofstream output(filename.c_str());
        output << "Step\t Fx \t Fy \t Fz \n";
        output.close();
    }

    double oldForce = dragData.dpX;

    if ( (step%dragData.computeInterval == 0 ) && (step>100) )
    {
        computeDrag(iniN,geo,dragData);
        ofstream output(filename.c_str(), ofstream::app );
        output << step << "\t" << dragData.dpX << "\t" << dragData.dpY << "\t"<< dragData.dpZ << endl;
        output.close();

        if (dragData.precision > 0 )
        {
            double error = abs( (oldForce - dragData.dpX)/dragData.dpX );
            if (error < dragData.precision ) keepRunning = false;
        }
    }

    return keepRunning;
}

void DragVelXBeg(double* iniN, stGeometry &geo,double u)
{
    int x = 0;

    #pragma omp parallel for
    for (int z=0; z < geo.nz; z++)
    {
	for (int y=0; y < geo.ny; y++)
	{
            int id = getIndex(geo,x, y, z);
	    if (id != 0)
	    {
        double* f = iniN + mapMemory(id-1);
		zouHeVelocityLeft( f , u, 0.0, 0.0);
	    }
	}
    }
}

void DragVelYBeg(double* iniN, stGeometry &geo,double u)
{
    int y = 0;

    #pragma omp parallel for
    for (int z=0; z < geo.nz; z++)
    {
	for (int x=0; x < geo.nx - 1; x++)
	{
            int id = getIndex(geo,x, y, z);
	   if (id != 0)
	   {
           double* f = iniN + mapMemory(id-1);
	       zouHeVelocityFront( f , u, 0.0, 0.0);
	   }
	}
    }
}

void DragVelYEnd(double* iniN, stGeometry &geo,double u)
{
    int y = geo.ny-1;

    #pragma omp parallel for
    for (int z=0; z < geo.nz; z++)
    {
	for (int x=0; x < geo.nx - 1; x++)
	{
           int id = getIndex(geo,x, y, z);
	   if (id != 0)
	   {
           double* f = iniN + mapMemory(id-1);
	       zouHeVelocityBack( f , u, 0.0, 0.0);
	   }
	}
    }
}

void DragVelZBeg(double* iniN, stGeometry &geo,double u)
{
    int z = 0;

    #pragma omp parallel for
    for (int y=0; y < geo.ny; y++)
    {
	for (int x=0; x < geo.nx - 1; x++)
	{
            int id = getIndex(geo,x, y, z);
	    if (id != 0)
	    {
        double* f = iniN + mapMemory(id-1);
		zouHeVelocityBottom( f , u, 0.0, 0.0);
	    }
	}
    }
}

void DragVelZEnd(double* iniN, stGeometry &geo,double u)
{
    int z = geo.nz-1;

    #pragma omp parallel for
    for (int y=0; y < geo.ny; y++)
    {
	for (int x=0; x < geo.nx - 1; x++)
	{
            int id = getIndex(geo,x, y, z);
	    if (id != 0)
	    {
        double* f = iniN + mapMemory(id-1);
		zouHeVelocityTop( f , u, 0.0, 0.0);
	    }
	}
    }
}

void setupDragComputation(stGeometry &geo, stProblem& dragData)
{
    for (int id = 1; id <= geo.numberOfPoints; id++)
    {
        int x, y ,z;
        getPosition(geo,id,x,y,z);
        for (int i=0; i < NUM_OF_VEL; i++)
        {
            int id_neigh = getIndex(geo, x + cx[i], y + cy[i], z + cz[i]);
            if (id_neigh == 0)
            {
                dragData.idList.push_back( id );
                dragData.iList. push_back( i  );
            }
        }
    }
}

void computeDrag(double* iniN, stGeometry &geo, stProblem& dragData)
{
    double dpX = 0;
    double dpY = 0;
    double dpZ = 0;

    #pragma omp parallel for reduction(+:dpX,dpY,dpZ)
    for (unsigned int n=0;n < dragData.idList.size();n++)
    {
                int i  = dragData.iList[n];
                int id = dragData.idList[n];
                double* f = iniN + mapMemory(id-1);
                dpX +=  f[i_op[i]] * cx[i];
                dpY +=  f[i_op[i]] * cy[i];
                dpZ +=  f[i_op[i]] * cz[i];
    }

    dragData.dpX = 2 * dpX;
    dragData.dpY = 2 * dpY;
    dragData.dpZ = 2 * dpZ;
}



}
