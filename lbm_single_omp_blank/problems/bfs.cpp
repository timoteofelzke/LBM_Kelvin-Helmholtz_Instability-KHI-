/**
 *  @file   bfs.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of functions used in the Backward Facing Step problem
 *
 *  Implementation of functions used to run and analise the Backward Facing Step (BFS) problem
 *
 *  @image html bfs.png
 *
 */


#include <fstream>
#include "../lbmbase/boundary/zouhe.h"
#include "../lbmbase/boundary/neumann.h"
#include "../lbmbase/geometry.h"
#include "../parameters.h"
#include "../lbmbase/lattice/lattice.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "bfs.h"

using namespace std;

namespace bfsProblem
{

stProblem problemParameters( string filename )
{
    stProblem bfs;
	using boost::property_tree::ptree;
	ptree pt;
	read_ini(filename.c_str(),pt);
    bfs.U   =  pt.get<double>("BFS.uavg");
    bfs.h1  =  pt.get<int>("BFS.h1");
    bfs.h2  =  pt.get<int>("BFS.h2");
    bfs.l1  =  pt.get<int>("BFS.l1");
    bfs.l2  =  pt.get<int>("BFS.l2");
    bfs.tavg = pt.get<int>("BFS.tavg");
    bfs.tstb = pt.get<int>("BFS.tstb");
    bfs.filename = pt.get<string>("BFS.filename","bfs.csv");

    return bfs;
}

void problemInitialize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force)
{

}

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& prm, stProblem& bfs, stForce& force, int step )
{
	static stBoundaryNeumann boundaryOutput = defineNeumannYZ(geo,geo.nx-1,-1);

    double y0 = bfs.h1 + 0.5*(bfs.h2+1.0 );

    int x = 0;
    for (int z = 0; z < geo.nz; z++)
    {
        for (int y = bfs.h1 + 1; y < ( geo.ny -1 ) ; y++)
        {
            int id = getIndex(geo,x,y,z);
            if (id != 0)
            {
                double* f = iniN + mapMemory(id-1);
                double v = bfs.U * parabolicProfile( (y - y0)/static_cast<double>(bfs.h2) );
                zouHeVelocityLeft(f,v,0,0);
            }
        }

    }

    neumann( boundaryOutput, iniN, prm.initialRho);
}

void problemPosCollision(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step )
{

}

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& prm, stProblem& bfs, stForce& force, int step )
{
	int y = 1;
	int z = 0;

    static double* vmed = new double[geo.nx - bfs.l1] ();   // Vector to store the time average of the velocity at bottom of the channel after the step
	static int t = 0;

    if ( (t < bfs.tavg) && (step >= bfs.tstb) )
	{
		t = t + 1;
        for (int x= bfs.l1; x < geo.nx ; x++)
		{
			double rho,vx,vy,vz;
			int id = getIndex(geo,x,y,z);
            double* f = iniN + mapMemory(id-1);
			calculateMacro(f,vx,vy,vz,rho);
            vmed[x - bfs.l1 ] += vx;
		}
	}

	if (t == bfs.tavg)
	{
        ofstream file( (bfs.filename).c_str() );
        for (int x= bfs.l1; x < geo.nx ; x++)
		{
            file << 0.5 + x - bfs.l1 << "\t" << vmed[x - bfs.l1] / float(t) << endl;
            vmed[x - bfs.l1] = 0.0;
		}
		file.close();
		t = 0;
	}
	return true;
}

double parabolicProfile(double _x)
{
    return -6.0*(_x-0.5)*(_x+0.5);
}

}
