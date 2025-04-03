 /**
 *  @file   lbm.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Main file of the Lattice Boltzmann Method code.
 */

#include <iostream>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <omp.h>
#include <cmath>
#include <iomanip>

#include <boost/timer/timer.hpp>

#include "lbmbase/lattice/lattice.h"
#include "lbmbase/streaming/streaming.h"
#include "lbmbase/collision/collision.h"
#include "lbmbase/force/force.h"
#include "lbmbase/initial.h"
#include "lbmbase/dataout.h"
#include "lbmbase/geometry.h"

#include "options.h"
#include "parameters.h"

#include "problems/problems.h" // modelo de problema a ser resolvido, ex: cavity.h (altere o arquivo problems.h para incluir o modelo de problema desejado)

using namespace std;
using namespace boost::timer;

/* Set namespace with specific functions to be executed for a particular problem */ 

void divergeCheck(double* ini_N,stGeometry geo);

int main(int argc, char** argv)
{
	printInfo();
    /* Read data from files and command line */

	stOptions opt = readCommandLine(argc,argv);
	stParameters prm = readDataINI(opt.filename);
	stGeometry geo = readGeometry( prm.geoFilename, prm.geoFiletype, prm.geoSizeX, prm.geoSizeY, prm.geoSizeZ );
	stProblem bfs = problemParameters(opt.filename);
        stForce force = setConstantForce(prm.accX,prm.accY,prm.accZ);
        stCollision col = collisionParameters(opt.filename);

    /* Print parameters and options */
	reportOptions(opt);
	reportParameters( prm );
	reportGeometry(geo);
        reportCollision(col);

    /* Create an tStreaming with information necessary for the
     * streaming step, for more information see streaming.h  */
	stStreaming swapBuffer = defineStreaming(geo);

	unsigned int initialStep = 0;
	unsigned int finalStep = prm.numberOfSteps;

	double *ini_N;

    /* Allocate memory and set initial the conditions */
	if (prm.setInitial) ini_N = setInitialFromVtk(geo, prm.iniRhoFilename, prm.iniVelFilename);
	else ini_N = setInitial(geo, prm.initialRho);

    /* Load the iniN with data backup from another process */
        if (opt.continueFlag)
	{
		loadRecovery("sim.bak",ini_N,geo.numberOfPoints,initialStep);
	}

    /* Create timers to mesure how much each step of the algorithm take to run */
	cpu_timer simulationTimer;
	cpu_timer collisionTimer , streamingTimer, boundaryTimer;

	collisionTimer.stop();
	streamingTimer.stop();
	boundaryTimer.stop();

	double nextBackup = opt.backupInterval * 60.;

    problemInitialize(ini_N,col,geo,prm,bfs,force);

	if  (prm.writeZero)
	{
		if( prm.fileDensity)  vtkDensity( geo,  ini_N, 0, prm.fileFormat);
        if (prm.fileVelocity) vtkVelocity( geo, ini_N, force , 0 , prm.fileFormat);
	}

	unsigned int step = initialStep;
	bool keepRunning = true;

	cout << "Step: " << flush;
	while ( (step < finalStep) && (keepRunning) )
	{
		boundaryTimer.resume();
	        problemBoundary(ini_N,col,geo,prm,bfs,force,step);
	
        	if  ( ( (step) % prm.recordInterval ==  0 ) && ( (step==0)?prm.writeZero:true ) )
	        {
        	        if( prm.fileDensity)  vtkDensity(  geo,  ini_N,  step, prm.fileFormat);
	                if (prm.fileVelocity) vtkVelocity( geo, ini_N,  force ,step, prm.fileFormat);
        	}

		boundaryTimer.stop();
		collisionTimer.resume();
	        collisionSwap(ini_N, force, col , geo.numberOfPoints);
		collisionTimer.stop();

        	problemPosCollision(ini_N,col,geo,prm,bfs,force,step);

		streamingTimer.resume();
		streaming(swapBuffer,ini_N);
		streamingTimer.stop();

		if ( (step+1) % ( 100 ) == 0 )
		{
	            cout << right << setw(12) << (step+1) << flush << " ";
        	    if ( (step+1)%500 == 0) cout << endl << setw(6) << " ";
                    divergeCheck(ini_N,geo);
		}

		double elapsedTime = static_cast<double>( simulationTimer.elapsed().wall ) *1E-9;

		if  ( (opt.backupFlag) && ( elapsedTime > nextBackup ) )
		{
			nextBackup  = nextBackup + opt.backupInterval * 60.;
			saveRecovery("sim.bak",ini_N,geo.numberOfPoints,step+1);
		}

	        keepRunning = problemOutput(ini_N,col,geo,prm,bfs,force,step+1);
		step++;
	}

    problemBoundary(ini_N,col,geo,prm,bfs,force,step);

    if  (prm.writeFinal)
    {
        if( prm.fileDensity)  vtkDensity( geo,  ini_N, step, prm.fileFormat);
        if (prm.fileVelocity) vtkVelocity( geo, ini_N,  force ,step, prm.fileFormat);
    }

	cout << endl << endl;
	cout << left << setw(40) << "Collision time: " << right << setw(30) << collisionTimer.format(4,"%ws") <<  endl;
	cout << left << setw(40) << "Streaming time: " << right << setw(30) <<  streamingTimer.format(4,"%ws") <<  endl;
	cout << left << setw(40) << "Boundary time: "  << right << setw(30) <<  boundaryTimer.format(4,"%ws") <<  endl << endl;
	cout << left << setw(40) <<  "Total time: " <<  right << setw(30) << simulationTimer.format(4,"%ws") <<  endl;

        double elapsedTime = static_cast<double>( simulationTimer.elapsed().wall ) *1E-9;
        cout << left << setw(40) <<  "Performance: " <<  right << setw(30) << (finalStep-initialStep) * static_cast<double>( geo.numberOfPoints ) / ( simulationTimer.elapsed().wall * 1E-3 ) << " MLPs" << endl;

	problemFinalize(ini_N,col,geo,prm,bfs,force,step+1);


	if  (opt.endBackupFlag) saveRecovery("sim.bak",ini_N,geo.numberOfPoints,finalStep);

}

void divergeCheck(double* ini_N,stGeometry geo)
{    
    for (int n = 0; n < geo.numberOfPoints; n += (geo.numberOfPoints/10 > 0 ? geo.numberOfPoints/10: 1) )
	{
		double rho,vx,vy,vz;
		double *f = ini_N + n*NUM_OF_VEL;
		calculateMacro(f,vx,vy,vz,rho);
		if (rho>10)
		{
			cout << "Simulations has diverged!!!" << endl;
			exit(EXIT_FAILURE);
		}
	}
}

