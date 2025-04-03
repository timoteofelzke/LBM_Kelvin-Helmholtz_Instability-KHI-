/**
 *  @file   streaming.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of functions necesssary for the streaming step.
 */

#include "streaming.h"
#include "../lattice/lattice.h"
#include <algorithm>

using namespace std;

void preStreaming(double* N, int numberOfPoints)
{
	#pragma omp parallel for
	for (int n =0 ; n < numberOfPoints; n++)
	{
		// Perform the swap between opposite directions in the same site.
		// The lattice vectors must be organize in such way that i and i+1 are opposite directions.
		for (int i = 1; i < NUM_OF_VEL; i = i + 2)
		{
			swap( N[ n * NUM_OF_VEL + i] ,  N[ n * NUM_OF_VEL + i + 1 ]  );
		}
	}
}

void streaming(stStreaming &streamInfo, double* N)
{
	#pragma omp parallel for schedule(static,10000)
	for (unsigned int k=0; k< streamInfo.cnt  ; k = k + 2)
	{
		std::swap( N[ streamInfo.buffer[k+1] ]   , N[  streamInfo.buffer[k] ] ) ;
	}
}

stStreaming defineStreaming(stGeometry& geo)
{
	stStreaming streamInfo;
    streamInfo.buffer =  new unsigned int[ geo.numberOfPoints *  (NUM_OF_VEL-1)] () ;
	streamInfo.cnt = 0;

	for (int idOrigin = 1 ; idOrigin <= geo.numberOfPoints; idOrigin++)
	{
		int x,y,z;
		getPosition(geo,idOrigin,x,y,z);

		for ( int i = 1; i < NUM_OF_VEL; i++)
		{
			// Computes x,y,z position of the target site
			int j = i_op[i];



			// Stores the position in the distribution function array where distribution for the direction i
			// of the local site is stored.

			int idTarget = getIndex( geo, x + cx[j] ,y + cy[j], z + cz[j] );

            unsigned int memoryOrigin = static_cast<unsigned int>(idOrigin - 1 )*NUM_OF_VEL + i;
            unsigned int memoryTarget = static_cast<unsigned int>(idTarget - 1 )*NUM_OF_VEL + j;

			// If the target is a solid the nothing is done.
			// Otherwise, if the target is a fluid site, record the position in the array
			// for the direction j (opposite to i) of the target site.

			if ( ( idTarget ) && ( memoryOrigin < memoryTarget) )
			{
				streamInfo.buffer[ streamInfo.cnt  ] = memoryOrigin;
				streamInfo.buffer[ streamInfo.cnt+1] = memoryTarget;
				streamInfo.cnt += 2;
			}
		}
	}

	return streamInfo;
}
