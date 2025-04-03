#include <vector>
#include "wetting.h"

#include "../lattice/lattice.h"
#include "../collision/collision.h"
#include "../geometry.h"

namespace shcTrt
{
    void wetting(stWetting& wet, double* N, stCollision& col)
    {
        #pragma omp parallel for schedule(static,10000)
        for (unsigned int k=0; k< wet.location.size()  ; k++)
        {
            N[ wet.location[k] ] = col.wrWall;
        }
    }
}

namespace shcMrt
{
    void wetting(stWetting& wet, double* N, stCollision& col)
    {
        #pragma omp parallel for schedule(static,10000)
        for (unsigned int k=0; k< wet.location.size()  ; k++)
        {
            N[ wet.location[k] ] = col.drWall;
        }
    }
}

stWetting defineWetting(stGeometry& geo)
{
    stWetting wet;

    for (int idOrigin = 1 ; idOrigin <= geo.numberOfPoints; idOrigin++)
    {
        int x,y,z;
        getPosition(geo,idOrigin,x,y,z);

        for ( int i = 1; i < NUM_OF_VEL; i++)
        {
            int j = i_op[i];
            int idTarget = getIndex( geo, x + cx[j] ,y + cy[j], z + cz[j] );
            if ( !idTarget ) wet.location.push_back( (idOrigin - 1 )*NUM_OF_VEL + i );
        }
    }

    return wet;
}
