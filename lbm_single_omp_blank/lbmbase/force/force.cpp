#include "force.h"
#include "../geometry.h"

stForce setConstantForce(double& forceX, double& forceY, double& forceZ )
{
    stForce force;

    force.constant = true;
    force.constantForce[0] = forceX;
    force.constantForce[1] = forceY;
    force.constantForce[2] = forceZ;

    return force;
}

stForce setVariableForcePump(stGeometry& geo,  double& forceX, double& forceY, double& forceZ, int axis, int p0, int p1)
{
    stForce force;
    force.constant = false;
    force.forceField = new double[ 3* geo.numberOfPoints ];

    int r[3];
    for (int id = 1; id <= geo.numberOfPoints ; id++)
    {
        getPosition( geo,id, r[0], r[1], r[2]);
        if ( (r[axis] >= p0) && (r[axis] < p1) )
        {
            force.forceField[3*(id-1) + 0] = forceX;
            force.forceField[3*(id-1) + 1] = forceY;
            force.forceField[3*(id-1) + 2] = forceZ;
        }
    }

    return force;
}

