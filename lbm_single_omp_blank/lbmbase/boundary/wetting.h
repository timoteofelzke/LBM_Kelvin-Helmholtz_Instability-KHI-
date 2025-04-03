#ifndef WETTING_H
#define WETTING_H

#include <vector>
#include "../collision/collision.h"

struct stGeometry;

struct stWetting
{
    std::vector<unsigned> location;
};

stWetting defineWetting(stGeometry& geo);

namespace shcTrt
{
    void wetting(stWetting& wet, double* N, stCollision& col);
}

namespace shcMrt
{
    void wetting(stWetting& wet, double* N, stCollision& col);
}

#endif
