/**
 *  @file   zouhe.cpp
 *  @author Diogo Nardelli Siebert
 *  @brief  Implementation of functions to impose Zou & He type boundary conditions.
 */

#include "zouhe.h"
#include "../geometry.h"
#include "../lattice/lattice.h"

namespace d3q15
{

void zouHePressureLeft(double *f, double rho)
{
    double Vx =(rho - f[0] - 2*f[2] - f[3] - f[4] - f[5] - f[6] - 2*f[8] - 2*f[9] - 2*f[12] - 2*f[14])/rho;
    double dx = (rho*Vx)/12.;
    double dy = (f[4]-f[3])/4.;
    double dz = (f[6]-f[5])/4.;

    f[1] = 8*dx + f[2];
    f[7] = dx + dy + dz + f[8];
    f[10] = dx - dy - dz + f[9];
    f[11] = dx - dy + dz + f[12];
    f[13] = dx + dy - dz + f[14];
}

void zouHePressureRight(double *f, double rho)
{
    double Vx = (-rho + f[0] + 2*f[1] + f[3] + f[4] + f[5] + f[6] + 2*f[7] + 2*f[10] + 2*f[11] + 2*f[13])/rho;
    double dx = (rho*Vx)/12.;
    double dy = (-f[3] + f[4])/4.;
    double dz = (-f[5] + f[6])/4.;

    f[2] = -8*dx + f[1];
    f[8] = -dx - dy - dz + f[7];
    f[9] = -dx + dy + dz + f[10];
    f[12] = -dx + dy - dz + f[11];
    f[14] = -dx - dy + dz + f[13];
}

void zouHePressureFront(double *f, double rho)
{
    double Vy =(rho - f[0] - f[1] - f[2] - 2*f[4] - f[5] - f[6] - 2*f[8] - 2*f[10] - 2*f[11] - 2*f[14])/rho;
    double dy = (rho*Vy)/12.;
    double dx = (-f[1] + f[2])/4.;
    double dz = (-f[5] + f[6])/4.;

    f[3] = 8*dy + f[4];
    f[7] = dx + dy + dz + f[8];
    f[9] = -dx + dy + dz + f[10];
    f[12] = -dx + dy - dz + f[11];
    f[13] = dx + dy - dz + f[14];
}


void zouHePressureBack(double *f, double rho)
{
    double Vy =(-rho + f[0] + f[1] + f[2] + 2*f[3] + f[5] + f[6] + 2*f[7] + 2*f[9] + 2*f[12] + 2*f[13])/rho;
    double dy = (rho*Vy)/12.;
    double dx = (-f[1] + f[2])/4.;
    double dz = (-f[5] + f[6])/4.;

    f[4] = -8*dy + f[3];
    f[8] = -dx - dy - dz + f[7];
    f[10] = dx - dy - dz + f[9];
    f[11] = dx - dy + dz + f[12];
    f[14] = -dx - dy + dz + f[13];
}

void zouHePressureBottom(double *f, double rho)
{
    double Vz = (rho - f[0] - f[1] - f[2] - f[3] - f[4] - 2*f[6] - 2*f[8] - 2*f[10] - 2*f[12] - 2*f[13])/rho;
    double dz = (rho*Vz)/12.;
    double dx = (-f[1] + f[2])/4.;
    double dy = (-f[3] + f[4])/4.;

    f[5] = 8*dz + f[6];
    f[7] = dx + dy + dz + f[8];
    f[9] = -dx + dy + dz + f[10];
    f[11] = dx - dy + dz + f[12];
    f[14] = -dx - dy + dz + f[13];
}

void zouHePressureTop(double *f, double rho)
{
    double Vz = (-rho + f[0] + f[1] + f[2] + f[3] + f[4] + 2*f[5] + 2*f[7] + 2*f[9] + 2*f[11] + 2*f[14])/rho;
    double dz = (rho*Vz)/12.;
    double dx = (-f[1] + f[2])/4.;
    double dy = (-f[3] + f[4])/4.;

    f[6] = -8*dz + f[5];
    f[8] = -dx - dy - dz + f[7];
    f[10] = dx - dy - dz + f[9];
    f[12] = -dx + dy - dz + f[11];
    f[13] = dx + dy - dz + f[14];
}

void zouHeVelocityBack(double* f, double Vx, double Vy, double Vz)
{
    double rho = (f[0] + f[1] + f[2] + f[5] + f[6] + 2 * (f[4] + f[8] + f[10] + f[11] + f[14])) / (1 - Vy);
    double dx = (rho * Vx - f[1] + f[2]) / 4.;
    double dz = (rho * Vz - f[5] + f[6]) / 4.;
    double dy = (rho * Vy) / 12.;

    f[3] = 8 * dy + f[4];
    f[7] = dx + dy + dz + f[8];
    f[9] = -dx + dy + dz + f[10];
    f[12] = -dx + dy - dz + f[11];
    f[13] = dx + dy - dz + f[14];
}

void zouHeVelocityFront(double* f, double Vx, double Vy, double Vz)
{
    double rho = (f[0] + f[1] + f[2] + f[5] + f[6] + +2 * (f[3] + f[7] + f[9] + f[12] + f[13])) / (1 + Vy);
    double dx = (rho * Vx - f[1] + f[2]) / 4.;
    double dz = (rho * Vz - f[5] + f[6]) / 4.;
    double dy = (rho * Vy) / 12.;

    f[4] = -8 * dy + f[3];
    f[8] = -dx - dy - dz + f[7];
    f[10] = dx - dy - dz + f[9];
    f[11] = dx - dy + dz + f[12];
    f[14] = -dx - dy + dz + f[13];
}

void ZouHeVelocityRight(double* f, double Vx, double Vy, double Vz)
{
    double rho = (f[0] + f[3] + f[4] + f[5] + f[6] + 2 * (f[1] + f[7] + f[10] + f[11] + f[13])) / (1 + Vx);
    double dy = (rho * Vy - f[3] + f[4]) / 4.;
    double dz = (rho * Vz - f[5] + f[6]) / 4.;
    double dx = (rho * Vx) / 12.;

    f[2] = -8 * dx + f[1];
    f[8] = -dx - dy - dz + f[7];
    f[9] = -dx + dy + dz + f[10];
    f[12] = -dx + dy - dz + f[11];
    f[14] = -dx - dz + dz + f[13];
}

void ZouHeVelocityLeft(double* f, double Vx, double Vy, double Vz)
{
    double rho = (f[0] + f[3] + f[4] + f[5] + f[6] + 2 * (f[2] + f[8] + f[9] + f[12] + f[14])) / (1 - Vx);
    double dx = (rho * Vx) / 12.;
    double dy = (rho * Vy - f[3] + f[4]) / 4.;
    double dz = (rho * Vz - f[5] + f[6]) / 4.;

    f[1] = 8 * dx + f[2];
    f[7] = dx + dy + dz + f[8];
    f[10] = dx - dy - dz + f[9];
    f[11] = dx - dy + dz + f[12];
    f[13] = dx + dy - dz + f[14];
}

void ZouHeVelocityTop(double* f, double Vx, double Vy, double Vz)
{
    double rho = (f[0] + f[1] + f[2] + f[3] + f[4] + 2 * (f[5] + f[7] + f[9] + f[11] + f[14])) / (1 + Vz);
    double dx = (f[2] - f[1] + rho * Vx) / 4.0;
    double dy = (f[4] - f[3] + rho * Vy) / 4.0;
    double dz = (rho * Vz) / 12.;

    f[6] = +f[5] - 8 * dz;
    f[8] = +f[7] - dx - dy - dz;
    f[10] = +f[9] + dx - dy - dz;
    f[12] = +f[11] - dx + dy - dz;
    f[13] = +f[14] + dx + dy - dz;
}

void ZouHeVelocityBottom(double* f, double Vx, double Vy, double Vz)
{
    double rho = (f[0] + f[1] + f[2] + f[3] + f[4] + 2 * (f[6] + f[8] + f[10] + f[12] + f[13])) / (Vz - 1);
    double dx = (f[2] - f[1] + rho * Vx) / 4.0;
    double dy = (f[4] - f[3] + rho * Vy) / 4.0;
    double dz = (rho * Vz) / 12.;

    f[5] = f[6] + 8 * dz;
    f[7] = f[8] + dx + dy + dz;
    f[9] = f[10] - dx + dy + dz;
    f[11] = f[12] + dx - dy + dz;
    f[14] = f[13] - dx - dy + dz;
}
}

namespace d3q19
{
void zouHeVelocityLeft(double *f, double vx, double vy, double vz)
{
	double rho = (f[0] + f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18] + 2*(f[2] + f[8] + f[10] + f[12] + f[14]) )/(1 - vx);
	double dx = (rho*vx)/6.;
	double dy = (rho*vy - f[3] + f[4] + f[15] - f[16] + f[17] - f[18])/2.;
	double dz = (rho*vz - f[5] + f[6] + f[15] - f[16] - f[17] + f[18])/2.;
	f[1] = 2*dx + f[2];
	f[7] = dx + dy + f[8];
	f[9] = dx - dy + f[10];
	f[11] = dx + dz + f[12];
	f[13] = dx - dz + f[14];
}

void zouHeVelocityRight(double* f, double vx, double vy, double vz)
{
	double rho =(f[0] +  f[3] + f[4] + f[5] + f[6] + f[15] + f[16] + f[17] + f[18] + 2*( f[7] + f[9] + f[1] + f[11] + f[13] ) )/(1 + vx);
	double dx = (rho*vx)/6.;
	double dy = (rho*vy - f[3] + f[4] + f[15] - f[16] + f[17] - f[18])/2.;
	double dz = (rho*vz - f[5] + f[6] + f[15] - f[16] - f[17] + f[18])/2.;

	f[2] = f[1] - 2*dx;
	f[8] = f[7] - dx - dy;
	f[10] = f[9] - dx + dy;
	f[12] = f[11] - dx - dz;
	f[14] = f[13] - dx + dz;
}

void zouHeVelocityFront(double* f, double vx, double vy, double vz)
{
	double rho = (f[0] + f[1] + f[2]  + f[5] + f[6] + f[11] + f[12] + f[13] + f[14]  + 2*(f[4] + f[9] +  f[8] + f[15] + f[17] ) )/( 1 - vy);
	double dy = (rho*vy)/6.;
	double dx = (rho*vx - f[1] + f[2] - f[11] + f[12] - f[13] + f[14])/2.;
	double dz = (rho*vz - f[5] + f[6] - f[11] + f[12] + f[13] - f[14])/2.;
	f[3] = f[4] + 2*dy;
	f[7] = f[8] + dx + dy;
	f[10] = f[9] - dx + dy;
	f[16] = f[15] + dy + dz;
	f[18] = f[17] + dy - dz;
}

void zouHeVelocityBack(double* f, double vx, double vy, double vz)
{
	double rho = (f[0] + f[1] + f[2] + f[5] + f[6] + f[11] + f[12] + f[13] + f[14]   + 2*(f[3] + f[7] + f[10] + f[16] + f[18] ))/ (1 + vy);
	double dy = (rho*vy)/6.;
	double dx = (rho*vx - f[1] + f[2] - f[11] + f[12] - f[13] + f[14])/2.;
	double dz = (rho*vz - f[5] + f[6] - f[11] + f[12] + f[13] - f[14])/2.;
	f[4] = f[3] - 2*dy;
	f[8] = f[7] - dx - dy;
	f[9] = f[10] + dx - dy;
	f[15] = f[16] - dy - dz;
	f[17] = f[18] - dy + dz;
}

void zouHeVelocityBottom(double* f, double vx, double vy, double vz)
{
	double rho =  (f[0] + f[1] + f[2] + f[3] + f[4] + f[7] + f[8] + f[9] + f[10] + 2*( f[6]  + f[12] + f[13] + f[15] + f[18] ) )/(1 - vz);
	double dy =   (rho*vy - f[3] + f[4] - f[7] + f[8] + f[9] - f[10])/2.;
	double dx =   (rho*vx  -f[1] + f[2] - f[7] + f[8] - f[9] + f[10])/2.;
	double dz =   (rho*vz)/6.;
	f[5] = f[6] + 2*dz;
	f[11] = f[12] + dx + dz;
	f[14] = f[13] - dx + dz;
	f[16] = f[15] + dy + dz;
	f[17] = f[18] - dy + dz;
}

void zouHeVelocityTop(double* f, double vx, double vy, double vz)
{
	double rho =( f[0] + f[1] + f[2] + f[3] + f[4] + f[7] + f[8] + f[9] + f[10] + 2*( f[11] + f[14] + f[16] + f[17] + f[5] ) )/(1 + vz);
	double dy = ( rho*vy - f[3] + f[4] - f[7] + f[8] + f[9] - f[10])/2.;
	double dx = ( rho*vx - f[1] + f[2] - f[7] + f[8] - f[9] + f[10])/2.;
	double dz = ( rho*vz )/6.;

	f[6] = f[5] - 2*dz;
	f[12] = f[11] - dx - dz;
	f[13] = f[14] + dx - dz;
	f[15] = f[16] - dy - dz;
	f[18] = f[17] + dy - dz;
}

void zouHePressureLeft(double *f, double rho)
{
    double Vx =(rho - f[0] - 2*f[2] - f[3] - f[4] - f[5] - f[6] - 2*f[8] - 2*f[10] - 2*f[12] - 2*f[14] - f[15] - f[16] - f[17] - f[18])/rho;
    double dx = (rho*Vx)/12.;
    double dy = (-f[3] + f[4] + f[15] - f[16] + f[17] - f[18])/2.;
    double dz = (-f[5] + f[6] + f[15] - f[16] - f[17] + f[18])/2.;

    f[1] = 4*dx + f[2];
    f[7] = 2*dx + dy + f[8];
    f[9] = 2*dx - dy + f[10];
    f[11] = 2*dx + dz + f[12];
    f[13] = 2*dx - dz + f[14];
}

void zouHePressureRight(double *f, double rho)
{
    double Vx =(-rho + f[0] + 2*f[1] + f[3] + f[4] + f[5] + f[6] + 2*f[7] + 2*f[9] + 2*f[11] + 2*f[13] + f[15] + f[16] + f[17] + f[18])/rho;
    double dx = (rho*Vx)/12.;
    double dy = (-f[3] + f[4] + f[15] - f[16] + f[17] - f[18])/2.;
    double dz = (-f[5] + f[6] + f[15] - f[16] - f[17] + f[18])/2.;

    f[2] = -4*dx + f[1];
    f[8] = -2*dx - dy + f[7];
    f[10] = -2*dx + dy + f[9];
    f[12] = -2*dx - dz + f[11];
    f[14] = -2*dx + dz + f[13];
}

void zouHePressureFront(double *f, double rho)
{
    double Vy =(rho - f[0] - f[1] - f[2] - 2*f[4] - f[5] - f[6] - 2*f[8] - 2*f[9] - f[11] - f[12] - f[13] - f[14] - 2*f[15] - 2*f[17])/rho;
    double dy = (rho*Vy)/12.;
    double dx = (-f[1] + f[2] - f[11] + f[12] - f[13] + f[14])/2.;
    double dz = (-f[5] + f[6] - f[11] + f[12] + f[13] - f[14])/2.;

    f[3] = 4*dy + f[4];
    f[7] = dx + 2*dy + f[8];
    f[10] = -dx + 2*dy + f[9];
    f[16] = 2*dy + dz + f[15];
    f[18] = 2*dy - dz + f[17];
}

void zouHePressureBack(double *f, double rho)
{
    double Vy =(-rho + f[0] + f[1] + f[2] + 2*f[3] + f[5] + f[6] + 2*f[7] + 2*f[10] + f[11] + f[12] + f[13] + f[14] + 2*f[16] + 2*f[18])/rho;
    double dy = (rho*Vy)/12.;
    double dx = (-f[1] + f[2] - f[11] + f[12] - f[13] + f[14])/2.;
    double dz = (-f[5] + f[6] - f[11] + f[12] + f[13] - f[14])/2.;

    f[4] = -4*dy + f[3];
    f[8] = -dx - 2*dy + f[7];
    f[9] = dx - 2*dy + f[10];
    f[15] = -2*dy - dz + f[16];
    f[17] = -2*dy + dz + f[18];
}

void zouHePressureBottom(double *f, double rho)
{
    double Vz =(rho - f[0] - f[1] - f[2] - f[3] - f[4] - 2*f[6] - f[7] - f[8] - f[9] - f[10] - 2*f[12] - 2*f[13] - 2*f[15] - 2*f[18])/rho;
    double dz = (rho*Vz)/12.;
    double dx = (-f[1] + f[2] - f[7] + f[8] - f[9] + f[10])/2.;
    double dy = (-f[3] + f[4] - f[7] + f[8] + f[9] - f[10])/2.;

    f[5] = 4*dz + f[6];
    f[11] = dx + 2*dz + f[12];
    f[14] = -dx + 2*dz + f[13];
    f[16] = dy + 2*dz + f[15];
    f[17] = -dy + 2*dz + f[18];
}

void zouHePressureTop(double *f, double rho)
{
    double Vz =(-rho + f[0] + f[1] + f[2] + f[3] + f[4] + 2*f[5] + f[7] + f[8] + f[9] + f[10] + 2*f[11] + 2*f[14] + 2*f[16] + 2*f[17])/rho;
    double dz = (rho*Vz)/12.;
    double dx = (-f[1] + f[2] - f[7] + f[8] - f[9] + f[10])/2.;
    double dy = (-f[3] + f[4] - f[7] + f[8] + f[9] - f[10])/2.;

    f[6] = -4*dz + f[5];
    f[12] = -dx - 2*dz + f[11];
    f[13] = dx - 2*dz + f[14];
    f[15] = -dy - 2*dz + f[16];
    f[18] = dy - 2*dz + f[17];
}
}
