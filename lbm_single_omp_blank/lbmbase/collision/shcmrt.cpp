#include "shcmrt.h"
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include "mrt.h"
#include "../lattice/lattice.h"
#include "../force/force.h"

using namespace std;

namespace shcMrt
{
    void generateWr(double *ini_R, double *ini_B, double *ini_W,int numberOfPoints)
    {
        #pragma omp parallel for
        for (int n = 0; n < numberOfPoints; n++)
        {
            double* R = ini_R + n*NUM_OF_VEL;
            double* B = ini_B + n*NUM_OF_VEL;
            double* W = ini_W + n*NUM_OF_VEL;

            double rhoR = calculateMass(R);
            double rhoB = calculateMass(B);

            double dr = (rhoR - rhoB);

            for (int i = 0; i < NUM_OF_VEL ; i++ )
            {
                W[i] = dr;
            }
        }
    }

    stCollision collisionParameters(std::string filename)
    {
        using boost::property_tree::ptree;
        stCollision prm;

        ptree pt;
        read_ini(filename.c_str(),pt);

        prm.delta = pt.get<double>("SHCMRT.delta");
        prm.Beta = pt.get<double>("SHCMRT.beta");
        prm.A = pt.get<double>("SHCMRT.a");
//        prm.surfaceTension = pt.get<double>("SHCMRT.surfaceTension");

        prm.tauRed = pt.get<double>("SHCMRT.tau_red");
        prm.kinematicViscosityRed = c_s2 * (prm.tauRed - 0.5);

        prm.tauBlue = pt.get<double>("SHCMRT.tau_blue");
        prm.kinematicViscosityBlue = c_s2 *(prm.tauBlue - 0.5);
        
        prm.alphaBlue = pt.get<double>("SHCMRT.alpha_blue");
        prm.alphaRed = pt.get<double>("SHCMRT.alpha_red");
        prm.drWall = pt.get<double>("SHCMRT.wall_difference");

        prm.s1 = 2 * prm.tauRed * prm.tauBlue /(prm.tauBlue + prm.tauRed);
        prm.s2 = 2 *(prm.tauRed - prm.s1)/prm.delta;
        prm.s3 = -0.5 * prm.s2/prm.delta;

        prm.t1 = prm.s1;
        prm.t2 = 2 * (prm.t1 - prm.tauBlue) / prm.delta;
        prm.t3 = 0.5 * prm.t2 / prm.delta;
        
        return prm;
    }

    void reportCollision(stCollision& prm)
    {
        cout << endl;
        cout << "Collision Model: SHC/MRT";
        cout << endl;
        cout << left << setw(40) << "Relaxation Time Red:";
        cout << right << setw(30) << prm.tauRed << endl;
        cout << left << setw(40) << "Relaxation Time Blue:";
        cout << right << setw(30) << prm.tauBlue << endl;
        cout << left << setw(40) << "Alpha Red:";
        cout << right << setw(30) << prm.alphaRed << endl;
        cout << left << setw(40) << "Alpha Blue:";
        cout << right << setw(30) << prm.alphaBlue << endl;
        cout << left << setw(40) << "A:";
        cout << right << setw(30) << prm.A << endl;
        cout << left << setw(40) << "Beta:";
        cout << right << setw(30) << prm.Beta << endl;
        cout << left << setw(40) << "Density Diffence at the wall";
        cout << right << setw(30) << prm.drWall << endl;
        cout << left << setw(40) << "Viscosity Red:";
        cout << right << setw(30) << prm.kinematicViscosityRed << endl;
        cout << left << setw(40) << "Viscosity Blue:";
        cout << right << setw(30) << prm.kinematicViscosityBlue << endl;
        cout << left << setw(40) << "Surface Tension:";
        cout << right << setw(30) <<"Not Implemented" << endl;
        cout << left << setw(40) << "Interface length";
        cout << right << setw(30) <<"Not Implemented" << endl;

        cout << endl;
    }

    void collisionSwap(double* ini_R, double* ini_B, double* ini_D,  stForce& forceRed, stForce& forceBlue,  stCollision& prm, int numberOfPoints)
    {
        #pragma omp parallel for
        for ( int n = 0; n< numberOfPoints;n++)
        {
            
            double* R = ini_R + n * NUM_OF_VEL;
            double* B = ini_B + n * NUM_OF_VEL;
            double* D = ini_D + n * NUM_OF_VEL;

            //Compute the gradient;
            double grad[3];
            calculateMomentum(D,grad);

            // Color Blind distribution :
            double f[NUM_OF_VEL];
            double m[NUM_OF_VEL];
            double meq[NUM_OF_VEL];
            
            for (int i = 0; i < NUM_OF_VEL; i ++){
                f[i] = R[i] + B[i] ;
            }
                
            double module = sqrt( grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2] );
            double normal[3];

            for (int i = 0; i < 3; i++)
            {
                if ( module < 1E-10) normal[i] = 0;
                else normal[i] = grad[i] / module;
            }

            double rhoR = calculateMass(R);
            double rhoB = calculateMass(B);

            double wR = rhoR / (rhoR + rhoB);
            double wB = 1 - wR;

            double tau = computeTau(prm,rhoR, rhoB);
            double* s;
            s = relaxationMRT(tau);

            calculateMomentsMRT(f,m);
            equilibriumMomentsMRT(meq,m);
            
            for(int i = 4; i < NUM_OF_VEL; i++)
            {
                m[i] = m[i] - s[i] *(m[i] - meq[i]);
            }
            
            reconstructedDistributionMRT(f,m);

            double* Fred     =  forceRed.constant ?  forceRed.constantForce :  forceRed.forceField + n * NUM_OF_DIM;
            double* Fblue    =  forceBlue.constant ? forceBlue.constantForce : forceBlue.forceField + n * NUM_OF_DIM;
            double  Fblind[3] = {  wR*Fred[0] + wB*Fblue[0] ,wR*Fred[1] + wB*Fblue[1]  , wR*Fred[2] + wB*Fblue[2]  };

            double vx, vy, vz, rho;
            calculateMacroEquilibrium( f, vx, vy, vz, rho, Fblind[0] ,  Fblind[1]   , Fblind[2] );

            double omegaF[NUM_OF_VEL];
            forceTerm(omegaF, Fblind[0] ,  Fblind[1] ,  Fblind[2]  ,vx,vy,vz, tau );

            //Add the second collision term :

            double omega[NUM_OF_VEL];
            shcMrtForce(omega,normal,module,prm.A);

            for(int i = 0; i < NUM_OF_VEL; i++ )
            {
               f[i] = f[i] +  omega[i] + omegaF[i];
            }

            shcMrtRecolloring(f,R,B,normal,module,prm.Beta,prm.alphaRed,prm.alphaBlue);

            double temp;
            for (int i = 1; i< NUM_OF_VEL; i = i + 2)
            {
                temp = R[i];
                R[i] = R[i+1];
                R[i+1] = temp;

                temp = B[i];
                B[i] = B[i+1];
                B[i+1] = temp;

            }

            delete[] s;

        }
    }

    double computeTau(stCollision& prm,double& rhoR,double& rhoB)//Está forçando o Tau a ser o do arquivo[DEBUG]
    {
            double Tau;
            double psi = (rhoR - rhoB)/(rhoR + rhoB);
    
            Tau = prm.tauRed;

            if (prm.delta >= psi && psi > 0) // então estaremos em uma zona de transição com mais vermelho que azul;
            {
                Tau = prm.s1 + prm.s2*psi + prm.s3*psi*psi;
            }

            if ( 0 >= psi && psi >= - prm.delta) //então estaremos em uma zona de transição com mais azul que vermelho;
            {

                Tau = prm.t1 + prm.t2*psi + prm.t3*psi*psi;
            }

            if(psi < - prm.delta) //então estaremos no dominio do fluido azul;
            {
                Tau = prm.tauBlue;
            }

        return Tau;
    }
}

namespace d3q19
{

    void shcMrtRecolloring(double *f, double *R, double *B, double* normal,double module, double beta, double alphaRed, double alphaBlue)
    {
        double rhoR = calculateMass(R);
        double rhoB = calculateMass(B);
        double rho = rhoR + rhoB;
        double wr = rhoR/rho;
        double wb = 1 - wr;

        double Cr = 1.5 * (1 - alphaRed);
        double Cb = 1.5 * (1 - alphaBlue);

        R[0] = wr * f[0];
        B[0] = wb * f[0];

        double coef = beta * rho * wr * wb;
        double df = w1 * coef * normal[0];
        R[1] = wr * f[1] - df * Cr;
        B[1] = wb * f[1] + df * Cb;
        R[2] = wr * f[2] + df * Cr;
        B[2] = wb * f[2] - df * Cb;

        df = w1 * coef * normal[1];
        R[3] = wr * f[3] - df * Cr;
        B[3] = wb * f[3] + df * Cb;
        R[4] = wr * f[4] + df * Cr;
        B[4] = wb * f[4] - df * Cb;

        df = w1 * coef *  normal[2];
        R[5] = wr * f[5] - df * Cr;
        B[5] = wb * f[5] + df * Cb;
        R[6] = wr * f[6] + df * Cr;
        B[6] = wb * f[6] - df * Cb;

        df = w2 * coef * (normal[0] + normal[1]);
        R[7] = wr * f[7] - df * Cr;
        B[7] = wb * f[7] + df * Cb;
        R[8] = wr * f[8] + df * Cr;
        B[8] = wb * f[8] - df * Cb;

        df = w2 * coef * (normal[0] - normal[1]);
        R[9] = wr * f[9] - df * Cr;
        B[9] = wb * f[9] + df * Cb;
        R[10] = wr * f[10] + df * Cr;
        B[10] = wb * f[10] - df * Cb;

        df = w2 * coef * (normal[0] + normal[2]);
        R[11] = wr * f[11] - df * Cr;
        B[11] = wb * f[11] + df * Cb;
        R[12] = wr * f[12] + df * Cr;
        B[12] = wb * f[12] - df * Cb;

        df = w2 * coef * (normal[0] - normal[2]);
        R[13] = wr * f[13] - df * Cr;
        B[13] = wb * f[13] + df * Cb;
        R[14] = wr * f[14] + df * Cr;
        B[14] = wb * f[14] - df * Cb;

        df = w2 * coef * (-normal[1] - normal[2]);
        R[15] = wr * f[15] - df * Cr;
        B[15] = wb * f[15] + df * Cb;
        R[16] = wr * f[16] + df * Cr;
        B[16] = wb * f[16] - df * Cb;

        df = w2 * coef * (-normal[1] + normal[2]);
        R[17] = wr * f[17] - df * Cr;
        B[17] = wb * f[17] + df * Cb;
        R[18] = wr * f[18] + df * Cr;
        B[18] = wb * f[18] - df * Cb;
    }

    void shcMrtForce(double* omg, double* normal,double module, double A)
    {
        double C  = 0.5 * A * module;

        omg[0] = C* w0;

        double C1 = C * w1;
        double nc = normal[0];
        omg[1] = C1*(nc*nc - 1.);
        omg[2] = omg[1];

        nc = normal[1];
        omg[3] = C1*(nc*nc - 1.);
        omg[4] = omg[3];

        nc = normal[2];
        omg[5] = C1*(nc*nc - 1.);
        omg[6] = omg[5];

        double C2 = C * w2;
        nc = normal[0] + normal[1];
        omg[7] = C2*(nc*nc - 1.);
        omg[8] = omg[7];

        nc = normal[0] - normal[1];
        omg[9] = C2*(nc*nc - 1.);
        omg[10] = omg[9];

        nc = normal[0] + normal[2];
        omg[11] = C2*(nc*nc - 1.);
        omg[12] = omg[11];

        nc = normal[0] - normal[2];
        omg[13] = C2*(nc*nc - 1.);
        omg[14] = omg[13];

        nc = -normal[1] - normal[2];
        omg[15] = C2*(nc*nc - 1.);
        omg[16] = omg[15];

        nc = -normal[1] + normal[2];
        omg[17] = C2*(nc*nc - 1.);
        omg[18] = omg[17];
    }

}
