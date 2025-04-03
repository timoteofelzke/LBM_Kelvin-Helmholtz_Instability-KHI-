/**
 *  @file   drag.h
 *  @author Damylle Cristina Xavier Donati
 *  @brief  Header of functions to compute the drag around an object
 *
 *  @image html drag.png
 */


#ifndef __DRAG_H_INCLUDED__   // if datain.h hasn't been included yet...
#define __DRAG_H_INCLUDED__

#include <string>
#include <vector>

#include "../lbmbase/collision/collision.h"

struct stGeometry;
struct stParameters;

namespace dragProblem
{

struct stProblem
{
    double U;                   /*!< Velocity at the boundary (far from the object). */
    std::string filename;       /*!< Name of the file with output data for the drag problem. */
    double dpX;                 /*!< Change in the x momentum due to interaction with solid object */
    double dpY;                 /*!< Change in the y momentum due to interaction with solid object */
    double dpZ;                 /*!< Change in the z momentum due to interaction with solid object */
    double precision;           /*!< Precision limit in the drag to stop simulation */
    int computeInterval;        /*!< Step interval to revaluate the drag */
    std::vector<int> idList;    /*!< A list with the IDs of the fluid nodes around the object */
    std::vector<int>  iList;    /*!< A list with the directions of the solid nodes in relation the fluid nodes around the object */
};

/**
 * @brief Reads the parameters from the INI file.for the drag around an object problem
 * @param Name of the INI file;
 * @return The stProblem structure with the drag problem parameters read.
 */

stProblem problemParameters( std::string filename );

/**
 * @brief Set up any additional initial configuration to the problem being studied.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 *
 * This function applies some initial conditions or setups to the problem. For the drag problem the implementation is blank and
 * no special initial condition is applied, i. e., the function inplementation is left blank.
 */

void problemInitialize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force);

/**
 * @brief Performe different tasks before the program exits.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 *
 * This function performes some operatiosn before the programs finally exists.  For the drag problem the implementation is blank and
 * no special operations is applied, i. e., the function inplementation is left blank.
 */

void problemFinalize(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

 /**
 * @brief Applies the boundary condition at each step for the drag problem
 * @param iniN      pointer to the distribution function
 * @param geo       pointer to the stGeometry containing all the geometry information
 * @param prm       pointer to the stParameters containing all the initial parameters of the simulation
 * @param drag      pointer to the stProblem containing all the information needed in the drag problem
 * @param step      iteration value
 *
 * Calls all the velocity boundary (DragVel*) functions at once and apply the Neumann boundary condition to the x=(nx-1) plane.
 **/

void problemBoundary(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief Function used to perform any operations between collision and streaming.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param generalParameters The stParameters structure with the simulations parameters
 * @param specificParameters The stProblem structure with problem parameters.
 * @param step The current simulation time step.
 *
 * This function performes some operatiosn before the collision and streaming at each step. The implementation
 * is blank for the drag problem.
 */

void problemPosCollision(double* iniN, stGeometry& geo, stParameters& generalParameters, stProblem& specificParameters,int step);

/**
 * @brief Computes and writes the force over the object in different steps to an CSV file.
 * @param iniN pointer to the memory region where the distribution function of all simulation is stored.
 * @param geo The stGeometry structure with information about the simulation geometry
 * @param prm The stParameters structure with the simulations parameters
 * @param dragData The stProblem structure with drag around object problem parameters .
 * @param step The corrent step of the simulation.
 *
 * Computes the force over the object in different steps and writes the data to an CSV file.
 */

bool problemOutput(double* iniN, stCollision& col , stGeometry &geo, stParameters& gp, stProblem& sp, stForce& force, int step );

/**
 * @brief           Apply the boundary condition to the y=0 plane
 * @param iniN     pointer to the distribution function
 * @param geo       pointer to the stGeometry containing all the geometry information
 * @param u         initial velocity value to the planes
**/

void DragVelYBeg(double* iniN, stGeometry &geo,double u);

/**
 * @brief           apply the boundary condition to the z=0 plane
 * @param iniN     pointer to the distribution function
 * @param geo       pointer to the stGeometry containing all the geometry information
 * @param u         initial velocity value to the planes
**/
void DragVelZBeg(double* iniN, stGeometry &geo,double u);

/**
 * @brief           apply the boundary condition to the x=0 plane
 * @param iniN     pointer to the distribution function
 * @param geo       pointer to the stGeometry containing all the geometry information
 * @param u         initial velocity value to the planes
**/
void DragVelXBeg(double* iniN, stGeometry &geo,double u);

/**
 * @brief           apply the boundary condition to the z=(nz-1) plane
 * @param iniN     pointer to the distribution function
 * @param geo       pointer to the stGeometry containing all the geometry information
 * @param u         initial velocity value to the planes
 **/
void DragVelZEnd(double* iniN, stGeometry &geo,double u);

 /**
 * @brief           apply the boundary condition to the y=(ny-1) plane
 * @param iniN     pointer to the distribution function
 * @param geo       pointer to the stGeometry containing all the geometry information
 * @param u         initial velocity value to the planes
 **/
void DragVelYEnd(double* iniN, stGeometry &geo,double u);

/**
* @brief           Computes some initial operations to speed up the drag computation
* @param geo       pointer to the stGeometry containing all the geometry information
* @param dragData  pointer to the stProblem containing all the information needed in the drag problem
*
* Record the indexes of the sites neighboring the solid. This list of indexes is used
* to compute the drag coefficient.
**/

void setupDragComputation(stGeometry &geo, stProblem& dragData);

/**
* @brief           Computes the force over the object and stores the result in stProblem structure
* @param iniN      pointer to the distribution function
* @param geo       pointer to the stGeometry containing all the geometry information
* @param dragData  pointer to the stProblem containing all the information needed in the drag problem
*
* Computes the force that the object exerts in the fluid and stores the result in stProblem structure. The force is computed by adding up the change in the momentum due the bounce back
* of the particles around the object.
*
* \f{equation}{
* \Delta P_\alpha = \sum_{\vec{x} \in S} \sum_{i \in I_o } 2 f_i(\vec{x}) c_{i,\alpha}
* \f}
*
* where \f$ S\f$ is the set of fluid nodes around the object and  \f$ I_o\f$ is the set of directions where \f$\vec{x} + \vec{c}_i\f$ is a solid node
*
*
*  @image html drag_evaluation.png
*
**/

void computeDrag(double* iniN, stGeometry &geo, stProblem& dragData);


}
#endif
