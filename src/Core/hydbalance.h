/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file hydbalance.h
//! \brief Describes the HydBalance class.
#include <utility>
#include <map>
#include "dualflowhistory.h"


#ifndef HYDBALANCE_H_
#define HYDBALANCE_H_

// Forward declarations

class Network;
class Pipe;
class DualFlowHistory;
class CGGASolver;
//! \class HydBalance
//! \brief Computes the degree to which a network solution is unbalanced.
//!
//! The HydBalance class determines the error in satisfying the head loss
//! equation across each link and the flow continuity equation at each node
//! of the network for an incremental change in nodal heads and link flows.

struct HydBalance
{
    double    maxFlowErr;         //!< max. flow error (cfs)
    double    maxHeadErr;         //!< max. head loss error (ft)
    double    maxCharErr;         //!< max. characteristic equation error 
	double    h1ini;
	double    h1;
	double    h2ini;
	double    h2;
	double    phloss;
	double    gHn;
	double    previousMaxHeadErr; //!< previous max. head loss error (ft)
    double    maxFlowChange;      //!< max. flow change (cfs)
	double    previousMaxFlowChange;  //!< previous max. flow change (cfs)
    double    totalFlowChange;    //!< (summed flow changes) / (summed flows)
    double maxPressureWaveErr;  // Maximum pressure wave propagation error
    double maxCharEqErr;        // Maximum characteristic equation error

    int    maxFlowErrNode;     // index of node with max flow balance error
    int    maxHeadErrLink;     // index of link with max. head loss error
    int    maxFlowChangeLink;  // index of link with max. flow change
    int    maxPressureWaveLink; // index of link with max pressure wave error
    int    maxCharEqLink;      // index of link with max characteristic equation error

    double    evaluate(
                  double lamda, double dH[], double dQ[], double xQ[], Network* nw, int currentTime, double tstep, int currentMode);
    double    findHeadErrorNorm(
		double lamda, double dH[], double dQ[], double xQ[], Network* nw, int currentTime, double tstep);
    double    findFlowErrorNorm(double xQ[], Network* nw);
    
    double HydBalance::evaluateUnsteady(
        double lamda,      // step size
        double dH[],       // change in nodal heads
        double dQ[],       // change in link flows
        double dQ_start[], // change in link start flows
        double dQ_end[],   // change in link end flows
        double xQ[],       // nodal inflow minus outflow
        Network* nw,       // network being analyzed
        double currentTime,   // current time step
        double tstep);      // time step size
    double HydBalance::findUnsteadyHeadErrorNorm(
        double lamda, double dH[], double dQ[], double dQ_start[], double dQ_end[], 
        double xQ[], Network* nw, double currentTime, double tstep);
    void HydBalance::findUnsteadyNodeOutflows(double lamda, double dH[], double xQ[], Network* nw, double tstep);
    double HydBalance::findUnsteadyFlowErrorNorm(double xQ[], Network* nw);
    double HydBalance::findUnsteadyTotalFlowChange(double lamda, double dQ[], double dQ_start[], 
        double dQ_end[], Network* nw);

    CGGASolver* solverPtr = nullptr;
    void HydBalance::identifyHydraulicDomains(std::vector<int>& nodeDomain, Network* nw);
};

#endif
