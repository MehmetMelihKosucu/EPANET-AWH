/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 */
//! \file awhsolver.h
//! \brief Description of the AWHSolver class.


#ifndef AWHSOLVER_H_
#define AWHSOLVER_H_

// Forward declarations
class Network;
class Node;
class Link;
class Junction;
class Pipe;
class PumpCurve;
class MatrixSolver;
class HydSolver;

#include <vector>
#include <string>
#include <unordered_map>
#include <set>
#include <map>
#include <unordered_set>
#include "hydsolver.h"
#include "Core/hydbalance.h"
#include "Core/network.h"
#include "link.h"
#include "node.h"
#include "pipe.h"
#include <unordered_map>
#include <string>
#include "Elements/junction.h"
#include "Elements/valve.h"
#include "Elements/pump.h"
#include "Core/options.h"
#include <set>
#include <map>
#include <unordered_set>
#include "Elements/tank.h"
#include "dualflowhistory.h"

using Vector = std::vector<double>;

//! \class awhSolver
//! \brief A hydraulic solver based on Comprehensive Global Gradient Algorithm.

class AWHSolver : public HydSolver
{
public:

    //! Constructor
    AWHSolver(Network* nw, MatrixSolver* ms);
		
    // Make sure we have a proper destructor
    ~AWHSolver();

    //! Solve hydraulic equations
    //! \param tstep Time step (sec)
    //! \param trials Number of trials for convergence
    //! \param currentTime Current simulation time
    //! \param Xm Auxiliary parameter for specific cases
    //! \return Status of the solver (SUCCESSFUL or FAILED)
    int solve(double& tstep, int& trials, double currentTime);

    // define thresholds:
    static constexpr double PHI_A_D = 0.8;   // dynamic compressibility threshold
    static constexpr double PHI_R_D = 0.2;  // relative dynamic threshold
    static constexpr double PHI_A_I = 0.1;   // inertial threshold
    static constexpr double PHI_R_I = 0.05;  // relative inertial threshold

    // define time steps:
    //static constexpr double Delta_WH = 0.002; // Water Hammer step
    static constexpr double Delta_RWC = 1.0;  // RWC step
    static constexpr double Delta_QS = 1.0; // Quasi-steady step

    // Define the SolverMode enumeration
    enum SolverMode {
        QUASI_STEADY,
        RWC,
        WATER_HAMMER
    };

    enum SystemType { LINEAR_SYSTEM, NONLINEAR_SYSTEM };
    
    SolverMode currentMode = SolverMode::QUASI_STEADY;

    std::vector<double> dQ_start;
    std::vector<double> dQ_end;

    double recommendedTimestep(double currentTime, bool significantHydraulicEvent, bool waterHammerMode) const;

    double minErrorNorm;      //!< Minimum observed error norm
    
private:

    double lamda;

    // Constants
    static const int ARF = 1;  // Adaptive Relaxation Factor
    static const int BRF = 2;  // Backward Relaxation Factor

    // General solver variables
    int    nodeCount;         //!< Number of network nodes
    int    linkCount;         //!< Number of network links
    int    hLossEvalCount;    //!< Number of head loss evaluations
    int    stepSizing;        //!< Newton step sizing method
    int    trialsLimit;       //!< Limit on the number of trials
    bool   reportTrials;      //!< Flag to report summary of each trial
    double tstep;             //!< Time step (sec)
    double errorNorm;         //!< Solution error norm
    double oldErrorNorm;      //!< Previous error norm
    
	double theta;             //!< Relaxation factor for step sizing
	double kappa;             //!< Temporal discretization parameter

    int AWHSolver::sign(double val);

    std::vector<bool> isCompressible;
    double maxWaveSpeed;
    double minSegmentLength;
    void AWHSolver::identifyHydraulicDomains(std::vector<int>& nodeDomain);
	// Network and solver components
    HydBalance hydBalance;  // Instance of the HydBalance class
	
    // Flow characterization and switching
    double phiATolerance;     //!< Threshold for absolute dynamic indicator
    double phiRTolerance;     //!< Threshold for relative dynamic indicator
    double inertiaTolerance;  //!< Threshold for inertia effects
    std::vector<double> phiA; //!< Absolute flow indicators for each link
    std::vector<double> phiR; //!< Relative flow indicators for each link
    std::vector<int>    SC;   //!< Compressibility indicator for each link
    std::vector<int>    SI;   //!< Inertial effects indicator for each link

    // Hydraulic variables
    std::vector<double> dH;   //!< Head change at each node (ft)
    std::vector<double> dQ;   //!< Flow change in each link (cfs)
    std::vector<double> xQ;   //!< Node flow imbalances (cfs)

    // Convergence parameters
    double headErrLimit;      //!< Allowable head error (ft)
    double flowErrLimit;      //!< Allowable flow error (cfs)
    double flowChangeLimit;   //!< Allowable flow change (cfs)
    double flowRatioLimit;    //!< Allowable total flow change / total flow
    
    // Water hammer specific error tracking
    double maxPressureWaveErr;
    double maxCharEqErr;
    double pressureWaveErrLimit;
    double charEqErrLimit;

    // Functions for solver initialization and convergence
    bool hasConverged();
    bool linksChangedStatus();
    
    // Functions to assemble and solve equations
	void setFixedGradeNodes();  //!< Adjust fixed grade status of specific nodes
   
    void setSteadyMatrixCoeffs();       //!< Set coefficients for hybrid equations
    void setSteadyLinkCoeffs(); 	 //!< Set link coefficients
    void setSteadyNodeCoeffs();	   //!< Set node coefficients
    void setSteadyValveCoeffs();	  //!< Set valve coefficients

    void setUnsteadyMatrixCoeffs(double currentTime);       //!< Set coefficients for hybrid equations
    void setUnsteadyLinkCoeffs(double currentTime); 	 //!< Set link coefficients
    void setUnsteadyNodeCoeffs(double currentTime);	   //!< Set node coefficients
    void setUnsteadyValveCoeffs(double currentTime);	  //!< Set valve coefficients
    
    void handleCompressibleFlowCoeffs(Link* link, Node* node1, Node* node2, double currentTime, int n1, int n2);
    void handleQuasiSteadyFlowCoeffs(Link* link, Node* node1, Node* node2, int n1, int n2);

    int findSteadyHeadChanges();       //!< Solve for nodal head changes in Steady State Condition
    void findSteadyFlowChanges();       //!< Solve for link flow changes in Steady State Condition

    int findUnsteadyHeadChanges(double currentTime);
    void findUnsteadyFlowChanges (int currentMode, double currentTime);
    double findUnsteadyErrorNorm(double lamda, double currentTime, double tstep);

    // Functions for solution updates and error handling
    double findStepSize(int trials, double currentTime);
    std::vector<double> Lambda;
	double findErrorNorm(double lamda, double currentTime, double tstep);
    void updateSolution(double lamda, bool useShadowMode = false );
    void reportTrial(int trials, double lamda);
    double AWHSolver::adjustTimeStep(double* currentTime, bool significantHydraulicEvent, bool waterHammerMode);
   
    void setConvergenceLimits();

    int determineSolverMode(SolverMode currentMode, double* currentTime);

    void AWHSolver::classifyNodesAndLinks();

    std::vector<double> currentHeads;  // Current head values at nodes
    std::vector<double> currentFlows;  // Current flow values in links

    int previousMode = 0; // Track previous solver mode

    bool waterHammerMode = false;  // Current solver mode 
    
    // Add these methods
    void switchToWaterHammerMode(double currentTime);
    void switchToStandardMode();
    
    bool AWHSolver::detectSignificantHydraulicEvent(double currentTime);
    
    void AWHSolver::handleQuasiSteadyFlowChange(Link* link, double H_A, double H_B, int linkIndex);
    
    void AWHSolver::updateWHSolution(double lamda);

    // Shadow calculation vectors
    std::vector<double> shadowHeads;
    std::vector<double> shadowFlows;
    double calculateMaxWaveTravel();

};
#endif // AWHSOLVER_H