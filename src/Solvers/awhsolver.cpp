/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 */
//! \file awhsolver.cpp
//! \brief Implementation of the AWHSolver class.


// Ensure that all necessary headers are included
#include "awhsolver.h"
#include "matrixsolver.h"
#include "Core/network.h"
#include "Core/constants.h"
#include "Elements/control.h"
#include "Core/hydbalance.h"
#include "Core/options.h" 
#include "Elements/junction.h"
#include "Elements/tank.h"
#include "Elements/valve.h"
#include "Elements/node.h"
#include "Elements/link.h"
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Elements/pumpcurve.h"
#include "Elements/curve.h"
#include "Elements/reservoir.h"
#include "Elements/emitter.h"
#include "Models/headlossmodel.h"
#include "Models/leakagemodel.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <numeric>  // Required for std::accumulate
#include "sparspak.h"
#include "Models/demandmodel.h"
#include "sparspaksolver.h"
#include "project.h"
#include <unordered_map>
#include <string>
#include <set>
#include <unordered_set>
#include <queue>
#include "dualflowhistory.h"
#include "constants.h"

#include <iomanip>
using namespace std;
static const string s_Trial = "    Trial ";
static const string s_IllConditioned = "  Hydraulic matrix ill-conditioned at node ";
static const string s_StepSize = "    Step Size   = ";
static const string s_TotalError = "    Error Norm  = ";
static const string s_HlossEvals = "    Head Loss Evaluations = ";
static const string s_HeadError = "    Head Error  = ";
static const string s_ForLink = " for Link ";
static const string s_FlowError = "    Flow Error  = ";
static const string s_AtNode = " at Node ";
static const string s_FlowChange = "    Flow Change = ";
static const string s_TotFlowChange = "    Total Flow Change Ratio = ";
static const string s_NodeLabel = "  Node ";
static const string s_FGChange = "    Fixed Grade Status changed to ";

//-----------------------------------------------------------------------------

// error norm threshold for status checking
static const double ErrorThreshold = 1e-6;
static const double Huge = numeric_limits<double>::max();

// step sizing enumeration
enum StepSizing { FULL, RELAXATION, LINESEARCH, BRF, ARF };

// Updated constructor with network state management initialization
AWHSolver::AWHSolver(Network* nw, MatrixSolver* ms) : HydSolver(nw, ms), 
    hydBalance(),  
    waterHammerMode(false),
    previousMode(0) 
{  

    if (!nw || !ms) {
        throw std::invalid_argument("Network or MatrixSolver cannot be null.");
    }

    // Store network and matrix solver pointers
    network = nw;
    matrixSolver = ms;

    nodeCount = nw->nodes.size();
    linkCount = nw->links.size();

    // Resize vectors to match network dimensions
    dH.resize(nodeCount, 0.0);    dQ.resize(linkCount, 0.0);    xQ.resize(nodeCount, 0.0);    
    dQ_start.resize(linkCount, 0.0);    dQ_end.resize(linkCount, 0.0);
    SC.resize(linkCount, 0);    SI.resize(linkCount, 0);    phiA.resize(linkCount, 0.0);
    phiR.resize(linkCount, 0.0);
 SI.resize(linkCount, 0);    // Inertial effects indicators
    phiA.resize(linkCount, 0.0);
    phiR.resize(linkCount, 0.0);

    // Initialize convergence parameters
    errorNorm = std::numeric_limits<double>::max();
    headErrLimit = 0.1;  // Default head error limit
    flowErrLimit = 0.5;  // Default flow error limit
    flowChangeLimit = 0.1; // Default flow change limit

    // Initialize tolerances for dynamic flow indicators
    phiATolerance = 0.05;   // Threshold for absolute indicator
    phiRTolerance = 0.1;    // Threshold for relative indicator
    inertiaTolerance = 0.2; // Threshold for inertial effects

    if (network->option(Options::STEP_SIZING) == "RELAXATION" )
        stepSizing = RELAXATION;
    else if (network->option(Options::STEP_SIZING) == "LINESEARCH" )
        stepSizing = LINESEARCH;
    else if (network->option(Options::STEP_SIZING) == "BRF")
        stepSizing = BRF;
    else if (network->option(Options::STEP_SIZING) == "ARF")
        stepSizing = ARF;
    else stepSizing = FULL;

    // Log initialization (optional)
    std::cout << "AWHSolver initialized with "
        << nodeCount << " active nodes and "
        << linkCount << " active links." << std::endl;
}

// Updated destructor to clean up all resources
AWHSolver::~AWHSolver() {
    // Clean up vectors
    dH.clear();
    dQ.clear();
    xQ.clear();
    SC.clear();
    SI.clear();
    phiA.clear();
    phiR.clear();
    
}

int AWHSolver::solve(double& tstep_, int& trials, double currentTime) {
    // Initialize variables
    tstep = tstep_;  
    trials = 1;

    // Initialize variables
    double lamda = 1.0;
    bool statusChanged = true;
    bool converged = false;
    double errorNorm = Huge;
    minErrorNorm = 1e9;
    double dl = 1.0;
    int currentMode = static_cast<SolverMode>(0);

    theta = network->option(Options::TIME_WEIGHT);
    theta = min(theta, 1.0);
    if (theta > 0.0) theta = max(theta, 0.5);

    kappa = network->option(Options::TEMP_DISC_PARA);
    kappa = min(kappa, 1.0);
    if (kappa < 0) kappa = 0; //

    // Initialize solver
    minErrorNorm = HUGE;

    double simulationEndTime = network->option(Options::TOTAL_DURATION);
    
    // CRITICAL: Check for hydraulic events FIRST
    bool significantHydraulicEvent = detectSignificantHydraulicEvent(currentTime);

    if (currentTime == 0) {
        significantHydraulicEvent = false; // No event at time zero
    }
    
    // If a significant hydraulic event is detected, skip steady-state and go straight to water hammer
    if (significantHydraulicEvent || waterHammerMode) {

        switchToWaterHammerMode(currentTime);
        //prepareAllPipesForawh();

        // Use the new function to get the recommended time step
        double stepBeforeAdjustment = tstep;
        tstep = adjustTimeStep(&currentTime, significantHydraulicEvent, waterHammerMode);
        tstep_ = tstep;  // Update the reference parameter

        simulationEndTime = network->option(Options::TOTAL_DURATION);
        if (currentTime + tstep > simulationEndTime) {
            // Adjust time step to exactly hit the end time
            tstep = simulationEndTime - currentTime;
            tstep_ = tstep;
                    
            if (tstep <= 0) {
                network->msgLog << "\nSimulation reached end time of " << simulationEndTime << " seconds.";
                    return HydSolver::SUCCESSFUL;
            }
                    
                std::cout << "Final time step adjusted to " << tstep 
                            << " to exactly reach simulation end time of " 
                            << simulationEndTime << " seconds" << std::endl;
        }

        // ==== NEWTON-RAPHSON ITERATION PARAMETERS ====
        const int maxIterations = 33;          // Maximum number of iterations
        const double convergenceTol = 1e-6;    // Convergence tolerance
        int iterations = 0;                      // Iteration counter
        double prevErrorNorm = HUGE_VAL;       // Previous error norm for comparison
        double errorNorm = HUGE_VAL;    // Current error norm
        //minErrorNorm = 10000000; 
        converged = false;                // Convergence flag
        double damping = 1.0;                  // Solution damping factor (can be adjusted for stability)
        statusChanged = true;            // Flag to check for link status changes
            
        // Classify nodes and links for awh implementation
        classifyNodesAndLinks();
            
        // Initialize history storage
        network->msgLog << "\n\nInitializing history storage...";
        double maxWaveTime = 0.0;
        for (Link* link : network->links) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe) {
                double waveTime = pipe->length / pipe->getWaveSpeed();
                maxWaveTime = std::max(maxWaveTime, waveTime);
            }
        }
            
        double historyDuration = 2.5 * maxWaveTime;
        network->msgLog << "\n  Max wave time: " << maxWaveTime << " s";
        network->msgLog << "\n  History duration: " << historyDuration << " s";
            
        // Initialize FlowHistoryManager
        double startTime = currentTime;  // Use current time as start time
            
        // Initial error calculation
        //errorNorm = findUnsteadyErrorNorm(1.0, currentTime, tstep);

        setConvergenceLimits();

        dl = 1.00000;
        double lambda = 1.0;
            
        // Enhanced Newton-Raphson loop with step sizing for water hammer
        while (!converged && iterations < maxIterations) {
            network->msgLog << "\n--- Newton-Raphson Iteration " << iterations+1 << " ---";
                
            // Store the error from the previous iteration - this is our baseline to beat
            prevErrorNorm = errorNorm;
                
            // Perform initial setup for this iteration
            
            setFixedGradeNodes();
                
            // Handle any status changes from previous iteration
            if (statusChanged) {
                // Recalculate error after status change
                prevErrorNorm = findUnsteadyErrorNorm(1.0, currentTime, tstep);
                network->msgLog << "\nError after status change: " << prevErrorNorm;
                lambda = 1.0;  // Reset to full step after status change
            }
            statusChanged = false;
                
            // Initialize step size search parameter
            dl = 1.0;  // This controls the granularity of lambda search
                
            // STEP 1: Solve for head and flow changes
            int unsteadyErrorCode = findUnsteadyHeadChanges(currentTime);
            if (unsteadyErrorCode >= 0) {
                network->msgLog << "\nUnsteady system failed at position: " << unsteadyErrorCode;
                return HydSolver::FAILED_ILL_CONDITIONED;
            }
            findUnsteadyFlowChanges(currentMode, currentTime);
                
            // Standard approach without step sizing
            double lambda = 1.0;
            errorNorm = findUnsteadyErrorNorm(lambda, currentTime, tstep);
            updateWHSolution(lambda);
                    
            converged = hasConverged();
            if (converged) {
                statusChanged = linksChangedStatus();
            }
                    
            if (converged && !statusChanged) break;
                
            iterations++;
                
            std::cout << "Iteration " << iterations << " completed with error=" << errorNorm << std::endl;
        }
            

        // ==== REPORT FINAL CONVERGENCE STATUS ====
        if (converged) {
            network->msgLog << "\nWater hammer solution converged after " << iterations << " iterations";
            network->msgLog << "\nFinal error norm: " << errorNorm;
        } 
        else {
            network->msgLog << "\nWarning: Water hammer solution did not converge after " << maxIterations << " iterations";
            network->msgLog << "\nFinal error norm: " << errorNorm;
                // Might still be usable if error is low enough
            if (errorNorm < 0.1) {
                network->msgLog << "\nContinuing with non-converged solution (error is small)";
            } else {
                    // Consider returning a warning or error code
                return HydSolver::FAILED_NO_CONVERGENCE;
            }
        }

        // ==== PHASE 6: TIME STEP ADJUSTMENT ====
        // Set appropriate time step for next iteration
        double simulationEndTime = network->option(Options::TOTAL_DURATION);
        double remainingTime = simulationEndTime - currentTime;
            
        if (Delta_WH > remainingTime && remainingTime > 0) {
            tstep = remainingTime;
        } else {
            tstep = Delta_WH;  // Use water hammer time step
        }
        tstep_ = tstep; // Update reference parameter
            
        // Log progress
        network->msgLog << "\nTime: " << currentTime << ", Next step: " << tstep
                        << ", Mode: Water Hammer";
                
        // Return result based on convergence
        return converged ? HydSolver::SUCCESSFUL : HydSolver::FAILED_NO_CONVERGENCE;    
       
    }
    else if (!significantHydraulicEvent) {
        // Switch back to standard mode
        switchToStandardMode();

        // THIS IS CRITICAL: Check if we need to switch modes
        // If we're in water hammer mode but need to switch back
        if (waterHammerMode) {
            // Check if we still need water hammer mode
            currentMode = determineSolverMode(static_cast<SolverMode>(previousMode), &currentTime);
            
            if (currentMode != WATER_HAMMER) {
                // Switch back to standard mode
                switchToStandardMode();
            }
        }

        // ... set values for convergence limits

        setConvergenceLimits();

        while (trials <= trialsLimit)
        {
            // ... save current error norm

            oldErrorNorm = errorNorm;

            // ... determine which nodes have fixed heads (e.g., PRVs)

            setFixedGradeNodes();

            // ... re-compute error norm if links change status

            if (statusChanged)
            {
                oldErrorNorm = findErrorNorm(0.0, currentTime, tstep);
                lamda = 1.0;
            }
            statusChanged = false;
            dl = 1.000000000000000;

            int errorCode = findSteadyHeadChanges();
            if (errorCode >= 0)
            {
                Node* node = network->node(errorCode);
                network->msgLog << endl << s_IllConditioned << node->name;
                return HydSolver::FAILED_ILL_CONDITIONED;
            }
            findSteadyFlowChanges();

            errorNorm = findStepSize(trials, currentTime);
            updateSolution(1.00, false);

            // ... check for convergence

            if (reportTrials) reportTrial(trials, lamda);
            converged = hasConverged();

            // ... if close to convergence then check for any link status changes

            if (converged && errorNorm < ErrorThreshold )
            {
                statusChanged = linksChangedStatus();
            }

            // ... check if the current solution can be accepted

            if (converged && !statusChanged) break;
            trials++;
        }
    
        return HydSolver::SUCCESSFUL;
    }

}

bool AWHSolver::detectSignificantHydraulicEvent(double currentTime) {
    // 1. Check for valve operations in this time step
    for (int i = 0; i < network->linkCount(); i++) {
        Link* link = network->link(i);
        if (link->type() == Link::VALVE ) {
            Valve* valve = dynamic_cast<Valve*>(link);
            if (!valve) continue;

            // Compare current setting/status with previous setting/status
            double previousSetting = valve->pastSetting;
            double currentSetting = valve->setting;
            
            // If valve setting changed significantly (closing or opening)
            if (fabs(currentSetting - previousSetting) > 0.1) {
                return true;
            }
            
            // Or if valve status changed (open/closed)
            if (link->status != link->previousStatus) {
                return true;
            } // */
        }

        if (link->type() == Link::PUMP) {
            
            Pump* pump = dynamic_cast<Pump*>(link);
            if (!pump) continue;

            // Compare current setting/status with previous setting/status
            double previousSetting = pump->pastSetting;
            double currentSetting = pump->setting;

            // If pump setting changed significantly (closing or opening)
            if (fabs(currentSetting - previousSetting) > 0.1) {
                return true;
            }

            // Or if pump status changed (open/closed)
            if (link->status != link->previousStatus) {
                return true;
            } // */
        }

    }
    
    return false;
}

//  Compute the error norm associated with a given step size for unsteady analysis.
double AWHSolver::findUnsteadyErrorNorm(double lamda, double currentTime, double tstep)
{
    hLossEvalCount++;
    return hydBalance.evaluateUnsteady(lamda, &dH[0], &dQ[0], &dQ_start[0], &dQ_end[0],
                                      &xQ[0], network, currentTime, tstep);
}

double AWHSolver::recommendedTimestep(double currentTime, bool significantHydraulicEvent, bool waterHammerMode) const 
{
    // Time zero or startup condition
    if (currentTime <= 1e-12) {
        return Delta_QS;  // 1.0
    }

    if (significantHydraulicEvent || waterHammerMode) {
        return Delta_WH;  
    }
    else {
        return Delta_QS;  
    }
}

void AWHSolver::setConvergenceLimits()
{
    // ... maximum trials
    trialsLimit = network->option(Options::MAX_TRIALS);

    // ... limit on relative flow change ratio (old EPANET2 criterion)
    flowRatioLimit = network->option(Options::RELATIVE_ACCURACY);

    // ... tolerance in satisfying hydraulic balances
    headErrLimit = network->option(Options::HEAD_TOLERANCE) /
        network->ucf(Units::LENGTH);
    flowErrLimit = network->option(Options::FLOW_TOLERANCE) / network->ucf(Units::FLOW);

    // ... limit on largest flow change
    flowChangeLimit = network->option(Options::FLOW_CHANGE_LIMIT) /
        network->ucf(Units::FLOW);

    // ... use a default head error limit if need be
    if (flowRatioLimit == 0.0 && headErrLimit == 0.0 &&
        flowErrLimit == 0.0 && flowChangeLimit == 0.0)
    {
        headErrLimit = 0.005;
    }

    // ... convert missing limits to a huge number
    if (flowRatioLimit == 0.0) flowRatioLimit = Huge;
    if (headErrLimit == 0.0) headErrLimit = Huge;
    if (flowErrLimit == 0.0) flowErrLimit = Huge;
    if (flowChangeLimit == 0.0) flowChangeLimit = Huge;
}

int AWHSolver::findUnsteadyHeadChanges(double currentTime)
{
    // ... setup the coeff. matrix of the GGA linearized system

    setUnsteadyMatrixCoeffs(currentTime);

    // ... temporarily use the head change array dH[] to store new heads

    double* h = &dH[0];

    int errorCode = matrixSolver->solve(nodeCount, h);
    if (errorCode >= 0) return errorCode;

    // ... save new heads as head changes

    for (int i = 0; i < nodeCount; i++)
    {
        dH[i] = h[i] - network->nodes[i]->head;
    }

    // ... return a negative number indicating that
    //     the matrix solver ran successfully

    return -1;
}

void AWHSolver::findUnsteadyFlowChanges(int currentMode, double currentTime) {
    // Main loop for all links in the network
    for (int i = 0; i < linkCount; i++) {
        Link* link = network->links[i];
        if (!link) continue;
        
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;
        
        // Reset flow changes
        dQ_start[i] = 0.0;
        dQ_end[i] = 0.0;
        dQ[i] = 0.0;
        
        // Get nodal heads (including changes from linear system solution)
        double H_A = link->fromNode->head + dH[n1];
        double H_B = link->toNode->head + dH[n2];

        // COMPRESSIBLE FLOW MODEL (WATER HAMMER)
        if (link->flowModel == Link::COMPRESSIBLE) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue; // Skip if not a pipe
            
            // === CALCULATE CHARACTERISTIC PARAMETERS ===
            
            // Physical constants
            double c = pipe->getWaveSpeed();
            double A = pipe->getArea();
            double L = pipe->length;
            double g = GRAVITY;
            double Bj = c / (g * A);  // Wave impedance term (B')
            double D = pipe->diameter;
            double epsilon = 0;  // Friction integration parameter
            
            // Calculate wave travel time and reach-back time
            double physicalWaveTime = L / c;
            double reachBackTime;
            
            // Determine reach-back index based on number of reaches
            int reachBackSteps = pipe->numReaches;
            double waveTravel = reachBackSteps * tstep;  // Wave travel distance

            // Retrieve reach-back flows and heads from history
            double Qa = pipe->pastStartFlow;  // Default fallback
            double Qb = pipe->pastEndFlow;    // Default fallback
            double Ha = pipe->fromNode->pastHead;  // Default fallback
            double Hb = pipe->toNode->pastHead;    // Default fallback
            
            // Use the FlowHistoryManager to get reach-back values
            FlowHistoryResult historyResult = FlowHistoryManager::getInstance().getReachBackValues(pipe, currentTime, physicalWaveTime, network);
            
            // If history was found, use the values from the result
            if (historyResult.found) {
                Qa = historyResult.startFlow;
                Qb = historyResult.endFlow;
                Ha = historyResult.startHead;
                Hb = historyResult.endHead;
            }
    
            // === FRICTION AND RESISTANCE CALCULATIONS (CORRECTED) ===
            
            double KA, KB;

            if (network->option(Options::HEADLOSS_MODEL) == "D-W")
            {
                double frictionFactor = pipe->getFrictionFactor(pipe->getRe(pipe->flow, network->option(Options::KIN_VISCOSITY)));
                double K = frictionFactor * pipe->length / (2 * g * D * A * A);

                double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(Qa, network->option(Options::KIN_VISCOSITY)));
                double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(Qb, network->option(Options::KIN_VISCOSITY)));
                // Split K for each end as needed (KA for upstream, KB for downstream)
                KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
                KB = frictionFactorB * pipe->length / (2 * g * D * A * A);

            }
            else
            {
                KA = pipe->resistance;
                KB = pipe->resistance;
            }
            
            const double MIN_RESISTANCE = 1e-10;
            KA = std::max(KA, MIN_RESISTANCE);
            KB = std::max(KB, MIN_RESISTANCE);
            
            // Unsteady friction terms
            double kUFa = pipe->computeShearDecayCoeff(Qa);
            double kUFb = pipe->computeShearDecayCoeff(Qb);
            double Ua = 0.0;
            double Ub = 0.0;
            
            if (kUFa > 0.0 && kUFb > 0.0) {
                auto sign = [](double val) -> double { 
                    return (val > 0.0) ? 1.0 : ((val < 0.0) ? -1.0 : 0.0); 
                };
                
                // CORRECTED: Unsteady friction terms
                Ua = 0.5 * Bj * kUFb * (Qb - 2.0 * sign(Qa) * std::abs(Qb - Qa));
                Ub = 0.5 * Bj * kUFa * (Qa - 2.0 * sign(Qb) * std::abs(Qa - Qb));
            }
            
            // CORRECTED: Characteristic equation coefficients
            double BA, BB, RA, RB;
            
            if (network->option(Options::HEADLOSS_MODEL) == "D-W") {
                // For boundary A: using C+ from point b
                BA = Bj + epsilon * KB * std::abs(Qb) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * std::abs(Qb)) * Qb - Hb + Ub;
                
                // For boundary B: using C- from point a
                BB = Bj + epsilon * KA * std::abs(Qa) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * std::abs(Qa)) * Qa + Ha + Ua;
            }
            else if (network->option(Options::HEADLOSS_MODEL) == "H-W") {
                // For boundary A: using C+ from point b
                BA = Bj + epsilon * KB * pow(std::abs(Qb), 0.852) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * pow(std::abs(Qb), 0.852)) * Qb - Hb + Ub;
                
                // For boundary B: using C- from point a
                BB = Bj + epsilon * KA * pow(std::abs(Qa), 0.852) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * pow(std::abs(Qa), 0.852)) * Qa + Ha + Ua;
            }
            
            pipe->charB_plus = BA;
            pipe->charB_minus = BB;
            pipe->charC_plus = RA;
            pipe->charC_minus = RB;
            
            // === CALCULATE NEW FLOWS ===
            
            double newStartFlow = (H_A + RA) / BA;
            double newEndFlow = (-H_B + RB) / BB;

            double avgFlow = (newStartFlow + newEndFlow) / 2;

            // Limit extreme values
            if(std::abs(newStartFlow) > 1e9) {
                network->msgLog << "\nWarning: Extreme flow at pipe " << pipe->name;
                newStartFlow = 0.0;
            }
            if(std::abs(newEndFlow) > 1e9) {
                newEndFlow = 0.0;
            }
            
            // Calculate and store flow changes
            dQ_start[i] = newStartFlow - pipe->startFlow;
            dQ_end[i] = newEndFlow - pipe->endFlow;
            dQ[i] = 0.5 * (dQ_start[i] + dQ_end[i]);
        }
        // QUASI-STEADY FLOW MODEL
        else {
            handleQuasiSteadyFlowChange(link, H_A, H_B, i);
        }
    }
}

// Helper function for quasi-steady flow calculations
void AWHSolver::handleQuasiSteadyFlowChange(Link* link, double H_A, double H_B, int linkIndex) {
    int n1 = link->fromNode->index;
    int n2 = link->toNode->index;
    
    // Handle pressure regulating valves
    if (link->hGrad == 0.0) {
        if (link->isPRV()) {
            dQ[linkIndex] = -xQ[n2] - link->flow;
        }
        else if (link->isPSV()) {
            dQ[linkIndex] = xQ[n1] - link->flow;
        }
        return;
    }

    // Regular flow change based on head difference
    double dh = (link->fromNode->head + dH[n1]) - (link->toNode->head + dH[n2]);

	double dq = -(dh - link->hLoss) / (link->hGrad);

    // Special case for constant HP pumps
    if (link->isHpPump() && link->status == Link::LINK_OPEN && dq > link->flow) {
        dq = link->flow / 2.0;
    }
        
    dQ[linkIndex] = -dq;

    // Save flow changes (same for start and end in quasi-steady flow)
    dQ_start[linkIndex] = -dq;
    dQ_end[linkIndex] = -dq; // */
}

//  Compute the coefficient matrix of the linearized set of equations for heads.
void AWHSolver::setUnsteadyMatrixCoeffs(double currentTime)
{
    // Ensure xQ is properly sized
    if (xQ.size() < nodeCount) {
        xQ.resize(nodeCount, 0.0);
    } else {
        // Zero out the vector using vector operations instead of memset
        std::fill(xQ.begin(), xQ.begin() + nodeCount, 0.0);
    }
    
    matrixSolver->reset();
    setUnsteadyLinkCoeffs(currentTime);
    setUnsteadyNodeCoeffs(currentTime);
    setUnsteadyValveCoeffs(currentTime);
}

//  Compute matrix coefficients for links in unsteady flow conditions.
void AWHSolver::setUnsteadyLinkCoeffs(double currentTime)
{
    // Clear node flow balances
    for (int i = 0; i < nodeCount; i++) {
        if (i < xQ.size()) {
            xQ[i] = 0.0;
        }
    }

    // Create and initialize the vector tracking nodes connected to incompressible links
    std::vector<bool> nodeHasIncompressibleLink(nodeCount, false);

    std::vector<int> nodeDomain(nodeCount, -1);
    identifyHydraulicDomains(nodeDomain);
    
    // Identify nodes connected to incompressible flow links
    for (int i = 0; i < nodeCount; i++) {
        Node* node = network->nodes[i];
        if (!node) continue;
        
        // Check all links connected to this node
        for (Link* link : network->getLinksConnectedToNode(node)) {
            // Verify the link is connected to this node (redundant but safe check)
            if (link->fromNode == node || link->toNode == node) {
                // Check if this is an incompressible flow link
                if (link->flowModel != Link::COMPRESSIBLE) // && link->status != Link::LINK_CLOSED) {
                {    
                    nodeHasIncompressibleLink[i] = true;
                    break; 
                }
            }
        }
    }
    
    // Process all links for flow contributions
    for (int j = 0; j < linkCount; j++)
    {
        Link* link = network->links[j];
        if (!link) continue;

        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        if (!node1 || !node2) continue;
        
        int n1 = node1->index;
        int n2 = node2->index;

        // Skip links with zero head gradient
        if (link->hGrad == 0.0) {
            continue;
        }

        // Update node flow balances based on flow model and node connectivity
        if (link->flowModel != Link::COMPRESSIBLE) {
            xQ[n1] -= link->flow;
            xQ[n2] += link->flow;
        }

        // Handle hydraulic coefficients based on flow model
        if (link->type() == Link::VALVE || link->type() == Link::PUMP) {
            handleQuasiSteadyFlowCoeffs(link, node1, node2, n1, n2);
        }
        else if (link->flowModel == Link::COMPRESSIBLE) {            
            handleCompressibleFlowCoeffs(link, node1, node2, currentTime, n1, n2);
        }
        else {
            bool hasInertia = (j < SI.size() && SI[j] > 0);
            
            if (hasInertia) {
                handleQuasiSteadyFlowCoeffs(link, node1, node2, n1, n2);
            }
            else {
                handleQuasiSteadyFlowCoeffs(link, node1, node2, n1, n2);
            }
        }
    }
}

void AWHSolver::identifyHydraulicDomains(std::vector<int>& nodeDomain) {
    int domainCount = 0;
    
    // Reset all domain assignments
    for (int i = 0; i < nodeDomain.size(); i++) {
        nodeDomain[i] = -1;
    }
    
    // Use a breadth-first search to identify connected domains
    for (int i = 0; i < nodeCount; i++) {
        if (nodeDomain[i] >= 0) continue;  // Already assigned
        
        // Start a new domain
        nodeDomain[i] = domainCount++;
        
        // Use BFS to find all connected nodes
        std::queue<int> nodeQueue;
        nodeQueue.push(i);
        
        while (!nodeQueue.empty()) {
            int currentNodeIdx = nodeQueue.front();
            nodeQueue.pop();
            
            Node* currentNode = network->nodes[currentNodeIdx];
            if (!currentNode) continue;
            
            // Check all connected links
            for (Link* link : network->getLinksConnectedToNode(currentNode)) {
                // Skip closed links - they create domain boundaries
                if (link->status == Link::LINK_CLOSED) continue;
                
                // Find the node on the other side
                Node* otherNode = (link->fromNode == currentNode) ? 
                                   link->toNode : link->fromNode;
                
                if (!otherNode) continue;
                
                int otherNodeIdx = otherNode->index;
                
                // If not yet assigned, add to this domain
                if (nodeDomain[otherNodeIdx] < 0) {
                    nodeDomain[otherNodeIdx] = nodeDomain[currentNodeIdx];
                    nodeQueue.push(otherNodeIdx);
                }
            }
        }
    }
}

//  Compute matrix coefficients for dynamic tanks and external node outflows.

void AWHSolver::setUnsteadyNodeCoeffs(double currentTime)
{
    for (int i = 0; i < nodeCount; i++)
    {
        // ... if node's head not fixed

        Node* node = network->node(i);
        if ( !node->fixedGrade )
        {
            // ... for dynamic tanks, add area terms to row i
            //     of the head solution matrix & r.h.s. vector

            if ( node->type() == Node::TANK && theta != 0.0 )
            {
                Tank* tank = static_cast<Tank*>(node);


                /*double a = tank->area / (theta * tstep);
				matrixSolver->addToDiag(i, a);

				a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
				matrixSolver->addToRhs(i, a); // */

				if (tank->head == tank->pastHead)
				{

					double a = tank->area / (theta * tstep);
					matrixSolver->addToDiag(i, a);

					a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
					matrixSolver->addToRhs(i, a); // 
				}
				
				else
				{
					double a = node->qGrad + (tank->area / (theta * tstep));
					matrixSolver->addToDiag(i, a);

					double b = (tank->pastArea) * tank->pastHead / (theta * tstep) + (1.0 - theta) * tank->pastOutflow / theta;
					matrixSolver->addToRhs(i, b); //

                    xQ[i] -= tank->outflow;
				} // */
            }

            // ... for junctions, add effect of external outflows

            else if ( node->type() == Node::JUNCTION )
            {
                // ... update junction's net inflow
                xQ[i] -= node->outflow;
                double pastOutflow = -(1.0 - theta) * node->pastOutflow / theta;
                matrixSolver->addToDiag(i, node->qGrad);
                matrixSolver->addToRhs(i, node->qGrad * node->head);
                matrixSolver->addToRhs(i, pastOutflow);
            }

            // ... add node's net inflow to r.h.s. row
            matrixSolver->addToRhs(i, (double)xQ[i]);
        }

        // ... if node has fixed head, force solution to produce it

        else
        {
            matrixSolver->setDiag(i, 1.0);
            matrixSolver->setRhs(i, node->head);
        }
    }
}

//-----------------------------------------------------------------------------

//  Compute matrix coefficients for pressure regulating valves.

void  AWHSolver::setUnsteadyValveCoeffs(double currentTime)
{
    for (Link* link : network->links)
    {
        // ... skip links that are not active pressure regulating valves

        if (link->hGrad > 0.0) continue;

        // ... determine end node indexes of link

        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... add net inflow of downstream node of a PRV to the
        //     r.h.s. row of its upstream node

        if (link->isPRV())
        {
            matrixSolver->addToRhs(n1, (double)xQ[n2]);
        }

        // ... add net inflow of upstream node of a PSV to the
        //     r.h.s. row of its downstream node

        if (link->isPSV())
        {
            matrixSolver->addToRhs(n2, (double)xQ[n1]);
        }
    }
}

//  Find changes in nodal heads by solving a linearized system of equations.
int AWHSolver::findSteadyHeadChanges()
{
    // ... setup the coeff. matrix of the GGA linearized system

    setSteadyMatrixCoeffs();

    double* h = &dH[0];

    int errorCode = matrixSolver->solve(nodeCount, h);
    if (errorCode >= 0) return errorCode;

    // ... save new heads as head changes

    for (int i = 0; i < nodeCount; i++)
    {
        dH[i] = h[i] - network->nodes[i]->head;
    }

    return -1;
}

//  Find the changes in link flows resulting from a set of nodal head changes.

void AWHSolver::findSteadyFlowChanges()
{
    for (int i = 0; i < linkCount; i++)
    {
        // ... get link object and its end node indexes

        dQ[i] = 0.0;
        Link* link = network->link(i);
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... flow change for pressure regulating valves

        if (link->hGrad == 0.0)
        {
            if (link->isPRV()) dQ[i] = -xQ[n2] - link->flow;
            if (link->isPSV()) dQ[i] = xQ[n1] - link->flow;
            continue;
        }

        // ... apply GGA flow change formula:

        double dh = (link->fromNode->head + dH[n1]) -
            (link->toNode->head + dH[n2]);
        double dq = (link->hLoss - dh) / link->hGrad;

        // ... special case to prevent negative flow in constant HP pumps

        if (link->isHpPump() &&
            link->status == Link::LINK_OPEN &&
            dq > link->flow) dq = link->flow / 2.0;

        // ... save flow change

        dQ[i] = -dq;
    }
}

//  Compute the coeffciient matrix of the linearized set of equations for heads.
void AWHSolver::setSteadyMatrixCoeffs()
{
    memset(&xQ[0], 0, nodeCount * sizeof(double));
    matrixSolver->reset();
    setSteadyLinkCoeffs();
    setSteadyNodeCoeffs();
    setSteadyValveCoeffs();
}

//  Compute matrix coefficients for link head loss gradients.
void AWHSolver::setSteadyLinkCoeffs()
{
    // CRITICAL: Double-check that our linkCount matches reality
    if (linkCount != network->links.size()) {
        network->msgLog << "\nWARNING: linkCount mismatch in setSteadyLinkCoeffs() - " 
                       << linkCount << " vs " << network->links.size();
        // Update to correct value
        linkCount = network->links.size();
    }
    
    // Clear node flow balances
    for (int i = 0; i < nodeCount; i++) {
        if (i < xQ.size()) {
            xQ[i] = 0.0;
        }
    }

    for (int j = 0; j < linkCount; j++)
    {
        Link* link = network->links[j];

        // Skip links with zero head gradient.
        if (link->hGrad == 0.0) continue;

        // Identify end nodes and their indices.
        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        int n1 = node1->index;
        int n2 = node2->index;

        // Update node flow balances.
        xQ[n1] -= link->flow;
        xQ[n2] += link->flow;

        // Compute coefficients: a is the link contribution factor, and b is head loss.
        double a = 1.0 / link->hGrad;
        double b = a * link->hLoss;

        //double b = a * (link->hLoss - link->hGrad * link->flow);

        // Update off-diagonal coefficient if both nodes are not fixed.
        if (!node1->fixedGrade && !node2->fixedGrade) {
            matrixSolver->addToOffDiag(j, -a);
        }

        // For the start node:
        if (node1->fixedGrade) {
            matrixSolver->addToRhs(n2, a * node1->head);
        }
        else {
            matrixSolver->addToDiag(n1, a);
            matrixSolver->addToRhs(n1, b);
        }

        // For the end node:
        if (node2->fixedGrade) {
            matrixSolver->addToRhs(n1, a * node2->head);
        }
        else {
            matrixSolver->addToDiag(n2, a);
            matrixSolver->addToRhs(n2, -b);
        }
    }
}

//  Compute matrix coefficients for dynamic tanks and external node outflows.

void  AWHSolver::setSteadyNodeCoeffs()
{
    for (int i = 0; i < nodeCount; i++)
    {
        // ... if node's head not fixed

        Node* node = network->nodes[i];
        if (!node->fixedGrade)
        {
            // ... for dynamic tanks, add area terms to row i
            //     of the head solution matrix & r.h.s. vector

            if (node->type() == Node::TANK && theta != 0.0)
            {
                Tank* tank = static_cast<Tank*>(node);
                double a = tank->area / (theta * tstep);
                matrixSolver->addToDiag(i, a);

                a = a * tank->pastHead + (1.0 - theta) * tank->pastOutflow / theta;
                matrixSolver->addToRhs(i, a);
            }

            // ... for junctions, add effect of external outflows

            else if (node->type() == Node::JUNCTION)
            {
                // ... update junction's net inflow
                xQ[i] -= node->outflow;
                matrixSolver->addToDiag(i, node->qGrad);
                matrixSolver->addToRhs(i, node->qGrad * node->head);
            }

            // ... add node's net inflow to r.h.s. row
            matrixSolver->addToRhs(i, (double)xQ[i]);
        }

        // ... if node has fixed head, force solution to produce it

        else
        {
            matrixSolver->setDiag(i, 1.0);
            matrixSolver->setRhs(i, node->head);
        }
    }
}

//  Compute matrix coefficients for pressure regulating valves.

void  AWHSolver::setSteadyValveCoeffs()
{
    for (Link* link : network->links)
    {
        // ... skip links that are not active pressure regulating valves

        if (link->hGrad > 0.0) continue;

        // ... determine end node indexes of link

        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        // ... add net inflow of downstream node of a PRV to the
        //     r.h.s. row of its upstream node

        if (link->isPRV())
        {
            matrixSolver->addToRhs(n1, (double)xQ[n2]);
        }

        // ... add net inflow of upstream node of a PSV to the
        //     r.h.s. row of its downstream node

        if (link->isPSV())
        {
            matrixSolver->addToRhs(n2, (double)xQ[n1]);
        }
    }
}

int AWHSolver::determineSolverMode(SolverMode currentMode, double* currentTime)
{
    // Threshold constants (can be static members or file-scope globals)
    static const double PHI_A_D = 0.5;  // “dynamic” threshold for compressibility
    static const double PHI_R_D = 0.25; // “dynamic” threshold for relative
    static const double PHI_A_I = 0.1;  // “inertial” threshold
    static const double PHI_R_I = 0.05; // “inertial” threshold

	double current_time = *currentTime;

    // 1) Gather maxPhiA, maxPhiR over all links
    double maxPhiA = 0.0;
    double maxPhiR = 0.0;
    for (Link* link : network->links) {
        if (!link) continue;
        if (link->phiA > maxPhiA) maxPhiA = link->phiA;
        if (link->phiR > maxPhiR) maxPhiR = link->phiR;
    }

    if (current_time <= 1e-12) {
        return 0;  // 0 => QUASI_STEADY
    }

    if ((maxPhiA < PHI_A_I) || (maxPhiR < PHI_R_I)) {
        return 0; // 0 => QUASI_STEADY
    }

    if (currentMode == 2) {
        // If below “dynamic” thresholds => switch to RWC
        // (the figure shows: if maxPhiA < phiA_D or maxPhiR < phiR_D => RWC)
        if ((maxPhiA < PHI_A_D) || (maxPhiR < PHI_R_D)) {
            return 1; // 1 => RWC
        }
        else {
            // stay Water Hammer
            return 2;
        }
    }
    //    If we are currently in RWC (1):
    else if (currentMode == 1) {
        // If above “dynamic” thresholds => switch to Water Hammer
        if ((maxPhiA > PHI_A_D) && (maxPhiR > PHI_R_D)) {
            return 2; // Water Hammer
        }
        else {
            // remain RWC
            return 1;
        }
    }

    if ((maxPhiA > PHI_A_D) && (maxPhiR > PHI_R_D)) {
        return 2;
    }
    else {
        return 1;
    }
}

void AWHSolver::switchToWaterHammerMode(double currentTime) {
    network->msgLog << "\nSwitching to water hammer mode at time " << currentTime;
    
    // Calculate initialization parameters locally
    double maxWaveTime = calculateMaxWaveTravel();
    double safetyMargin = 5.0 * network->option(Options::HYD_STEP);
    double requiredHistoryDuration = maxWaveTime + safetyMargin;
    double historyStartTime = currentTime - requiredHistoryDuration;
    
    // Update mode flags
    waterHammerMode = true;
    currentMode = WATER_HAMMER;
    
    network->msgLog << "\n  Water hammer mode activated";
}

double AWHSolver::calculateMaxWaveTravel() {
    double maxWaveTime = 0.0;
    
    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (!pipe) continue;
        
        double waveSpeed = pipe->getWaveSpeed();
        if (waveSpeed > 0) {
            double pipeWaveTime = pipe->length / waveSpeed;
            maxWaveTime = std::max(maxWaveTime, pipeWaveTime);
        }
    }
    
    return maxWaveTime;
}

void AWHSolver::switchToStandardMode() {
    if (!waterHammerMode) return;  // Already in standard mode
    
    network->msgLog << "\nSwitching back to Standard mode from Water Hammer mode";
    
    // Reset flags
    waterHammerMode = false;
}

// Check convergence of the solution
bool AWHSolver::hasConverged() {
    return
        (hydBalance.maxHeadErr < headErrLimit) &&
        (hydBalance.maxFlowErr < flowErrLimit) &&
        (hydBalance.maxFlowChange < flowChangeLimit) &&
        (hydBalance.totalFlowChange < flowRatioLimit);
}

void AWHSolver::handleCompressibleFlowCoeffs(Link* link, Node* node1, Node* node2, double currentTime, int n1, int n2) {
    // Skip if link isn't a pipe
    Pipe* pipe = dynamic_cast<Pipe*>(link);
    if (!pipe) return;        
    
    // === CALCULATE CHARACTERISTIC PARAMETERS ===
    
    // Physical constants
    double c = pipe->getWaveSpeed();
    double A = pipe->getArea();
    double L = pipe->length;
    double g = GRAVITY;
    double Bj = c / (g * A);  // Wave impedance term (B')
    double D = pipe->diameter;
    double epsilon = 0;  // Friction integration parameter
    
   // Calculate wave travel time and reach-back time
    double physicalWaveTime = L / c;
            
    // Determine reach-back index based on number of reaches
    int reachBackSteps = pipe->numReaches;
    double waveTravel = reachBackSteps * tstep;  // Wave travel distance

    // Retrieve reach-back flows and heads from history
    double Qa = pipe->pastStartFlow;  // Default fallback
    double Qb = pipe->pastEndFlow;    // Default fallback
    double Ha = pipe->fromNode->pastHead;  // Default fallback
    double Hb = pipe->toNode->pastHead;    // Default fallback
            
    // Use the FlowHistoryManager to get reach-back values
    FlowHistoryResult historyResult = FlowHistoryManager::getInstance().getReachBackValues(pipe, currentTime, physicalWaveTime, network);
            
    // If history was found, use the values from the result
    if (historyResult.found) {
        Qa = historyResult.startFlow;
        Qb = historyResult.endFlow;
        Ha = historyResult.startHead;
        Hb = historyResult.endHead;
    } 

    // Convert gradients to resistance coefficients for characteristic equations
    double KA, KB;

    if (network->option(Options::HEADLOSS_MODEL) == "D-W")
    {
        double frictionFactor = pipe->getFrictionFactor(pipe->getRe(pipe->flow, network->option(Options::KIN_VISCOSITY)));
        double K = frictionFactor * pipe->length / (2 * g * D * A * A);

        double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(Qa, network->option(Options::KIN_VISCOSITY)));
        double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(Qb, network->option(Options::KIN_VISCOSITY)));
                // Split K for each end as needed (KA for upstream, KB for downstream)
        KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
        KB = frictionFactorB * pipe->length / (2 * g * D * A * A);

    }
    else
    {
        KA = pipe->resistance;
        KB = pipe->resistance;
    }

    // === ENSURE NUMERICAL STABILITY ===
    const double MIN_RESISTANCE = 1e-10;
    KA = std::max(KA, MIN_RESISTANCE);
    KB = std::max(KB, MIN_RESISTANCE);
    
    // === UNSTEADY FRICTION TERMS ===
    double kUFa = pipe->computeShearDecayCoeff(Qa);
    double kUFb = pipe->computeShearDecayCoeff(Qb);
    double Ua = 0.0;
    double Ub = 0.0;
    
    if (kUFa > 0.0 && kUFb > 0.0) {
        auto sign = [](double val) -> double { 
            return (val > 0.0) ? 1.0 : ((val < 0.0) ? -1.0 : 0.0); 
        };
        
        // CORRECTED: Unsteady friction terms
        Ua = 0.5 * Bj * kUFb * (Qb - 2.0 * sign(Qa) * std::abs(Qb - Qa));
        Ub = 0.5 * Bj * kUFa * (Qa - 2.0 * sign(Qb) * std::abs(Qa - Qb));
    }
    
    // CORRECTED: Characteristic equation coefficients
    double BA, BB, RA, RB;

    if (network->option(Options::HEADLOSS_MODEL) == "D-W") {
        // For boundary A: using C+ from point b
        BA = Bj + epsilon * KB * std::abs(Qb) + 0.5 * Bj * kUFb;
        RA = (Bj - (1.0 - epsilon) * KB * std::abs(Qb)) * Qb - Hb + Ub;

        // For boundary B: using C- from point a
        BB = Bj + epsilon * KA * std::abs(Qa) + 0.5 * Bj * kUFa;
        RB = (Bj - (1.0 - epsilon) * KA * std::abs(Qa)) * Qa + Ha + Ua;
    }
    else if (network->option(Options::HEADLOSS_MODEL) == "H-W") {
        // For boundary A: using C+ from point b
        BA = Bj + epsilon * KB * pow(std::abs(Qb), 0.852) + 0.5 * Bj * kUFb;
        RA = (Bj - (1.0 - epsilon) * KB * pow(std::abs(Qb), 0.852)) * Qb - Hb + Ub;
        
        // For boundary B: using C- from point a
        BB = Bj + epsilon * KA * pow(std::abs(Qa), 0.852) + 0.5 * Bj * kUFa;
        RB = (Bj - (1.0 - epsilon) * KA * pow(std::abs(Qa), 0.852)) * Qa + Ha + Ua;

    }

    // === APPLY WAVE ATTENUATION SCHEME IF ENABLED ===
    if (pipe->impingingWaveCount > 0) {
        //updateWaveAttenuationValues(pipe, BA, BB, RA, RB);
    }
    
    // Store for reference and debugging
    pipe->charB_plus = BA;  // Coefficient for boundary A (using C+)
    pipe->charB_minus = BB; // Coefficient for boundary B (using C-)
    pipe->charC_plus = RA;  // Reach-back term for boundary A
    pipe->charC_minus = RB; // Reach-back term for boundary B

    // === APPLY BOUNDARY CONDITIONS TO MATRIX ===

    if (!node1->fixedGrade) 
    {
        double diagContribA = 1.0 / BA;
        
        matrixSolver->addToRhs(n1, -RA / BA);
        matrixSolver->addToDiag(n1, 1.0 / BA);
    }

    if (!node2->fixedGrade) 
    {
        double diagContribB = 1.0 / BB;
        
        matrixSolver->addToRhs(n2, RB / BB);
        matrixSolver->addToDiag(n2, 1.0 / BB);
    }
}

// Handle quasi-steady flow
void AWHSolver::handleQuasiSteadyFlowCoeffs(Link* link, Node* node1, Node* node2, int n1, int n2)
{
    // Calculate coefficient a = 1/hGrad
    double a = 1.0 / (link->hGrad) ;
	double b = a * (link->hLoss); // - link->hGrad * link->flow);

    //double b = a * ((link->hLoss) + ((1 - kappa) / kappa) * (link->pastHloss - (node1->pastHead - node2->pastHead)) * (link->pastFlow) - link->hGrad * link->flow); // + link->flow;

    if (link->hLoss == 0.0) {
        b = MIN_GRADIENT; // Avoid division by zero if hLoss is zero
    }
    // Process upstream node (node1)
    if (!node1->fixedGrade) {
        // Add to diagonal
        matrixSolver->addToDiag(n1, a);
        
        // Add to RHS
        matrixSolver->addToRhs(n1, b);
    } else {
        matrixSolver->addToRhs(n2, a * node1->head);
    }
    
    // Process downstream node (node2)
    if (!node2->fixedGrade) {
        // Add to diagonal
        matrixSolver->addToDiag(n2, a);
        
        // Add to RHS
        matrixSolver->addToRhs(n2, -b);
    } else {
        matrixSolver->addToRhs(n1, a * node2->head);
    }
    
    // Add off-diagonal terms if both nodes are not fixed
    if (!node1->fixedGrade && !node2->fixedGrade) {
        matrixSolver->addToOffDiag(link->index, -a);
    }
}

void AWHSolver::classifyNodesAndLinks()
{
    // Initialize arrays for tracking compressibility and inertia
    SC.resize(network->links.size(), 0);
    SI.resize(network->links.size(), 0);

    // Log the beginning of classification
    //std::cout << "Classifying links using simplified COMPRESSIBLE/INCOMPRESSIBLE model..." << std::endl;

    // Count for reporting
    int compressibleCount = 0;
    int incompressibleCount = 0;

    // Classify each link based on its dynamic properties
    for (size_t i = 0; i < network->links.size(); i++)
    {
        Link* link = network->links[i];
        if (!link) continue;

        Pipe* pipe = dynamic_cast<Pipe*>(link);

        // Default all links to incompressible flow
        link->flowModel = Link::INCOMPRESSIBLE;
        SC[i] = 0;
        SI[i] = 0;

        // Valves and pumps are always incompressible
        if (link->type() == Link::VALVE || link->type() == Link::PUMP || link->status == Link::LINK_CLOSED)
        {
            link->flowModel = Link::INCOMPRESSIBLE;
            incompressibleCount++;
            continue;
        }
        else
        {
            link->flowModel = Link::COMPRESSIBLE;
            compressibleCount++;    
        }
    }
    // Log classification results
    std::cout << "Classification complete: "
              << compressibleCount << " compressible links, "
              << incompressibleCount << " incompressible links" << std::endl;
}

void AWHSolver::updateSolution(double lamda, bool useShadowMode) {
    // First, ensure vectors are properly sized
    if (currentHeads.size() != nodeCount) {
        currentHeads.resize(nodeCount, 0.0);
    }

    if (currentFlows.size() != linkCount) {
        currentFlows.resize(linkCount, 0.0);
    }

    // Update nodes
    for (int i = 0; i < nodeCount && i < dH.size(); i++) {
        if (network->nodes[i]) {
            if (useShadowMode) {
                // Update shadow values only
                shadowHeads[i] += lamda * dH[i];
                currentHeads[i] = shadowHeads[i]; // For flow indicators
            } else {
                // Normal update
                network->nodes[i]->head += lamda * dH[i];
                currentHeads[i] = network->nodes[i]->head;
            }
        }
    }

    // Update links
    for (int i = 0; i < linkCount && i < dQ.size(); i++) {
        if (network->links[i]) {
            Link* link = network->links[i];
            
            if (useShadowMode) {
                // Update shadow flows only
                shadowFlows[i] += lamda * dQ[i];
                currentFlows[i] = shadowFlows[i];
            } else {
                // Normal updates based on flow model
                if (link->flowModel == Link::COMPRESSIBLE) {
                    Pipe* pipe = dynamic_cast<Pipe*>(link);
                    if (pipe) {
                        pipe->startFlow += lamda * dQ[i];
                        pipe->endFlow += lamda * dQ[i];
                        currentFlows[i] = 0.5 * (pipe->startFlow + pipe->endFlow);
                    } else {
                        link->flow += lamda * dQ[i];
                        currentFlows[i] = link->flow;
                    }
                } else {
                    link->flow += lamda * dQ[i];
                    Pipe* pipe = dynamic_cast<Pipe*>(link);
                    if (pipe) {
                        pipe->startFlow = link->flow;
                        pipe->endFlow = link->flow;
                    }
                    currentFlows[i] = link->flow;
                }
            }
        }
    }
}

void AWHSolver::updateWHSolution(double lamda) {
    // Track maximum solution changes for monitoring convergence
    double maxHeadChange = 0.0;
    double maxFlowChange = 0.0;
    
    // Resize solution tracking vectors if needed
    if (currentHeads.size() != nodeCount) currentHeads.resize(nodeCount);
    if (currentFlows.size() != linkCount) currentFlows.resize(linkCount);

    // ==== STEP 1: UPDATE NODE HEADS ====
    for (int i = 0; i < nodeCount && i < dH.size(); i++) {
        Node* node = network->nodes[i];
        if (!node) continue;  // Skip null nodes
        
        // Skip fixed-grade nodes (reservoirs) - their heads don't change
        if (node->fixedGrade) continue;

        // Store previous head for tracking changes
        double prevHead = node->head;
        
        // Apply change with lambda damping factor
        node->head += lamda * dH[i];
        
        // Track maximum head change
        double headChange = std::abs(node->head - prevHead);
        if (headChange > maxHeadChange) {
            maxHeadChange = headChange;
        }
        
        // Store updated value
        currentHeads[i] = node->head;
    }

    // ==== STEP 2: UPDATE LINK FLOWS ====
    // First update compressible links (water hammer pipes)
    for (int i = 0; i < linkCount; i++) {
        Link* link = network->links[i];
        if (!link) continue;  // Skip null links
        
        // ==== AWH HANDLING FOR PIPES ====
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (pipe && pipe->flowModel == Link::COMPRESSIBLE) {
            // Store previous values for stability control and change tracking
            double prevStartFlow = pipe->startFlow;
            double prevEndFlow = pipe->endFlow;
            double prevFlow = pipe->flow;
            
            // Special handling for reservoir boundary conditions
            bool hasReservoirStart = (pipe->fromNode->type() == Node::RESERVOIR);
            bool hasReservoirEnd = (pipe->toNode->type() == Node::RESERVOIR);
            
            // If connected to reservoirs, recalculate flows directly from compatibility equations
            if (hasReservoirStart || hasReservoirEnd) {
                // Get characteristic parameters
                double BA = pipe->charB_minus;
                double BB = pipe->charB_plus;
                double RA = pipe->charC_minus;
                double RB = pipe->charC_plus;
                
                // For upstream reservoir - directly determine startFlow
                if (hasReservoirStart) {
                    double fixedHead = pipe->fromNode->head;
                    pipe->startFlow += lamda * dQ_start[i]; // = (fixedHead + RA) / BA;
                } else {
                    // Normal update for startFlow
                    pipe->startFlow += lamda * dQ_start[i];
                }
                
                // For downstream reservoir - directly determine endFlow
                if (hasReservoirEnd) {
                    double fixedHead = pipe->toNode->head;
                    pipe->endFlow += lamda * dQ_end[i]; // = (-fixedHead + RB) / BB;
                } else {
                    // Normal update for endFlow
                    pipe->endFlow += lamda * dQ_end[i];
                }
            } else {
                // Normal pipes - apply flow changes with lambda damping
                if (i < dQ_start.size()) pipe->startFlow += lamda * dQ_start[i];
                if (i < dQ_end.size()) pipe->endFlow += lamda * dQ_end[i];
            }
            
            // Set reference flow (average for reporting)
            pipe->flow = 0.5 * (pipe->startFlow + pipe->endFlow);
            currentFlows[i] = pipe->flow;
            
            // Track maximum flow change
            double flowChange = std::max(
                std::abs(pipe->startFlow - prevStartFlow),
                std::abs(pipe->endFlow - prevEndFlow)
            );
            if (flowChange > maxFlowChange) {
                maxFlowChange = flowChange;
            }
        }
        // Handle all other links after water hammer pipes
    }
    
    // Now update incompressible flow links
    for (int i = 0; i < linkCount; i++) {
        Link* link = network->links[i];
        if (!link) continue;
        
        // Skip water hammer pipes (already updated)
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (pipe && pipe->flowModel == Link::COMPRESSIBLE) continue;
        
        // Process all other links (incompressible)
        // Store previous value for stability control
        double prevFlow = link->flow;
            
        // Apply changes with lambda factor
        link->flow += lamda * dQ[i];
            
        // Basic stability control
        double maxChange = 0.5 * std::abs(prevFlow) + 0.001;  // Add small constant for zero flows
        
        // For pipes with incompressible flow, set start/end flows equal
        if (pipe) {
            pipe->startFlow = link->flow;
            pipe->endFlow = link->flow;
            pipe->flow = link->flow ;  // Use average flow for reporting
        }
        
        currentFlows[i] = link->flow;
        
        // Track maximum flow change
        double flowChange = std::abs(link->flow - prevFlow);
        if (flowChange > maxFlowChange) {
            maxFlowChange = flowChange;
        }
    }
}

// Dynamically check for link status changes
bool AWHSolver::linksChangedStatus() {
    bool statusChanged = false;

    for (int i = 0; i < linkCount; i++) {
        Link* link = network->link(i);
        if (!link) continue;

        double q = link->flow;                  // Current flow
        double h1 = link->fromNode->head;       // Upstream head
        double h2 = link->toNode->head;         // Downstream head

        // Update link status based on flow and head conditions
        link->updateStatus(q, h1, h2);

        // ... check for flow into full or out of empty tanks

        if ( link->status > Link::LINK_CLOSED )
        {
            if ( link->fromNode->isClosed(q) || link->toNode->isClosed(-q) )
            {
                link->status = Link::TEMP_CLOSED;
                link->flow = ZERO_FLOW;
                
            }
        }

        // Check if the link status changed
        if (link->status != link->previousStatus) {
            statusChanged = true;
            link->previousStatus = link->status;
            if (link->status == Link::TEMP_CLOSED || link->status == Link::LINK_CLOSED) {
                link->flow = ZERO_FLOW;
                dQ[i] = 0.0;
            }
        }

        // Adaptive flow regime switching
        if (phiA[i] > phiATolerance && phiR[i] > phiRTolerance) {
            SC[i] = 1;  // Enable compressible model
        }
        else if (phiA[i] > inertiaTolerance) {
            SI[i] = 1;  // Enable inertial model
        }
        else {
            SC[i] = 0;  // Default to quasi-steady model
            SI[i] = 0;
        }
    }

    return statusChanged;
}

// Report the results of each trial
void AWHSolver::reportTrial(int trials, double stepSize) {
	std::cout << "Trial " << trials << ": Step size = " << stepSize << std::endl;
	std::cout << "Total error norm = " << errorNorm << std::endl;
}


// Find the error norm for a given step size
double AWHSolver::findErrorNorm(double lamda, double currentTime, double tstep) {
	
    // Safety check for array sizes
    if (xQ.size() < nodeCount || dH.size() < nodeCount || dQ.size() < linkCount) {
        std::cerr << "Error: Vector sizes too small in findErrorNorm" << std::endl;
        return std::numeric_limits<double>::max();
    }

    hLossEvalCount++;
    return hydBalance.evaluate(lamda, (double*)&dH[0], (double*)&dQ[0],
        (double*)&xQ[0], network, currentTime, tstep, currentMode);
}

// Find the optimal step size for the next trial
double AWHSolver::findStepSize(int trials, double currentTime) {
	double lamda = 1.0;
	errorNorm = findErrorNorm(lamda, currentTime, tstep);
	return errorNorm;
}

double AWHSolver::adjustTimeStep(double* currentTime, bool significantHydraulicEvent, bool waterHammerMode) {
    //std::cout << "Adjusting time step with virtual segmentation..." << std::endl;

    double current_Time = *currentTime;
    double tstep = 0.0; // Default value

    // Use the new function to get recommended time step
    double recommendedStep = recommendedTimestep(*currentTime, significantHydraulicEvent, waterHammerMode);

    return recommendedStep;
}

void AWHSolver::setFixedGradeNodes()
{
    // Loop through each active link to check for PRV/PSV valves.
    for (Link* link : network->links)
    {
        Node* node = nullptr;  // Declare a local pointer for this iteration.

        // For a PRV, the control node is the toNode;
        // for a PSV, the control node is the fromNode.
        if (link->isPRV())
            node = link->toNode;
        else if (link->isPSV())
            node = link->fromNode;
        else
            continue;  // Skip other types.

        // Set the fixed grade status based on link status.
        if (link->status == Link::VALVE_ACTIVE)
        {
            node->fixedGrade = true;
            node->head = link->setting + node->elev;
        }
        else {
            node->fixedGrade = false;
        }
    }

    // After time 0, if theta and tstep are nonzero, ensure that tank nodes are not fixed.
    if (theta > 0.0 && tstep > 0.0)
    {
        for (Node* tankNode : network->nodes)
        {
            if (tankNode->type() == Node::TANK)
                tankNode->fixedGrade = false;
        }
    }
}

// Sign function helper
int AWHSolver::sign(double val) {
    return (val > 0.0) ? 1 : ((val < 0.0) ? -1 : 0);
}







    
    