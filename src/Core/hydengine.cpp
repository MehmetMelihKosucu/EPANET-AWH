/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

 //////////////////////////////////////////////
 //  Implementation of the HydEngine class.  //
 //////////////////////////////////////////////

 // TO DO:
 // - add support for Rule-Based controls
 // - add support for saving hydraulics results to a binary file

#include "hydengine.h"
#include "network.h"
#include "error.h"
#include "Solvers/hydsolver.h"
#include "Solvers/matrixsolver.h"
#include "link.h"
#include "Elements/tank.h"
#include "Elements/pattern.h"
#include "Elements/control.h"
#include "Utilities/utilities.h"
#include "awhsolver.h"
#include "pipe.h"
#include "project.h"
#include "Solvers/awhsolver.h"

#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>
#include <vector>
using namespace std;

//static const string s_Balancing  = " Balancing the network:";
static const string s_Unbalanced =
    "  WARNING - network is unbalanced. Flows and pressures may not be correct.";
static const string s_UnbalancedHalted =
    "  Network is unbalanced. Simulation halted by user.";
static const string s_IllConditioned =
    "  Network is numerically ill-conditioned. Simulation halted.";
static const string s_Balanced   = "  Network balanced in ";
static const string s_Trials     = " trials.";
static const string s_Deficient  = " nodes were pressure deficient.";
static const string s_ReSolve1   =
    "\n    Re-solving network with these made fixed grade.";
static const string s_Reductions1 =  " nodes require demand reductions.";
static const string s_ReSolve2    =
    "\n    Re-solving network with these reductions made.";
static const string s_Reductions2 =
    " nodes require further demand reductions to 0.";

//-----------------------------------------------------------------------------

//  Constructor

HydEngine::HydEngine() : 
    engineState(HydEngine::CLOSED),
    network(nullptr),
    hydSolver(nullptr),
    matrixSolver(nullptr),
    saveToFile(false),
    halted(false),
    startTime(0),
    rptTime(0),
    hydStep(0),
    currentTime(0),
    timeOfDay(0),
    peakKwatts(0.0)
{
}

//-----------------------------------------------------------------------------

//  Destructor

HydEngine::~HydEngine()
{
    close();

//    cout << "\nHydEngine destructed.\n";
}

//-----------------------------------------------------------------------------

//  Opens the hydraulic engine.

void HydEngine::open(Network* nw)
{
    // Close a currently opened engine if necessary
    if (engineState != HydEngine::CLOSED) close();

    // Assign the main network pointer
    network = nw;

    // 1) Create sub-models in the main network
    network->createHeadLossModel();
    network->createDemandModel();
    network->createLeakageModel();


    // 3) Create and initialize a matrix solver
    matrixSolver = MatrixSolver::factory(
        network->option(Options::MATRIX_SOLVER),
        network->msgLog
    );
    if (matrixSolver == nullptr)
    {
        throw SystemError(SystemError::MATRIX_SOLVER_NOT_OPENED);
    }
    initMatrixSolver();

    // 4) Create a hydraulic solver
    hydSolver = HydSolver::factory(
        network->option(Options::HYD_SOLVER),
        network,
        matrixSolver
    );
    if (hydSolver == nullptr)
    {
        throw SystemError(SystemError::HYDRAULIC_SOLVER_NOT_OPENED);
    }

    // 5) Mark engine as open
    engineState = HydEngine::OPENED;
}

/*void HydEngine::open(Network* nw)
{
    // Close an already opened engine
    if (engineState != HydEngine::CLOSED) close();
    network = nw;

    // Initialize network segments
    network->initializeSegments();

    // Create hydraulic sub-models (can throw exceptions)
    network->createHeadLossModel();
    network->createDemandModel();
    network->createLeakageModel();

    // Calculate total number of nodes including interior points
    int totalNodes = network->nodes.size() + network->getInteriorNodes().size();
    int totalLinks = network->getTotalSegmentCount();

    std::cout << "Initializing hydraulic engine:\n"
        << "Total nodes (including interior): " << totalNodes << "\n"
        << "Total segments: " << totalLinks << std::endl;

    // Create and initialize matrix solver with enhanced dimensions
    matrixSolver = MatrixSolver::factory(
        network->option(Options::MATRIX_SOLVER), network->msgLog);
    if (!matrixSolver) {
        throw SystemError(SystemError::MATRIX_SOLVER_NOT_OPENED);
    }

    initMatrixSolver();

    // Create hydraulic solver
    hydSolver = HydSolver::factory(
        network->option(Options::HYD_SOLVER), network, matrixSolver);
    if (!hydSolver) {
        throw SystemError(SystemError::HYDRAULIC_SOLVER_NOT_OPENED);
    }

    engineState = HydEngine::OPENED;
} // */

//-----------------------------------------------------------------------------

//  Initializes the hydraulic engine.

void HydEngine::init(bool initFlows)
{
    if (engineState == HydEngine::CLOSED) return;

    // Initialize regular network elements first
    for (Link* link : network->links)
    {
        link->initialize(initFlows);
        link->setResistance(network);
    }

    for (Node* node : network->nodes)
    {
        node->initialize(network);
    }

    int patternStep = network->option(Options::PATTERN_STEP);
    int patternStart = network->option(Options::PATTERN_START);
    for (Pattern* pattern : network->patterns)
    {
        pattern->init(patternStep, patternStart);
    }

    // ====== CRITICAL FIX: Proper Flow History Initialization ======
    
    // Get the actual simulation parameters
    startTime = network->option(Options::START_TIME);
    hydStep = network->option(Options::HYD_STEP); // Make sure this gets the actual time step
    currentTime = 0; // This is simulation time, separate from absolute time
    
    // Calculate how much historical data we need
    double maxWaveTime = calculateMaxWaveTravel();
    double safetyMargin = 5.0 * hydStep; // Add safety margin
    double requiredHistoryDuration = maxWaveTime + safetyMargin;
    
    // Calculate when history should start (before simulation start)
    double historyStartTime = currentTime - requiredHistoryDuration;
    
    // Log this information for debugging
    network->msgLog << "\n=== Flow History Initialization ===";
    network->msgLog << "\nSimulation start time: " << currentTime;
    network->msgLog << "\nMax wave travel time: " << maxWaveTime << " seconds";
    network->msgLog << "\nHistory start time: " << historyStartTime;
    network->msgLog << "\nTime step: " << hydStep << " seconds";
    
    // Initialize flow history with proper time alignment
    FlowHistoryManager::getInstance().initializeHistory(
        network, 
        currentTime,      // Current simulation time
        hydStep,          // Time step size
        historyStartTime  // How far back to initialize history
    );

    // Initialize other engine parameters
    halted = 0;
    rptTime = network->option(Options::REPORT_START);
    peakKwatts = 0.0;
    engineState = HydEngine::INITIALIZED;
    timeStepReason = "";
}


//-----------------------------------------------------------------------------

//  Solves network hydraulics at the current point in time.

int HydEngine::solve(double* t) {
    
    if (engineState != HydEngine::INITIALIZED) return 0;

    // Log current time and conditions
    if (network->option(Options::REPORT_STATUS)) {
        network->msgLog << "\n  Hour " << Utilities::getTime(currentTime) << timeStepReason;
    }

    *t = currentTime;

    timeOfDay = static_cast<int>((currentTime + startTime)) % 86400;

    // Update current hydraulic conditions
    updateCurrentConditions();

    if (network->option(Options::REPORT_TRIALS))  network->msgLog << endl;
    int trials = 0;
    int statusCode = hydSolver->solve(hydStep, trials, currentTime);

    if (statusCode == HydSolver::SUCCESSFUL && isPressureDeficient())
    {
        statusCode = resolvePressureDeficiency(trials);
    }

    for (Link* link : network->links)
    {
        link->previousStatus = link->status;
    }

    // Check if this is a new time point (not a re-solve of the same time)
    if (std::abs(currentTime - lastHistoryUpdateTime) > 1e-9) {
        updateFlowHistory(currentTime);
        lastHistoryUpdateTime = currentTime;

        if (network->option(Options::REPORT_STATUS)) {
            network->msgLog << "\nUpdated flow history at time " << currentTime;
        }
    } else {
        // This is a re-solve at the same time point
        if (network->option(Options::REPORT_STATUS)) {
            network->msgLog << "\nSkipped history update (re-solve) at time " << currentTime;
        }
    }

    reportDiagnostics(statusCode, trials);
    if (halted) throw SystemError(SystemError::HYDRAULICS_SOLVER_FAILURE); // */
    return statusCode;

}

double HydEngine::calculateMaxWaveTravel() {
    double maxWaveTime = 0.0;
    
    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (!pipe) continue;
        
        // Calculate the actual physical wave travel time for this pipe
        double waveSpeed = pipe->getWaveSpeed();
        if (waveSpeed > 0) {
            double pipeWaveTime = pipe->length / waveSpeed;
            maxWaveTime = std::max(maxWaveTime, pipeWaveTime);
        }
    }
    
    return maxWaveTime;
}

// ============ Enhanced updateFlowHistory for Building Real History ============
void HydEngine::updateFlowHistory(double currentTime) {
    // Get reference to the history manager
    FlowHistoryManager& historyManager = FlowHistoryManager::getInstance();

    // We ALWAYS record history, regardless of solver mode. This ensures we have continuous data for reachback calculations.

    double actualTimeStep = hydStep;
    
    // Update history for EVERY pipe, ALWAYS
    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (!pipe) continue;
        
        // Get the current hydraulic state
        // These are the actual values from your hydraulic solver
        double actualStartFlow = pipe->startFlow;
        double actualEndFlow = pipe->endFlow;
        double actualFromHead = pipe->fromNode->head;
        double actualToHead = pipe->toNode->head;
        
        // Basic validation - skip only if data is corrupted
        if (!std::isfinite(actualStartFlow) || !std::isfinite(actualEndFlow) ||
            !std::isfinite(actualFromHead) || !std::isfinite(actualToHead)) {
            network->msgLog << "\n[ERROR] Invalid hydraulic values for pipe " << pipe->name 
                           << " at t=" << currentTime << " - skipping this pipe";
            continue;
        }
        
        // Create the history entry vectors
        std::vector<double> startFlows(1, actualStartFlow);
        std::vector<double> endFlows(1, actualEndFlow);
        std::vector<double> startHeads(1, actualFromHead);
        std::vector<double> endHeads(1, actualToHead);
        std::vector<double> junctionHeads;  // Empty for standard pipes
        
        // Initialize wave points vector - it starts empty
        std::vector<WaveAttenuationPoint> wavePoints;
        
        // CRITICAL FIX: Validate impingingWaveCount before using it
        // Check if we should handle wave attenuation points
        if (pipe->impingingWaveCount > 0 && pipe->impingingWaveCount < 100) {
            // Only process if we have a reasonable number of wave points
            // 100 is a sanity check - adjust based on your system's actual needs
            
            try {
                // Safely resize the vector
                wavePoints.resize(pipe->impingingWaveCount);
                
                // Calculate wave point values
                for (int i = 0; i < pipe->impingingWaveCount; i++) {
                    double position = (i + 1.0) / (pipe->impingingWaveCount + 1.0);
                    
                    // Set position (0 to 1 along the pipe)
                    wavePoints[i].position = position;
                    
                    // Linear interpolation of flow along the pipe
                    wavePoints[i].flow = actualStartFlow * (1.0 - position) + 
                                        actualEndFlow * position;
                    
                    // Linear interpolation of head along the pipe
                    wavePoints[i].head = actualFromHead * (1.0 - position) + 
                                        actualToHead * position;
                    
                    // Initialize amplitude to 0 (can be calculated elsewhere if needed)
                    wavePoints[i].amplitude = 0.0;
                }
            } catch (const std::exception& e) {
                // If resize fails, log the error and continue with empty wave points
                network->msgLog << "\n[ERROR] Failed to create wave points for pipe " 
                               << pipe->name << ": " << e.what();
                network->msgLog << "\n  impingingWaveCount = " << pipe->impingingWaveCount;
                wavePoints.clear();  // Ensure vector is empty
            }
        }
        // If impingingWaveCount is 0 or negative, wavePoints remains empty (which is fine)
        
        // ALWAYS update the history, regardless of mode or wave point success
        historyManager.updateHistoryForPipe(
            pipe, 
            currentTime,
            actualTimeStep,
            startFlows, 
            endFlows, 
            startHeads, 
            endHeads, 
            junctionHeads, 
            wavePoints);
    }
}


//-----------------------------------------------------------------------------

int HydEngine::determineSolverMode(int currentMode, double currentTime)
{
    // Threshold constants (can be static members or file-scope globals)
    static const double PHI_A_D = 0.5;  // “dynamic” threshold for compressibility
    static const double PHI_R_D = 0.25; // “dynamic” threshold for relative
    static const double PHI_A_I = 0.1;  // “inertial” threshold
    static const double PHI_R_I = 0.05; // “inertial” threshold

    // 1) Gather maxPhiA, maxPhiR over all links
    double maxPhiA = 0.0;
    double maxPhiR = 0.0;
    for (Link* link : network->links) {
        if (!link) continue;
        if (link->phiA > maxPhiA) maxPhiA = link->phiA;
        if (link->phiR > maxPhiR) maxPhiR = link->phiR;
    }

    // 2) Front-End check (time = 0 => Quasi-Steady)
    if (currentTime <= 1e-12) {
        return 0;  // 0 => QUASI_STEADY
    }

    // 3) Back-End check:
    //    If below smaller “inertial” thresholds => Quasi-Steady
    if ((maxPhiA < PHI_A_I) && (maxPhiR < PHI_R_I)) {
        return 0; // 0 => QUASI_STEADY
    }

    // 4) Otherwise, we remain in “midway.” Distinguish Water Hammer vs RWC:
    //    If we are currently in Water Hammer (2):
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

    // 5) If we were Quasi-Steady but the dynamic thresholds are above 
    //    or some other corner case => pick whichever is appropriate
    //    e.g. if maxPhiA > PHI_A_D && maxPhiR > PHI_R_D => Water Hammer
    if ((maxPhiA > PHI_A_D) && (maxPhiR > PHI_R_D)) {
        return 2;
    }
    else {
        return 1;
    }
}

//  Advances the simulation to the next point in time.

void HydEngine::advance(double* dt)
{
    tstep = 0;

    if ( engineState != HydEngine::INITIALIZED ) return;

    // ... save current results to hydraulics file
    //if ( saveToFile ) errcode = hydWriter.writeResults(hydState, hydTime);

    // ... if time remains, find time (hydStep) until next hydraulic event

    hydStep = 0;
    double timeLeft = network->option(Options::TOTAL_DURATION) - currentTime;

    if ( halted ) timeLeft = 0;
    if (timeLeft > 0)
    {
        // Get base time step from scheduling/controls
        double scheduledStep = getTimeStep();

        // Get recommended time step based on hydraulic conditions
        AWHSolver* awhSolver = dynamic_cast<AWHSolver*>(hydSolver);
        double recommendedStep = scheduledStep; 

        bool waterHammerModeOn = false;

        if (awhSolver) {
            // Get the current solver mode
            bool significantHydraulicEvent = detectSignificantHydraulicEvent(currentTime, waterHammerModeOn) ;

            bool waterHammerMode = waterHammerModeOn;

            int currentMode = awhSolver->currentMode;
            
            // Get recommended time step based on current hydraulic conditions
            recommendedStep = awhSolver->recommendedTimestep(currentTime, significantHydraulicEvent, currentMode);
        }
        
        // Use the smaller of the two time steps (hydraulic needs may require smaller steps)
        //hydStep = std::min(scheduledStep, recommendedStep);

        hydStep = recommendedStep;

        // Ensure we don't go beyond total duration
        if (hydStep > timeLeft) hydStep = timeLeft;

        // Log the final selected time step
        network->msgLog << "\nSelected time step: " << hydStep;
    }
    *dt = hydStep; //*/

    // ... update energy usage and tank levels over the time step

    updateEnergyUsage();
    updateTanks();

	pastJunction();
	pastLink();

    // ... advance time counters

    currentTime += hydStep;
    if ( currentTime >= rptTime )
    {
        rptTime += network->option(Options::REPORT_STEP);
    }

    // ... advance time patterns

    updatePatterns();
}

bool HydEngine::detectSignificantHydraulicEvent(double currentTime, bool waterHammerModeOn) {
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
                waterHammerModeOn = true;
                return true;
            }
            
            // Or if valve status changed (open/closed)
            if (link->status != link->previousStatus) {
                waterHammerModeOn = true;
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
                waterHammerModeOn = true;
                return true;
            }

            // Or if pump status changed (open/closed)
            if (link->status != link->previousStatus) {
                waterHammerModeOn = true;
                return true;
            } // */
        }

    }
    
    // Check for rapid demand changes at nodes
    // (optional, depending on your system)
    
    return false;
}

//-----------------------------------------------------------------------------

//  Closes the hydraulic solver.

void HydEngine::close()
{
    if ( engineState == HydEngine::CLOSED ) return;
    delete matrixSolver;
    matrixSolver = nullptr;
    delete hydSolver;
    hydSolver = nullptr;
    engineState = HydEngine::CLOSED;

    //... Other objects created in HydEngine::open() belong to the
    //    network object and are deleted by it.
}

//-----------------------------------------------------------------------------

void HydEngine::initMatrixSolver()
{
    int nodeCount = network->count(Element::NODE);
    int linkCount = network->count(Element::LINK);
    try
    {
        // ... place the start/end node indexes of each network link in arrays

        vector<int> node1(linkCount);
        vector<int> node2(linkCount);
        for (int k = 0; k < linkCount; k++)
        {
            node1[k] = network->link(k)->fromNode->index;
            node2[k] = network->link(k)->toNode->index;
        }

        // ...  initialize the matrix solver

        matrixSolver->init(nodeCount, linkCount, (int*)&node1[0], (int*)&node2[0]);
    }
    catch (...)
    {
        throw;
    }
}

//  Updates network conditions at start of current time step.

void HydEngine::updateCurrentConditions()
{
    // Get global multiplier and pattern factor from options.
    double multiplier = network->option(Options::DEMAND_MULTIPLIER);
    double patternFactor = 1.0;
    int p = network->option(Options::DEMAND_PATTERN);
    if (p >= 0)
        patternFactor = network->pattern(p)->currentFactor();
    if (patternFactor < 0)
        patternFactor = 0;

    // Update each active node’s demand and head values.
    // (Assuming that activeNodes is accessible via your Project instance.)
    for (Node* node : network->nodes)
    {
        node->findFullDemand(multiplier, patternFactor);
        node->setFixedGrade();
    } 

    // Update each active link’s hydraulic state.
    for (Link* link : network->links)
    {
        // Apply any control patterns.
        link->applyControlPattern(network->msgLog);
    }

    // Update controls (these may still be defined on the base network).
    for (Control* control : network->controls)
    {
        control->apply(network, currentTime, timeOfDay);
    }
}

bool HydEngine::isPressureDeficient()
{
    int count = 0;
    for (Node* node : network->nodes)
    {
        // ... This only gets evaluated for the CONSTRAINED demand model
        if ( node->isPressureDeficient(network) ) count++;
    }
    if ( count > 0 && network->option(Options::REPORT_TRIALS) )
    {
        network->msgLog << "\n\n    " << count << s_Deficient;
    }
    return (count > 0);
}

//-----------------------------------------------------------------------------

int HydEngine::resolvePressureDeficiency(int& trials)
{
    int trials2 = 0;
    int trials3 = 0;
    int trials4 = 0;
    int count1 = 0;
    int count2 = 0;
    bool reportTrials = ( network->option(Options::REPORT_TRIALS) );

    // ... re-solve network hydraulics with the pressure deficient junctions
    //     set to fixed grade (which occurred in isPressureDeficient())

    if ( reportTrials ) network->msgLog << s_ReSolve1;
    int statusCode = hydSolver->solve(hydStep, trials2, currentTime);
    if ( statusCode == HydSolver::FAILED_ILL_CONDITIONED ) return statusCode;

    // ... adjust actual demands for the pressure deficient junctions

    for (Node* node : network->nodes)
    {
        if ( node->type() == Node::JUNCTION && node->fixedGrade )
        {
            node->actualDemand = min(node->actualDemand, node->fullDemand);
            node->actualDemand = max(0.0, node->actualDemand);
            if ( node->actualDemand < node->fullDemand ) count1++;
            node->fixedGrade = false;
        }
    }

    // ... re-solve once more with the reduced demands at the affected junctions

    if (reportTrials )
    {
        network->msgLog << "\n\n    " << count1 << s_Reductions1;
        network->msgLog << s_ReSolve2;
    }
    statusCode = hydSolver->solve(hydStep, trials3, currentTime);

    // ... check once more for any remaining pressure deficiencies

    for (Node* node : network->nodes)
    {
        if ( node->isPressureDeficient(network) )
        {
            count2++;

            // ... remove fixed grade status set in isPressureDeficient
            //     and make actual demand 0
            node->fixedGrade = false;
            node->actualDemand = 0.0;
        }
    }

     // ... if there are any, then re-solve once more

    if ( count2 > 0 )
    {
        if ( reportTrials )
        {
            network->msgLog << "\n    " << count2 << s_Reductions2;
            network->msgLog << s_ReSolve2 << "\n";
        }
        statusCode = hydSolver->solve(hydStep, trials4, currentTime);
    }

    trials += trials2 + trials3 + trials4;
    return statusCode;
}

//-----------------------------------------------------------------------------

//  Report diagnostics on current hydraulics run.

void HydEngine::reportDiagnostics(int statusCode, int trials)
{
    if ( statusCode == HydSolver::FAILED_ILL_CONDITIONED ||
       ( statusCode == HydSolver::FAILED_NO_CONVERGENCE  &&
         network->option(Options::IF_UNBALANCED) == Options::STOP ))
        halted = true;

    if ( network->option(Options::REPORT_TRIALS) ) network->msgLog << endl;
    if ( network->option(Options::REPORT_STATUS) )
    {
        network->msgLog << endl;
        switch (statusCode)
        {
        case HydSolver::SUCCESSFUL:
            network->msgLog <<	s_Balanced << trials << s_Trials;
            break;
        case HydSolver::FAILED_NO_CONVERGENCE:
            if ( halted ) network->msgLog << s_UnbalancedHalted;
            else          network->msgLog << s_Unbalanced;
            break;
        case HydSolver::FAILED_ILL_CONDITIONED:
            network->msgLog << s_IllConditioned;
            break;
        }
        network->msgLog << endl;
    }
}

//-----------------------------------------------------------------------------

//  Determines the next time step to advance hydraulics.

int HydEngine::getTimeStep()
{
    // ... normal time step is user-supplied hydraulic time step

    string reason ;
    int tstep = network->option(Options::HYD_STEP);
    int n = currentTime / tstep + 1;
    tstep = n * tstep - currentTime;
    timeStepReason = "";

    // ... adjust for time until next reporting period

    int t = rptTime - currentTime;
    if ( t > 0 && t < tstep )
    {
        tstep = t;
        timeStepReason = "";
    }

    // ... adjust for time until next time pattern change

    tstep = timeToPatternChange(tstep);

    // ... adjust for shortest time to fill or drain a tank

    tstep = timeToCloseTank(tstep);

    // ... adjust for shortest time to activate a simple control

    tstep = timeToActivateControl(tstep);
    return tstep;
}

//-----------------------------------------------------------------------------

//  Finds shortest time until next change for all time patterns.

int HydEngine::timeToPatternChange(double tstep)
{
    Pattern* changedPattern = nullptr;
    for (Pattern* pattern : network->patterns)
    {
        int t = pattern->nextTime(currentTime) - currentTime;
        if ( t > 0 && t < tstep )
        {
            tstep = t;
            changedPattern = pattern;
        }
    }
    if ( changedPattern )
    {
        timeStepReason = "  (change in Pattern " + changedPattern->name + ")";
    }
    return tstep;
}

//-----------------------------------------------------------------------------

//  Finds the shortest time to completely fill or empty all tanks.

int HydEngine::timeToCloseTank(double tstep)
{
    Tank* closedTank = nullptr;
    for (Node* node : network->nodes)
    {
        // ... check if node is a tank

        if ( node->type() == Node::TANK )
        {
            // ... find the time to fill (or empty) the tank

            Tank* tank = static_cast<Tank*>(node);
            int t = tank->timeToVolume(tank->minVolume);
            if ( t <= 0 ) t = tank->timeToVolume(tank->maxVolume);

            // ... compare this time with current time step

            if ( t > 0 && t < tstep )
            {
                tstep = t;
                closedTank = tank;
            }
        }
    }
    if ( closedTank )
    {
        timeStepReason = "  (Tank " + closedTank->name + " closed)";
    }
    return tstep;
}

//-----------------------------------------------------------------------------

//  Finds the shortest time to activate a simple control.

int HydEngine::timeToActivateControl(double tstep)
{
    bool activated = false;
    for (Control* control : network->controls)
    {
        int t = control->timeToActivate(network, currentTime, timeOfDay);
        if ( t > 0 && t < tstep )
        {
            tstep = t;
            activated = true;
        }
    }
    if ( activated ) timeStepReason = "  (control activated)";
    return tstep;
}

//-----------------------------------------------------------------------------

//  Updates energy usage over the current time step.

void HydEngine::updateEnergyUsage()
{
    // ... use a nominal time step of 1 day if running a single period analysis

    int dt = hydStep;
    if ( network->option(Options::TOTAL_DURATION) == 0 ) dt = 86400;
    if ( dt == 0 ) return;

    // ... update energy usage for each pump link over the time step

    double totalKwatts = 0.0;
    for (Link* link : network->links)
    {
        totalKwatts += link->updateEnergyUsage(network, dt);
    }

    // ... update peak energy usage over entire simulation

    peakKwatts = max(peakKwatts, totalKwatts);
}

//-----------------------------------------------------------------------------

//  Updates tank area and volume over the current time step.

void HydEngine::updateTanks()
{
    for (Node* node : network->nodes)
    {
    	if ( node->type() == Node::TANK )
        {
            Tank* tank = static_cast<Tank*>(node);
            tank->pastHead = tank->head;
			tank->ph = tank->head;
            tank->pastVolume = tank->volume;
			tank->pastArea = tank->area;
            tank->pastOutflow = tank->outflow;
            node->fixedGrade = true;
            tank->updateVolume(hydStep);
            tank->updateArea();
        }
    }
}

void HydEngine::pastJunction()
{
    // Save previous head values for all regular nodes
    for (Node* node : network->nodes)
    {
        if (node->type() == Node::JUNCTION)
        {
            node->pastHead = node->head;
            node->ph = node->head;
            node->pastOutflow = node->outflow;
        }
        else if (node->type() == Node::RESERVOIR)
        {
            node->pastHead = node->head;
            node->ph = node->head;
            node->pastOutflow = node->outflow;
        }
    }
}

void HydEngine::pastLink()
{
    for (Link* link : network->links)
    {
        if (link->type() == Link::PIPE || link->type() == Link::PUMP || link->type() == Link::VALVE)
        {
            // Save past values for all links
            link->pastFlow = link->flow;
            link->pastHloss = link->hLoss;
            link->pastSetting = link->setting;

            // Handle pipe-specific data
            if (link->type() == Link::PIPE)
            {
                Pipe* pipe = static_cast<Pipe*>(link);
                
                // Save past flow values for pipe ends
                pipe->pastStartFlow = pipe->startFlow;
                pipe->pastEndFlow = pipe->endFlow;
                
            }
        }
    }
}

//-----------------------------------------------------------------------------

//  Advances all time patterns.

void HydEngine::updatePatterns()
{
    for (Pattern* pattern : network->patterns)
    {
        pattern->advance(currentTime);
    }
}

void HydEngine::computeFlowIndicators() {
    Epanet::Project* project;
    for (size_t i = 0; i < network->links.size(); i++) {
        
        Link* link = network->links[i];
        if (!link) continue;
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (!pipe) {
            // Log the type of link or handle it differently
            std::cerr << "Link is not a pipe, type: " << typeid(*link).name() << std::endl;
            continue; // Skip non-pipe links instead of throwing an exception
        }

        double dh = std::fabs(link->fromNode->head - link->toNode->head);
        double dq = std::fabs(link->flow - link->pastFlow);

        if (pipe->getWaveSpeed() > 0 && pipe->getArea() > 0) {
            double inertia = link->inertialTerm * dq / link->getArea();
            double compressibility = dh / pipe->getWaveSpeed();

            pipe->phiA = inertia + compressibility; // Absolute indicator
            pipe->phiR = pipe->phiA / (pipe->phiA + std::fabs(pipe->hLoss)); // Relative indicator
        }
        else {
            pipe->phiA = 0.0;
            pipe->phiR = 0.0;
        }
    }
}






