/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION */

// dualflowhistory.cpp
#include "dualflowhistory.h"
#include "pipe.h"  // For Pipe class
#include "network.h"  // For Network class
#include "constants.h"  // For Constants class
#include "Models/headlossmodel.h"
#include "Solvers/awhsolver.h"

// DualFlowHistory methods
void DualFlowHistory::clear() {
    timestamps.clear();
    startFlowHistory.clear();
    endFlowHistory.clear();
    startHeadHistory.clear();
    endHeadHistory.clear();
    junctionHeadHistory.clear();
    pointAFlowHistory.clear();
    pointBFlowHistory.clear();
    pointAHeadHistory.clear();
    pointBHeadHistory.clear();
    waveAttenuationHistory.clear();
}

bool DualFlowHistory::hasEnoughHistory(int reachBackSteps) const {
    return reachBackSteps < timestamps.size();
}


/**
 * Initialize history with physically consistent data and proper time alignment
 */
void FlowHistoryManager::initializeHistory(
    Network* network, 
    double currentTime, 
    double timeStep, 
    double historyStartTime
) {
    if (!network) return;
    
    network->msgLog << "\n=== Initializing Flow History ===";
    network->msgLog << "\n  Current time: " << currentTime;
    network->msgLog << "\n  Time step: " << timeStep;
    network->msgLog << "\n  History start: " << historyStartTime;
    
    clear();
    
    // Calculate required history depth
    int historyDepth = static_cast<int>(
        std::ceil((currentTime - historyStartTime) / timeStep)
    ) + 1;
    
    // Initialize each pipe with TRUE steady-state conditions
    for (Link* link : network->links) {
        Pipe* pipe = dynamic_cast<Pipe*>(link);
        if (!pipe) continue;
        
        auto history = std::make_unique<DualFlowHistory>();
        
        // Get the actual steady-state conditions
        double steadyFlow = pipe->flow;
        double steadyFromHead = pipe->fromNode->head;
        double steadyToHead = pipe->toNode->head;
        
        // Ensure minimum flow for numerical stability
        const double MIN_FLOW = 1e-6;
        if (std::abs(steadyFlow) < MIN_FLOW) {
            steadyFlow = (steadyFlow >= 0) ? MIN_FLOW : -MIN_FLOW;
        }
        
        // Fill history with CONSTANT steady-state values
        // No artificial variation - let the physics create real variation
        for (int t = 0; t < historyDepth; t++) {
            double timestamp = historyStartTime + (t * timeStep);
            history->timestamps.push_back(timestamp);
            
            // Use constant values - no artificial variation
            history->startFlowHistory.push_back(std::vector<double>(1, steadyFlow));
            history->endFlowHistory.push_back(std::vector<double>(1, steadyFlow));
            history->startHeadHistory.push_back(std::vector<double>(1, steadyFromHead));
            history->endHeadHistory.push_back(std::vector<double>(1, steadyToHead));
            history->junctionHeadHistory.push_back(std::vector<double>());
        }
        
        historyData[pipe] = std::move(history);
    }
    
    network->msgLog << "\n  Initialized " << historyData.size() 
                   << " pipes with " << historyDepth << " time points each";
}

void FlowHistoryManager::addHistory(Pipe* pipe, const DualFlowHistory& history) {
    if (!pipe) return;
    
    // Use make_unique to create a heap-allocated copy
    historyData[pipe] = std::make_unique<DualFlowHistory>(history);
}

/**
 * Retrieves hydraulic values from a previous time point based on wave travel characteristics
 * 
 * @param pipe The pipe to retrieve history for
 * @param currentTime The current simulation time
 * @param waveTravel The wave travel time to reach back (calculated from pipe length/wave speed)
 * @param network Optional pointer to network for logging warnings
 * @return FlowHistoryResult containing the reach-back values and status flags
 */

 /**
 * Enhanced method to retrieve hydraulic values with better error handling and validation
 */

FlowHistoryResult FlowHistoryManager::getReachBackValues(
    Pipe* pipe, 
    double currentTime, 
    double waveTravel, 
    Network* network) 
{
    FlowHistoryResult result;
    result.found = false;

    // Safety check for null pointer
    if (!pipe) {
        if (network) {
            network->msgLog << "\nError: Null pipe pointer passed to getReachBackValues";
        }
        return result;
    }

    // Check if history data container exists
    if (historyData.empty()) {
        if (network) {
            network->msgLog << "\nError: History data container is empty";
        }
        return result;
    }
    
    // Find the pipe's history
    auto historyIt = historyData.find(pipe);
    if (historyIt == historyData.end()) {
        if (network) {
            network->msgLog << "\nError: No history data found for pipe " << pipe->name;
        }
        return result;
    }

    // Validate the data pointer
    if (!historyIt->second) {
        if (network) {
            network->msgLog << "\nError: History data pointer is null for pipe " << pipe->name;
        }
        return result;
    }
    
    // Access the history object
    const auto& history = *(historyIt->second);
    
    // Safety check for empty history
    if (history.timestamps.empty()) {
        if (network) {
            network->msgLog << "\nError: History exists but contains no data for pipe " << pipe->name;
        }
        return result;
    }
    
    // Validate that flow/head vectors aren't empty
    if (history.startFlowHistory.empty() || history.endFlowHistory.empty() ||
        history.startHeadHistory.empty() || history.endHeadHistory.empty()) {
        if (network) {
            network->msgLog << "\nError: History has timestamps but missing flow/head data for pipe " 
                           << pipe->name;
        }
        return result;
    }

    // CRITICAL: Calculate the PHYSICAL wave travel time
    // This is the actual time it takes for a pressure wave to traverse the pipe
    double waveSpeed = pipe->getWaveSpeed();
    double physicalWaveTravel = pipe->length / waveSpeed;
    
    // Check if the provided waveTravel makes physical sense
    // If it's very different from the physical value, we have a problem
    double deviation = std::abs(waveTravel - physicalWaveTravel) / physicalWaveTravel;
    
    if (deviation > 0.1) {  // More than 10% deviation
        // This is a critical issue - log it for debugging
        if (network) {
            network->msgLog << "\n[CRITICAL] Wave travel time mismatch for pipe " << pipe->name
                           << "\n  Provided waveTravel: " << waveTravel << " seconds"
                           << "\n  Physical wave time: " << physicalWaveTravel << " seconds"
                           << "\n  Pipe length: " << pipe->length << " m"
                           << "\n  Wave speed: " << waveSpeed << " m/s"
                           << "\n  Deviation: " << (deviation * 100) << "%"
                           << "\n  Using physical value instead.";
        }
        waveTravel = physicalWaveTravel;
    }
    
    // Calculate the exact timestamp we want to reach back to
    double targetTime = currentTime - waveTravel;
    
    // Get the time range of our history
    double oldestTime = history.timestamps.front();
    double newestTime = history.timestamps.back();
    
    // Enhanced helper function with better error checking
    auto safeExtractValues = [&](size_t timeIndex) -> bool {
        // First, verify the time index is valid
        if (timeIndex >= history.timestamps.size()) {
            if (network) {
                network->msgLog << "\nError: Time index " << timeIndex 
                               << " out of bounds for pipe " << pipe->name;
            }
            return false;
        }
        
        // Check that all history vectors have data at this time index
        bool hasValidData = true;
        
        // Check each history vector individually for better diagnostics
        if (timeIndex >= history.startFlowHistory.size() || 
            history.startFlowHistory[timeIndex].empty()) {
            hasValidData = false;
            if (network) {
                network->msgLog << "\nError: Missing startFlow data at index " << timeIndex
                               << " for pipe " << pipe->name;
            }
        }
        
        if (timeIndex >= history.endFlowHistory.size() || 
            history.endFlowHistory[timeIndex].empty()) {
            hasValidData = false;
            if (network) {
                network->msgLog << "\nError: Missing endFlow data at index " << timeIndex
                               << " for pipe " << pipe->name;
            }
        }
        
        if (timeIndex >= history.startHeadHistory.size() || 
            history.startHeadHistory[timeIndex].empty()) {
            hasValidData = false;
            if (network) {
                network->msgLog << "\nError: Missing startHead data at index " << timeIndex
                               << " for pipe " << pipe->name;
            }
        }
        
        if (timeIndex >= history.endHeadHistory.size() || 
            history.endHeadHistory[timeIndex].empty()) {
            hasValidData = false;
            if (network) {
                network->msgLog << "\nError: Missing endHead data at index " << timeIndex
                               << " for pipe " << pipe->name;
            }
        }
        
        if (!hasValidData) {
            return false;
        }
        
        // Extract the values (using index 0 since you have single segments)
        result.startFlow = history.startFlowHistory[timeIndex][0];
        result.endFlow = history.endFlowHistory[timeIndex][0];
        result.startHead = history.startHeadHistory[timeIndex][0];
        result.endHead = history.endHeadHistory[timeIndex][0];
        result.found = true;
        
        // Add diagnostic output for pipes connected to J874
        /*if (pipe->toNode && pipe->toNode->name.find("874") != std::string::npos) {
            network->msgLog << "\n[J874 Debug] Retrieved history for pipe " << pipe->name
                           << " at time index " << timeIndex 
                           << " (timestamp: " << history.timestamps[timeIndex] << ")"
                           << "\n  Start Flow: " << result.startFlow
                           << "\n  End Flow: " << result.endFlow
                           << "\n  Start Head: " << result.startHead
                           << "\n  End Head: " << result.endHead;
        } // */
        
        return true;
    };
    
    // Determine which case we're in and extract/interpolate accordingly
    if (targetTime <= oldestTime) {
        // We need older data than we have - this is a problem!
        /*if (network) {
            network->msgLog << "\n[WARNING] Insufficient history for pipe " << pipe->name
                           << "\n  Current time: " << currentTime
                           << "\n  Wave travel: " << waveTravel
                           << "\n  Need data from: " << targetTime
                           << "\n  Oldest available: " << oldestTime
                           << "\n  Using oldest available data (may cause errors!)";
        } // */
        
        if (!safeExtractValues(0)) {
            return result;
        }
    }
    else if (targetTime >= newestTime) {
        
        /*if (network) {
            network->msgLog << "\n[ERROR] Target time too recent for pipe " << pipe->name
                           << "\n  Target time: " << targetTime
                           << "\n  Newest available: " << newestTime
                           << "\n  This suggests a logic error!";
        } // */
        
        size_t newestIndex = history.timestamps.size() - 1;
        if (!safeExtractValues(newestIndex)) {
            return result;
        }
    }
    else {
        // Normal case: interpolate between two time points
        
        // Binary search to find bracketing indices
        size_t lowerIndex = 0;
        size_t upperIndex = history.timestamps.size() - 1;
        
        while (upperIndex - lowerIndex > 1) {
            size_t midIndex = (lowerIndex + upperIndex) / 2;
            if (history.timestamps[midIndex] <= targetTime) {
                lowerIndex = midIndex;
            } else {
                upperIndex = midIndex;
            }
        }
        
        // Check for exact or near-exact match
        const double TIME_TOLERANCE = 1e-9;
        
        if (std::abs(history.timestamps[lowerIndex] - targetTime) < TIME_TOLERANCE) {
            if (!safeExtractValues(lowerIndex)) {
                return result;
            }
        }
        else if (std::abs(history.timestamps[upperIndex] - targetTime) < TIME_TOLERANCE) {
            if (!safeExtractValues(upperIndex)) {
                return result;
            }
        }
        else {
            // Need to interpolate between the two indices
            
            // First verify both indices have valid data
            bool lowerValid = (lowerIndex < history.startFlowHistory.size() && 
                              !history.startFlowHistory[lowerIndex].empty() &&
                              lowerIndex < history.endFlowHistory.size() && 
                              !history.endFlowHistory[lowerIndex].empty() &&
                              lowerIndex < history.startHeadHistory.size() && 
                              !history.startHeadHistory[lowerIndex].empty() &&
                              lowerIndex < history.endHeadHistory.size() && 
                              !history.endHeadHistory[lowerIndex].empty());
                              
            bool upperValid = (upperIndex < history.startFlowHistory.size() && 
                              !history.startFlowHistory[upperIndex].empty() &&
                              upperIndex < history.endFlowHistory.size() && 
                              !history.endFlowHistory[upperIndex].empty() &&
                              upperIndex < history.startHeadHistory.size() && 
                              !history.startHeadHistory[upperIndex].empty() &&
                              upperIndex < history.endHeadHistory.size() && 
                              !history.endHeadHistory[upperIndex].empty());
            
            if (!lowerValid || !upperValid) {
                if (network) {
                    network->msgLog << "\nError: Cannot interpolate - missing data at indices "
                                   << lowerIndex << " or " << upperIndex 
                                   << " for pipe " << pipe->name;
                }
                return result;
            }
            
            // Calculate interpolation factor
            double t1 = history.timestamps[lowerIndex];
            double t2 = history.timestamps[upperIndex];
            double alpha = (targetTime - t1) / (t2 - t1);
            
            // Ensure alpha is in valid range
            alpha = std::max(0.0, std::min(1.0, alpha));
            
            // Perform linear interpolation
            result.startFlow = linearInterpolate(
                history.startFlowHistory[lowerIndex][0], 
                history.startFlowHistory[upperIndex][0], 
                alpha
            );
            result.endFlow = linearInterpolate(
                history.endFlowHistory[lowerIndex][0], 
                history.endFlowHistory[upperIndex][0], 
                alpha
            );
            result.startHead = linearInterpolate(
                history.startHeadHistory[lowerIndex][0], 
                history.startHeadHistory[upperIndex][0], 
                alpha
            );
            result.endHead = linearInterpolate(
                history.endHeadHistory[lowerIndex][0], 
                history.endHeadHistory[upperIndex][0], 
                alpha
            );
            result.found = true;
            
        }
    }
    
    // Final validation
    if (result.found) {
        // Check for non-physical values
        const double MAX_REASONABLE_FLOW = 1000.0;  // Adjust based on your system
        const double MAX_REASONABLE_HEAD = 10000.0;  // Adjust based on your system
        
        if (std::abs(result.startFlow) > MAX_REASONABLE_FLOW ||
            std::abs(result.endFlow) > MAX_REASONABLE_FLOW ||
            std::abs(result.startHead) > MAX_REASONABLE_HEAD ||
            std::abs(result.endHead) > MAX_REASONABLE_HEAD ||
            !std::isfinite(result.startFlow) || !std::isfinite(result.endFlow) ||
            !std::isfinite(result.startHead) || !std::isfinite(result.endHead)) {
            
            result.found = false;
        
        }
    }
    
    return result;
}

/**
 * Linear interpolation helper function
 */
double FlowHistoryManager::linearInterpolate(double v1, double v2, double alpha) {
    return v1 + alpha * (v2 - v1);
}


/**
 * Calculate wave attenuation points with proper head loss model integration
 */
void FlowHistoryManager::calculateWaveAttenuationPoints(
    Pipe* pipe,
    std::vector<WaveAttenuationPoint>& wavePoints,
    double startFlow,
    double endFlow,
    double startHead,
    double endHead,
    const DualFlowHistory& history,
    Network* network)
{
    // Make sure the vector is properly sized
    if (wavePoints.size() != pipe->impingingWaveCount) {
        wavePoints.resize(pipe->impingingWaveCount);
    }
    
    // Calculate and store intermediate wave attenuation points
    for (int i = 0; i < pipe->impingingWaveCount; i++) {
        // Relative position along the pipe (evenly distributed)
        double position = (i + 1.0) / (pipe->impingingWaveCount + 1.0);
        
        // Linear interpolation of current values for initial conditions
        wavePoints[i].flow = startFlow * (1.0 - position) + endFlow * position;
        wavePoints[i].head = startHead * (1.0 - position) + endHead * position;
        wavePoints[i].position = position;
        
        // If this isn't the first time step, perform wave tracking calculations
        if (history.timestamps.size() > 1) {
            // Calculate the reach-back time index based on position
            int reachBackIndex = static_cast<int>(std::round(position * pipe->numReaches));
            if (reachBackIndex < history.timestamps.size() - 1) {
                // Apply advanced wave attenuation calculations
                applyWaveTracking(pipe, wavePoints[i], reachBackIndex, history, network);
            }
        }
    }
}

/**
 * Advanced wave tracking method with proper head loss model integration
 */
void FlowHistoryManager::applyWaveTracking(
    Pipe* pipe,
    WaveAttenuationPoint& point,
    int reachBackIndex,
    const DualFlowHistory& history,
    Network* network)
{
    // Only proceed if we have enough history
    if (reachBackIndex >= history.timestamps.size()) return;
    
    // Get wave characteristic parameters
    double c = pipe->getWaveSpeed();
    double A = pipe->getArea();
    double g = GRAVITY;
    double B = c / (g * A);
    double D = pipe->diameter;
    double epsilon = 0.5;  // Friction integration parameter
    
    // Get past flow values from opposite ends (for wave crossing)
    double pastStartFlow = history.startFlowHistory[reachBackIndex][0];
    double pastEndFlow = history.endFlowHistory[reachBackIndex][0];
    double pastStartHead = history.startHeadHistory[reachBackIndex][0];
    double pastEndHead = history.endHeadHistory[reachBackIndex][0];
    
    // CRITICAL FIX: Calculate resistance coefficients using the same head loss model
    // as the other functions to ensure consistency
    double KA = 0.0;  // Resistance coefficient for upstream end
    double KB = 0.0;  // Resistance coefficient for downstream end
    
    // Get the head loss model type from network options
    std::string headlossModel = network->option(Options::HEADLOSS_MODEL);
    
    // Calculate resistance coefficients based on selected head loss model
    if (headlossModel == "D-W") {
        // Darcy-Weisbach model
        double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(pastStartFlow, network->option(Options::KIN_VISCOSITY)));
        double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(pastEndFlow, network->option(Options::KIN_VISCOSITY)));
        
        KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
        KB = frictionFactorB * pipe->length / (2 * g * D * A * A);
    }
    else if (headlossModel == "H-W") {
        // Hazen-Williams model
        HW_HeadLossModel* hwModel = new HW_HeadLossModel(network->option(Options::KIN_VISCOSITY));
        hwModel->setResistance(pipe);
        
        double headLossA, gradientA, headLossB, gradientB;
        hwModel->findHeadLoss(pipe, pastStartFlow, headLossA, gradientA);
        hwModel->findHeadLoss(pipe, pastEndFlow, headLossB, gradientB);
        
        const double MIN_FLOW = 1e-6;
        KA = gradientA / (2.0 * std::max(std::abs(pastStartFlow), MIN_FLOW));
        KB = gradientB / (2.0 * std::max(std::abs(pastEndFlow), MIN_FLOW));
        
        delete hwModel;
    }
    else if (headlossModel == "C-M") {
        // Chezy-Manning model
        CM_HeadLossModel* cmModel = new CM_HeadLossModel(network->option(Options::KIN_VISCOSITY));
        cmModel->setResistance(pipe);
        
        double headLossA, gradientA, headLossB, gradientB;
        cmModel->findHeadLoss(pipe, pastStartFlow, headLossA, gradientA);
        cmModel->findHeadLoss(pipe, pastEndFlow, headLossB, gradientB);
        
        const double MIN_FLOW = 1e-6;
        KA = gradientA / (2.0 * std::max(std::abs(pastStartFlow), MIN_FLOW));
        KB = gradientB / (2.0 * std::max(std::abs(pastEndFlow), MIN_FLOW));
        
        delete cmModel;
    }
    else {
        // Default to Darcy-Weisbach
        double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(pastStartFlow, network->option(Options::KIN_VISCOSITY)));
        double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(pastEndFlow, network->option(Options::KIN_VISCOSITY)));
        
        KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
        KB = frictionFactorB * pipe->length / (2 * g * D * A * A);
    }
    
    // Ensure numerical stability
    const double MIN_RESISTANCE = 1e-10;
    KA = std::max(KA, MIN_RESISTANCE);
    KB = std::max(KB, MIN_RESISTANCE);
    
    // Calculate the wave characteristic values based on position
    // These formulas implement the GAWH approach for tracking waves
    
    // Positive characteristic (wave from upstream)
    double posCharResistance = (B - (1 - epsilon) * KA * std::abs(pastStartFlow)) * pastStartFlow;
    
    // Negative characteristic (wave from downstream)  
    double negCharResistance = (B - (1 - epsilon) * KB * std::abs(pastEndFlow)) * pastEndFlow;
    
    // Combine characteristics based on position to simulate wave interaction
    double weight = point.position;
    double totalResistance = epsilon * (KA * std::abs(pastStartFlow) + KB * std::abs(pastEndFlow));
    
    point.flow = (posCharResistance + negCharResistance + B * (pastStartHead - pastEndHead)) / 
                (2.0 * B + totalResistance);
    
    // Also update head based on characteristic relationships
    point.head = (1.0 - weight) * (pastStartHead - B * point.flow + posCharResistance / B) + 
                weight * (pastEndHead + B * point.flow - negCharResistance / B);
}

void FlowHistoryManager::clear() {
    historyData.clear();
}

bool FlowHistoryManager::hasHistoryForPipe(Pipe* pipe) const {
    return historyData.find(pipe) != historyData.end();
}


/**
 * Update history for a specific pipe with validation and error handling
 */
void FlowHistoryManager::updateHistoryForPipe(
    Pipe* pipe, 
    double currentTime,
    double actualTimeStep,
    const std::vector<double>& startFlows,
    const std::vector<double>& endFlows,
    const std::vector<double>& startHeads,
    const std::vector<double>& endHeads,
    const std::vector<double>& junctionHeads,
    const std::vector<WaveAttenuationPoint>& wavePoints)
{
    // === Input Validation ===
    if (startFlows.empty() || endFlows.empty() || 
        startHeads.empty() || endHeads.empty()) {
        return;
    }
    
    // === Find or Create History ===
    auto historyIt = historyData.find(pipe);
    
    if (historyIt == historyData.end()) {
        // First time - create new history
        auto newHistory = std::make_unique<DualFlowHistory>();
        
        newHistory->timestamps.push_back(currentTime);
        newHistory->startFlowHistory.push_back(startFlows);
        newHistory->endFlowHistory.push_back(endFlows);
        newHistory->startHeadHistory.push_back(startHeads);
        newHistory->endHeadHistory.push_back(endHeads);
        newHistory->junctionHeadHistory.push_back(junctionHeads);
        
        if (!wavePoints.empty()) {
            newHistory->waveAttenuationHistory.push_back(wavePoints);
        }
        
        historyData[pipe] = std::move(newHistory);
        return;
    }
    
    // === Access Existing History ===
    DualFlowHistory& history = *(historyIt->second);
    
    actualTimeStep = Delta_WH;

    // === Calculate Maximum History Length ===
    double waveSpeed = pipe->getWaveSpeed();
    double waveTime = pipe->length / waveSpeed;
    
    // CRITICAL FIX: Use the actual timestep, not a hardcoded value
    // We need enough history to cover at least 2 wave travel times for safety
    size_t requiredLength = static_cast<size_t>(
        std::ceil(2.0 * waveTime / actualTimeStep)  // 2x for safety margin
    ) + 10;  // Additional buffer
    
    // Set reasonable bounds
    const size_t MIN_HISTORY = 100;
    const size_t MAX_HISTORY = 10000;  // Increase this - memory is cheap, accuracy is precious
    
    size_t maxHistoryLength = std::max(MIN_HISTORY, requiredLength);
    maxHistoryLength = std::min(maxHistoryLength, MAX_HISTORY);
    
    // Add diagnostic output for debugging
    static bool firstTime = true;
    
    // === Check for Duplicate Timestamp ===
    if (!history.timestamps.empty()) {
        double lastTime = history.timestamps.back();
        if (std::abs(lastTime - currentTime) < 1e-9) {
            // Update the last entry instead of adding duplicate
            size_t lastIndex = history.timestamps.size() - 1;
            
            history.startFlowHistory[lastIndex] = startFlows;
            history.endFlowHistory[lastIndex] = endFlows;
            history.startHeadHistory[lastIndex] = startHeads;
            history.endHeadHistory[lastIndex] = endHeads;
            history.junctionHeadHistory[lastIndex] = junctionHeads;
            
            if (!history.waveAttenuationHistory.empty() && 
                history.waveAttenuationHistory.size() > lastIndex) {
                history.waveAttenuationHistory[lastIndex] = wavePoints;
            }
            
            return;
        }
    }
    
    // ============================================
    // SIMPLE FIX: Use rotate or pop_front approach instead of erase
    // ============================================
    
    // Option 1: SIMPLEST - Just clear and rebuild if too large
    if (history.timestamps.size() >= maxHistoryLength) {
        // Keep only the most recent half of the data
        size_t keepCount = maxHistoryLength / 2;
        size_t startIdx = history.timestamps.size() - keepCount;
        
        // Create temporary containers with recent data only
        std::vector<double> newTimestamps(
            history.timestamps.begin() + startIdx, 
            history.timestamps.end()
        );
        std::vector<std::vector<double>> newStartFlows(
            history.startFlowHistory.begin() + startIdx,
            history.startFlowHistory.end()
        );
        std::vector<std::vector<double>> newEndFlows(
            history.endFlowHistory.begin() + startIdx,
            history.endFlowHistory.end()
        );
        std::vector<std::vector<double>> newStartHeads(
            history.startHeadHistory.begin() + startIdx,
            history.startHeadHistory.end()
        );
        std::vector<std::vector<double>> newEndHeads(
            history.endHeadHistory.begin() + startIdx,
            history.endHeadHistory.end()
        );
        std::vector<std::vector<double>> newJunctionHeads(
            history.junctionHeadHistory.begin() + startIdx,
            history.junctionHeadHistory.end()
        );
        
        // Replace old data with new (smaller) data
        history.timestamps = std::move(newTimestamps);
        history.startFlowHistory = std::move(newStartFlows);
        history.endFlowHistory = std::move(newEndFlows);
        history.startHeadHistory = std::move(newStartHeads);
        history.endHeadHistory = std::move(newEndHeads);
        history.junctionHeadHistory = std::move(newJunctionHeads);
        
        // Handle wave attenuation if present
        if (!history.waveAttenuationHistory.empty() && 
            history.waveAttenuationHistory.size() >= history.timestamps.size()) {
            std::vector<std::vector<WaveAttenuationPoint>> newWaveHistory(
                history.waveAttenuationHistory.begin() + startIdx,
                history.waveAttenuationHistory.end()
            );
            history.waveAttenuationHistory = std::move(newWaveHistory);
        }
    }
    
    // === Add New Entry ===
    history.timestamps.push_back(currentTime);
    history.startFlowHistory.push_back(startFlows);
    history.endFlowHistory.push_back(endFlows);
    history.startHeadHistory.push_back(startHeads);
    history.endHeadHistory.push_back(endHeads);
    history.junctionHeadHistory.push_back(junctionHeads);
    
    if (!wavePoints.empty() || !history.waveAttenuationHistory.empty()) {
        history.waveAttenuationHistory.push_back(wavePoints);
    }
}