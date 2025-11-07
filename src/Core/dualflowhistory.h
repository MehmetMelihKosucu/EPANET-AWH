/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION */


// dualflowhistory.h
#ifndef DUALFLOWHISTORY_H
#define DUALFLOWHISTORY_H

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Core/network.h"
#include <set>

// Forward declarations
class Pipe;
class Network;
class Junction;

/**
 * @struct DualFlowHistory
 * @brief Stores historical flow and head data for water hammer analysis
 * 
 * This structure maintains time-series data for pipe flows and heads at both
 * ends of a pipe, allowing for reach-back operations needed in algebraic
 * water hammer (AWH) calculations.
 */


/**
 * @struct FlowHistoryResult
 * @brief Simple structure to return reach-back flow and head values
 * 
 * This structure encapsulates the results of a reach-back operation,
 * providing a clean interface for retrieving historical values.
 */
struct FlowHistoryResult {
    double startFlow = 0.0;      // Q_a
    double endFlow = 0.0;        // Q_b
    double startHead = 0.0;      // H_a
    double endHead = 0.0;        // H_b
    bool found = false;          // Whether the history was found
};

struct WaveAttenuationPoint {
    double position;    // Position along pipe (0.0 to 1.0)
    double flow;       // Flow at this position
    double head;       // Head at this position
    double amplitude;  // Wave amplitude
    
    // Default constructor
    WaveAttenuationPoint() : position(0.0), flow(0.0), head(0.0), amplitude(0.0) {}
    
    // Parameterized constructor
    WaveAttenuationPoint(double pos, double f, double h, double amp = 0.0) 
        : position(pos), flow(f), head(h), amplitude(amp) {}
};


/**
 * @class FlowHistoryManager
 * @brief Manages a collection of flow histories for multiple pipes
 * 
 * This class handles the creation, storage, and retrieval of flow history
 * data for all pipes in a hydraulic network. It provides an interface for
 * both initializing and accessing historical data needed for water hammer
 * calculations.
 */

class DualFlowHistory {
private:
    // Fixed-size storage vectors (pre-allocated)
    std::vector<double> timestamps;
    std::vector<std::vector<double>> startFlowHistory;
    std::vector<std::vector<double>> endFlowHistory;
    std::vector<std::vector<double>> startHeadHistory;
    std::vector<std::vector<double>> endHeadHistory;
    std::vector<std::vector<double>> junctionHeadHistory;
    std::vector<std::vector<double>> pointAFlowHistory;
    std::vector<std::vector<double>> pointBFlowHistory;
    std::vector<std::vector<double>> pointAHeadHistory;
    std::vector<std::vector<double>> pointBHeadHistory;
    std::vector<std::vector<WaveAttenuationPoint>> waveAttenuationHistory;
    
    // Circular buffer management
    size_t capacity;      // Maximum number of entries
    size_t currentSize;   // Current number of valid entries
    size_t headIndex;     // Index where next entry will be written
    size_t tailIndex;     // Index of oldest entry
    bool isCircular;      // Whether we're in circular mode (buffer is full)
    
    /**
     * Convert logical index (0 = oldest) to actual storage index
     * This is the key to making the circular buffer work transparently
     */
    size_t getStorageIndex(size_t logicalIndex) const {
        if (!isCircular) {
            // Still in linear mode - direct mapping
            return logicalIndex;
        }
        // In circular mode - map through tail position
        return (tailIndex + logicalIndex) % capacity;
    }
    
    /**
     * Safe version with bounds checking for debugging
     */
    size_t getStorageIndexSafe(size_t logicalIndex, const std::string& caller) const {
        // First, check if the logical index is valid
        if (logicalIndex >= currentSize) {
            std::cerr << "\n[ERROR] DualFlowHistory bounds check failed!" << std::endl;
            std::cerr << "  Caller: " << caller << std::endl;
            std::cerr << "  Requested index: " << logicalIndex << std::endl;
            std::cerr << "  Current size: " << currentSize << std::endl;
            throw std::out_of_range("Index out of range in " + caller);
        }
        return getStorageIndex(logicalIndex);
    }
    
public:
    /**
     * Constructor with specified capacity
     * @param maxCapacity Maximum number of time points to store (default 1000)
     */
    explicit DualFlowHistory(size_t maxCapacity = 1000) : 
        capacity(maxCapacity), 
        currentSize(0), 
        headIndex(0), 
        tailIndex(0),
        isCircular(false) {
        
        // Pre-allocate all storage to avoid reallocations
        timestamps.resize(capacity);
        startFlowHistory.resize(capacity);
        endFlowHistory.resize(capacity);
        startHeadHistory.resize(capacity);
        endHeadHistory.resize(capacity);
        junctionHeadHistory.resize(capacity);
        pointAFlowHistory.resize(capacity);
        pointBFlowHistory.resize(capacity);
        pointAHeadHistory.resize(capacity);
        pointBHeadHistory.resize(capacity);
        waveAttenuationHistory.resize(capacity);
    }
    
    // Copy constructor for compatibility
    DualFlowHistory(const DualFlowHistory& other) = default;
    
    // Assignment operator
    DualFlowHistory& operator=(const DualFlowHistory& other) = default;
    
    /**
     * Clear all history data and reset indices
     */
    void clear();
    
    /**
     * Check if we have enough history for the requested reach-back
     */
    bool hasEnoughHistory(int reachBackSteps) const;
    
    /**
     * Get the current number of valid entries
     */
    size_t size() const { return currentSize; }
    
    /**
     * Check if the buffer is empty
     */
    bool empty() const { return currentSize == 0; }
    
    /**
     * Get the buffer capacity
     */
    size_t getCapacity() const { return capacity; }
    
    /**
     * Add a new entry to the history
     * This automatically overwrites the oldest entry if the buffer is full
     */
    void addEntry(double timestamp,
                  const std::vector<double>& startFlows,
                  const std::vector<double>& endFlows,
                  const std::vector<double>& startHeads,
                  const std::vector<double>& endHeads,
                  const std::vector<double>& junctionHeads,
                  const std::vector<WaveAttenuationPoint>& wavePoints = {}) {
        
        // Store data at the head position
        timestamps[headIndex] = timestamp;
        startFlowHistory[headIndex] = startFlows;
        endFlowHistory[headIndex] = endFlows;
        startHeadHistory[headIndex] = startHeads;
        endHeadHistory[headIndex] = endHeads;
        junctionHeadHistory[headIndex] = junctionHeads;
        
        // Handle optional wave points
        if (!wavePoints.empty()) {
            waveAttenuationHistory[headIndex] = wavePoints;
        } else {
            waveAttenuationHistory[headIndex].clear();
        }
        
        // Update indices
        if (currentSize < capacity) {
            // Still filling up the buffer
            currentSize++;
            headIndex++;
            if (headIndex >= capacity) {
                headIndex = 0;  // Wrap around
                isCircular = true;
            }
        } else {
            // Buffer is full - we're in circular mode
            isCircular = true;
            headIndex = (headIndex + 1) % capacity;
            tailIndex = (tailIndex + 1) % capacity;
        }
    }
    
    /**
     * Get timestamp at logical index (0 = oldest)
     * This is what replaces direct access to timestamps[idx]
     */
    double getTimestamp(size_t logicalIndex) const {
        if (logicalIndex >= currentSize) {
            throw std::out_of_range("Timestamp index out of range");
        }
        size_t idx = getStorageIndex(logicalIndex);
        return timestamps[idx];
    }
    
    /**
     * Get flow/head values at logical index
     * These replace direct access to the history arrays
     */
    const std::vector<double>& getStartFlow(size_t logicalIndex) const {
        if (logicalIndex >= currentSize) {
            throw std::out_of_range("Flow history index out of range");
        }
        size_t idx = getStorageIndex(logicalIndex);
        return startFlowHistory[idx];
    }
    
    const std::vector<double>& getEndFlow(size_t logicalIndex) const {
        if (logicalIndex >= currentSize) {
            throw std::out_of_range("Flow history index out of range");
        }
        size_t idx = getStorageIndex(logicalIndex);
        return endFlowHistory[idx];
    }
    
    const std::vector<double>& getStartHead(size_t logicalIndex) const {
        if (logicalIndex >= currentSize) {
            throw std::out_of_range("Head history index out of range");
        }
        size_t idx = getStorageIndex(logicalIndex);
        return startHeadHistory[idx];
    }
    
    const std::vector<double>& getEndHead(size_t logicalIndex) const {
        if (logicalIndex >= currentSize) {
            throw std::out_of_range("Head history index out of range");
        }
        size_t idx = getStorageIndex(logicalIndex);
        return endHeadHistory[idx];
    }
    
    /**
     * Find the best indices for interpolating at a target time
     * Returns a pair of (lowerIndex, upperIndex) in logical space
     */
    std::pair<size_t, size_t> findInterpolationIndices(double targetTime) const {
        if (currentSize == 0) {
            throw std::runtime_error("Cannot interpolate in empty history");
        }
        
        // Handle edge cases
        if (currentSize == 1 || targetTime <= getTimestamp(0)) {
            return {0, 0};
        }
        if (targetTime >= getTimestamp(currentSize - 1)) {
            return {currentSize - 1, currentSize - 1};
        }
        
        // Binary search in logical index space
        size_t left = 0;
        size_t right = currentSize - 1;
        
        while (left < right) {
            size_t mid = left + (right - left) / 2;
            if (getTimestamp(mid) < targetTime) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        // Return interpolation indices
        if (left == 0) return {0, 0};
        return {left - 1, left};
    }
    
    /**
     * Update the most recent entry if timestamps match
     * This is useful when the solver re-evaluates the same time point
     */
    bool updateLastEntry(double timestamp,
                        const std::vector<double>& startFlows,
                        const std::vector<double>& endFlows,
                        const std::vector<double>& startHeads,
                        const std::vector<double>& endHeads,
                        const std::vector<double>& junctionHeads,
                        const std::vector<WaveAttenuationPoint>& wavePoints = {}) {
        
        if (currentSize == 0) return false;
        
        // Get the index of the most recent entry
        size_t lastIndex = (headIndex == 0) ? capacity - 1 : headIndex - 1;
        
        // Check if timestamps match (within tolerance)
        if (std::abs(timestamps[lastIndex] - timestamp) < 1e-9) {
            // Update the existing entry
            startFlowHistory[lastIndex] = startFlows;
            endFlowHistory[lastIndex] = endFlows;
            startHeadHistory[lastIndex] = startHeads;
            endHeadHistory[lastIndex] = endHeads;
            junctionHeadHistory[lastIndex] = junctionHeads;
            
            if (!wavePoints.empty()) {
                waveAttenuationHistory[lastIndex] = wavePoints;
            }
            
            return true;
        }
        
        return false;
    }
    
    // Friend class for any special access needs
    friend class FlowHistoryManager;
};



class FlowHistoryManager {
public:
    static FlowHistoryManager& getInstance() {
        static FlowHistoryManager instance;
        return instance;
    }
    
    void FlowHistoryManager::validateReachBackValues(
    Pipe* pipe,
    FlowHistoryResult& result,
    double currentTime,
    Network* network);

    // Delete copy constructor and assignment operator
    FlowHistoryManager(const FlowHistoryManager&) = delete;
    FlowHistoryManager& operator=(const FlowHistoryManager&) = delete;
    
    // Core functionality
    void initializeHistory(Network* network, double currentTime, double timeStep, double historyStartTime);
    void addHistory(Pipe* pipe, const DualFlowHistory& history);
    FlowHistoryResult FlowHistoryManager::getReachBackValues(
    Pipe* pipe, 
    double currentTime, 
    double waveTravel, 
    Network* network)   ;
    
    void FlowHistoryManager::updateHistoryForPipe(
    Pipe* pipe, 
    double currentTime,
    double actualTimeStep,
    const std::vector<double>& startFlows,
    const std::vector<double>& endFlows,
    const std::vector<double>& startHeads,
    const std::vector<double>& endHeads,
    const std::vector<double>& junctionHeads,
    const std::vector<WaveAttenuationPoint>& wavePoints);
    
    // Wave calculations
    void calculateWaveAttenuationPoints(
        Pipe* pipe,
        std::vector<WaveAttenuationPoint>& wavePoints,
        double startFlow,
        double endFlow,
        double startHead,
        double endHead,
        const DualFlowHistory& history,
        Network* network);
    
    void applyWaveTracking(
        Pipe* pipe,
        WaveAttenuationPoint& point,
        int reachBackIndex,
        const DualFlowHistory& history,
        Network* network);
    
    // Utility functions
    void clear();
    bool hasHistoryForPipe(Pipe* pipe) const;
    double linearInterpolate(double y1, double y2, double alpha);

    const std::map<Pipe*, std::unique_ptr<DualFlowHistory>>& getHistoryData() const {
        return historyData;
    }
    
private:
    FlowHistoryManager() {}
    std::map<Pipe*, std::unique_ptr<DualFlowHistory>> historyData;

    // Singleton instance
    static FlowHistoryManager* instance;

    // Map to track which pipes still need history building from incompressible phase
    std::map<Pipe*, bool> pipeHistoryStatus;
    
    // Track when we transitioned to compressible flow
    double compressibleFlowStartTime;
    bool hasTransitionedToCompressible;

    // Track pump shutdown events
    struct PumpShutdownEvent {
        double time;
        Pump* pump;
        Node* dischargeNode;
        double flowBefore;
        std::set<Pipe*> affectedPipes;
    };
    
    std::vector<PumpShutdownEvent> pumpShutdowns;
    std::map<Pipe*, std::vector<PumpShutdownEvent*>> pipeShutdownEvents;
};
#endif // DUALFLOWHISTORY_H