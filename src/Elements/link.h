/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

 //! \file link.h
 //! \brief Describes the Link class.

#ifndef LINK_H_
#define LINK_H_

#include "Elements/element.h"
#include "node.h"
#include "network.h"
#include <string>
#include <iostream>
#include <vector>
#include "segpool.h"
#include <memory>
#include "mempool.h"

// Forward declarations
class Node;
class Network;
class MemPool;

class Node;

class Link : public Element {
public:
    // Enums
    enum LinkType { PIPE, PUMP, VALVE };
    enum LinkStatus { LINK_CLOSED, LINK_OPEN, VALVE_ACTIVE, TEMP_CLOSED };
    enum LinkReaction { BULK, WALL };
    enum FlowModel { INCOMPRESSIBLE, COMPRESSIBLE };

    FlowModel flowModel;
    double wallThickness;
    double typeCoefficient;

    // Constructor and Destructor
    Link(std::string name_);
    virtual ~Link();// = default;

    // Factory Method
    static Link* factory(int type_, std::string name_, MemPool* memPool, Network* network);

    // Virtual methods (abstract for derived classes)
    virtual int type() const = 0;
    virtual std::string typeStr() = 0;
    virtual void convertUnits(Network* nw) = 0;
    virtual void findHeadLoss(Network* nw, double q) = 0;
    
    bool rptFlag;       //!< Indicates if results should be reported
    double initSetting; //!< Initial setting of the link (e.g., valve opening)
    virtual double      convertSetting(Network* nw, double s) { return s; }

    // Hydraulic initialization and computation
    virtual void initialize(bool initFlow);
    virtual void setInitFlow() {}
    virtual void setInitStatus(int s) {}
    virtual void setInitSetting(double s) {}
    virtual void setResistance(Network* nw) {}
    virtual void setLossFactor() {}

    // Getters for computed variables
    virtual double getVelocity() { return 0.0; }
    virtual double getRe(const double q, const double viscos) { return 0.0; }
    virtual double getResistance() { return 0.0; }
    virtual double getUnitHeadLoss();
    virtual double getSetting(Network* nw) { return setting; }
    
    double area;
    void setIndex(int idx) { index = idx; }
	double getArea() const { return area; }
	virtual double getVolume() { return 0.0; }
        
    // Energy usage and leakage
    virtual double updateEnergyUsage(Network* nw, int dt) { return 0.0; }
    virtual bool canLeak() { return false; }
    virtual double findLeakage(Network* nw, double h, double& dqdh) { return 0.0; }

    // Status and setting adjustment
    virtual void updateStatus(double q, double h1, double h2) {}
    virtual bool changeStatus(int newStatus, bool makeChange, const std::string reason, std::ostream& msgLog) { return false; }
    virtual bool changeSetting(double newSetting, bool makeChange, const std::string reason, std::ostream& msgLog) { return false; }
    virtual void validateStatus(Network* nw, double qTol) {}
    virtual void applyControlPattern(std::ostream& msgLog) {}
	virtual void validate(Network* nw) {}
    // Specialized link checks
    virtual bool isPRV() { return false; }
    virtual bool isPSV() { return false; }
    virtual bool isHpPump() { return false; }
	virtual bool isReactive() { return false; }

    // Miscellaneous
    std::string writeStatusChange(int oldStatus);

    // Properties
    Node* fromNode;            //!< Pointer to the link's start node
    Node* toNode;              //!< Pointer to the link's end node
    int initStatus;            //!< Initial Open/Closed status
    double diameter;           //!< Link diameter (ft)
    double lossCoeff;          //!< Minor head loss coefficient
    double setting;            //!< Current setting
    double quality;            //!< Avg. quality concentration (mass/ft�)
    
    int status;                //!< Current status
    int previousStatus;        //!< Status on the previous time step
    double flow;               //!< Flow rate (cfs)
    double pastFlow;           //!< Flow rate in the previous timestep
    double leakage;            //!< Leakage rate (cfs)
    double hLoss;              //!< Head loss (ft)
    double pastHloss;          //!< Head loss in the previous timestep
    double hGrad;              //!< Head loss gradient (ft/cfs)
    double pastSetting;        //!< Setting on the previous timestep
    double inertialTerm;       //!< Term for inertia calculations
    double phiA;               //!< Absolute flow indicator
    double phiR;               //!< Relative flow indicator
    Network* network;          // Pointer to the network
	
private:
    
protected:
    
    MemPool* memPool;          // Pointer to the memory pool
};


#endif
