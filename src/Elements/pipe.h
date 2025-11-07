/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file pipe.h
//! \brief Describes the Pipe class.

#ifndef PIPE_H_
#define PIPE_H_

#include "Elements/link.h"
#include "Elements/junction.h"

class Network;
class Node;
class Junction;


//! \class Pipe
//! \brief A circular conduit Link through which water flows.

class Pipe: public Link
{
  public:

    // Constructor/Destructor

    Pipe(std::string name);

    ~Pipe();
    
    // Methods
    int         type() const override { return Link::PIPE; }
    std::string typeStr() override { return "Pipe"; }
    void        convertUnits(Network* nw);
    bool        isReactive();
    void        setInitFlow();
    void        setInitStatus(int s);
    void        setInitSetting(double s);
    void        setResistance(Network* nw);
    void		setLossFactor();
    double      getRe(const double q, const double viscos);
    double      getResistance() {return resistance;}
    double      getVelocity();
    double      getUnitHeadLoss();
    double      getSetting(Network* nw) { return roughness; }
    double      getVolume() { return 0.785398 * length * diameter * diameter; }
    
    double getArea() const 
    {
        // π × diameter²/4
        if (diameter <= 0.0) {
            // Return a small but non-zero area for invalid diameters
            return 0.001;  // Some minimum value to prevent division by zero
        }
        return 0.7853981634 * diameter * diameter; // π/4 = 0.7853981634
    }

    double      getWaveSpeed() const;
    void        findHeadLoss(Network* nw, double q);
   
    bool        canLeak() { return leakCoeff1 > 0.0; }
    double      findLeakage(Network* nw, double h, double& dqdh);
    bool        changeStatus(int s, bool makeChange,
                            const std::string reason,
                            std::ostream& msgLog);
    void        validateStatus(Network* nw, double qTol);
    double getSegmentLength() const { return length / numSegments; }

    // Properties

    bool   hasCheckValve;    //!< true if pipe has a check valve
    double length;           //!< pipe length (ft)
    double roughness;        //!< roughness parameter (units depend on head loss model)
    double resistance;       //!< resistance factor (units depend head loss model)
    double lossFactor;       //!< minor loss factor (ft/cfs^2)
    double leakCoeff1;       //!< leakage coefficient (user units)
    double leakCoeff2;       //!< leakage coefficient (user units)
    double bulkCoeff;        //!< bulk reaction coefficient (mass^n/sec)
    double wallCoeff;        //!< wall reaction coefficient (mass^n/sec)
    double massTransCoeff;   //!< mass transfer coefficient (mass^n/sec)
	double diameter;         //!< pipe diameter (ft)
	double waveSpeed;       //!< Wave propagation speed
    //double flow;             //!< flow (cfs)
    double startFlow;      // Flow at upstream end (NodeA)
    double pastStartFlow;  // Previous startFlow
    double endFlow;        // Flow at downstream end (NodeB)
    double pastEndFlow;    // Previous endFlow
    double dStartFlow;     // Flow change at upstream end (NodeA)
    double dEndFlow;       // Flow change at downstream end (NodeB)
    // Interior flow values for GCM
    double interiorQa;      // Current flow at first interior point
    double interiorQb;      // Current flow at last interior point
    double pastInteriorQa;  // Past flow at first interior point
    double pastInteriorQb;  // Past flow at last interior point
    // Interior head values for GCM
    double interiorHa;      // Head at first interior point
    double interiorHb;      // Head at last interior point
    
	//int segments;         //!< number of segments used to model pipe
    
    double getFrictionFactor(double Re) const;

    int numSegments;                       // Number of segments.

    double courantNumber; 

    double computeShearDecayCoeff(double flow) const;
    double getLength() const;

    double Pipe::calcSwameeJain(double Re) const;
    // Characteristic method (MOC) properties
    double charC_plus;     // Positive characteristic value (C+)
    double charC_minus;    // Negative characteristic value (C-)
    double charB_plus;     // Positive characteristic coefficient (B+)
    double charB_minus;    // Negative characteristic coefficient (B-)
    int _numDiscretizedReaches;  // Number of theoretical reaches for Courant condition

    int getNumDiscretizedReaches() const;

    int numReaches; // Number of reaches
    int reachBackSteps; // Number of reach-back steps
    int impingingWaveCount = 0; // Number of impinging waves
    double waveAttenuation; // Wave attenuation factor
    std::vector<std::vector<double>>    wavePathStorage;
    //std::vector<WaveAttenuationPoint> wavePathStorage;
    double characteristicB;
    double resistanceCoeff;

    double startFlowCoeff;     // Coefficient of head term for startFlow calculation
    double startFlowConstant;  // Constant term for startFlow calculation  
    double endFlowCoeff;       // Coefficient of head term for endFlow calculation
    double endFlowConstant;    // Constant term for endFlow calculation

  private:
    
    
 };


#endif
