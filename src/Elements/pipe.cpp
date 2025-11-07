/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "pipe.h"
#include "Core/network.h"
#include "Core/constants.h"
#include "Models/headlossmodel.h"
#include "Models/leakagemodel.h"
#include "junction.h"

#include <iostream>
#include <cmath>
using namespace std;

//-----------------------------------------------------------------------------

Pipe::Pipe(std::string name)
      : Link(name), hasCheckValve(false), length(0), roughness(0),
        resistance(0), lossFactor(0), leakCoeff1(0), leakCoeff2(0),
        bulkCoeff(0), wallCoeff(0), massTransCoeff(0), waveSpeed(0)
{ 
    // Initialize area based on diameter
    diameter = 0.0; // Default initialization
    area = 0.0;     // Default initialization
} //*/


Pipe::~Pipe() {}

void Pipe::convertUnits(Network* nw)
{
    diameter /= nw->ucf(Units::DIAMETER);
    length   /= nw->ucf(Units::LENGTH);

    // ... convert minor loss coeff. from V^2/2g basis to Q^2 basis

    lossFactor = 0.02517 * lossCoeff / pow(diameter, 4);

    // ... convert roughness length units of Darcy-Weisbach headloss model
    //     (millifeet or millimeters to feet)

    if ( nw->option(Options::HEADLOSS_MODEL ) == "D-W")
    {
        roughness = roughness / nw->ucf(Units::LENGTH) / 1000.0;
    }

    // ... apply global default leakage coeffs.

    if ( leakCoeff1 == MISSING ) leakCoeff1 = nw->option(Options::LEAKAGE_COEFF1);
    if ( leakCoeff2 == MISSING ) leakCoeff2 = nw->option(Options::LEAKAGE_COEFF2);

    // ... apply global default reaction coeffs.

    if ( bulkCoeff == MISSING ) bulkCoeff = nw->option(Options::BULK_COEFF);
    if ( wallCoeff == MISSING ) wallCoeff = nw->option(Options::WALL_COEFF);
}

//-----------------------------------------------------------------------------

bool Pipe::isReactive()
{
    if ( bulkCoeff != 0.0 ) return true;
    if ( wallCoeff != 0.0 ) return true;
    return false;
}

//-----------------------------------------------------------------------------

void Pipe::setResistance(Network* nw)
{
    nw->headLossModel->setResistance(this);
}

//-----------------------------------------------------------------------------

void Pipe::setInitStatus(int s)
{
    initStatus = s;
}

//-----------------------------------------------------------------------------

void Pipe::setInitSetting(double s)
{
    if ( s == 0.0 ) initStatus = LINK_CLOSED;
    else            initStatus = LINK_OPEN;
}

//-----------------------------------------------------------------------------

void Pipe::setInitFlow()
{
    // ... flow at velocity of 1 ft/s
    double initialFlow = PI * diameter * diameter / 4.0;
    startFlow = PI * diameter * diameter / 4.0;
    endFlow = PI * diameter * diameter / 4.0;
    pastStartFlow = PI * diameter * diameter / 4.0;
    pastEndFlow = PI * diameter * diameter / 4.0;
    // Set overall flow as average of start and end flows
    flow = PI * diameter * diameter / 4.0;
    pastFlow = 0;
    pastHloss = 0;
}

//-----------------------------------------------------------------------------

void Pipe::setLossFactor()
{
    lossFactor = 0.02517 * lossCoeff / pow(diameter, 4);
}

double Pipe::getVelocity()
{
    // Use average of startFlow and endFlow for velocity calculation
    double avgFlow = (startFlow + endFlow) / 2.0;
    double area = PI * diameter * diameter / 4.0;
    return abs(flow) / area;
}

//-----------------------------------------------------------------------------

double Pipe::getUnitHeadLoss()
{
    if ( length > 0.0 ) return abs(hLoss) * 1000.0 / length;
    return 0.0;
}

//-----------------------------------------------------------------------------

double Pipe::getRe(const double q, const double viscos)
{
    return abs(q) / (PI * diameter * diameter / 4.0) * diameter / viscos;
}

//-----------------------------------------------------------------------------

void Pipe::findHeadLoss(Network* nw, double q) {
    // Use average of start and end flows for head loss calculation
    //double avgFlow = (startFlow + endFlow) / 2.0;
    
    // Handle closed links
    if (status == LINK_CLOSED || status == TEMP_CLOSED) {
        HeadLossModel::findClosedHeadLoss(q, hLoss, hGrad);
        inertialTerm = MIN_GRADIENT;
        return;
    }
    else {
        // Standard head loss calculation (no segmentation)
        nw->headLossModel->findHeadLoss(this, q, hLoss, hGrad);
        if (hasCheckValve) HeadLossModel::addCVHeadLoss(q, hLoss, hGrad);
        inertialTerm = (length * 4) / (32.174 * PI * diameter * diameter);
    }
}


//-----------------------------------------------------------------------------

double Pipe::findLeakage(Network* nw, double h, double& dqdh)
{
    return nw->leakageModel->findFlow(leakCoeff1, leakCoeff2, length, h, dqdh);
}


//-----------------------------------------------------------------------------

bool Pipe::changeStatus(int s, bool makeChange, const string reason, ostream& msgLog)
{
    if ( status != s )
    {
        if ( makeChange )
        {
            msgLog << "\n    " << reason;
            status = s;
        }
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------

// For debugging only

void Pipe::validateStatus(Network* nw, double qTol)
{
    if ( hasCheckValve && flow < -qTol )
    {
        nw->msgLog << "\nCV " << name << " flow = " << flow*nw->ucf(Units::FLOW);
    }
}

// Compute the unsteady friction coefficient for the link
double Pipe::computeShearDecayCoeff(double flow) const {
    const double D = diameter;    // Pipe diameter (m)
    const double A = PI * D * D / 4;       // Cross-sectional area (m?)
    const double nu = 1e-6;      // Kinematic viscosity of water (m?/s)

    // Compute Reynolds number
    double Re = fabs(flow) * D / (A * nu);

    // Compute the shear decay coefficient (k) based on flow regime
    if (Re < 2300) {
        // Laminar flow
        return 0.5 * sqrt(4.76e-3);
    }
    else {
        // Turbulent flow
        return 0.5 * sqrt(7.41 * pow(Re, -log10(14.3 * pow(Re, -0.05))));
    }
}

double Pipe::getWaveSpeed() const {
    // Constants in FPS (Imperial) units
    double E = 3.05e8;    // Young's modulus of ductile iron (psi)
    double rho = 1.94;    // Water density (slugs/ft^3)
    double K = 3.19e8;    // Bulk modulus of water (psi)
    double D = diameter;  // Outer diameter in feet

    // Wall thickness for 6-inch ductile iron pipe (Pressure Class 150)
    double e = 0.007 * diameter; // Wall thickness in feet

    // Safety check to prevent division by zero
    if (e <= 0.0) {
        throw std::invalid_argument("Wall thickness must be greater than zero.");
    }

    // Type coefficient (C), accounting for Poisson's ratio (? = 0.25)
    double typeCoefficient = 1;

    double waveSpeed = 0;

    /*if( index == 4 ) {
        typeCoefficient = 12;
    }
    else if (index == 7 || index == 8) {
        typeCoefficient = 1;
    }
    else{
        typeCoefficient = 1;
    } // */

    /*if( index == 4 ) {
        waveSpeed = 1000 / 0.3048;
    }
    else if (index == 7 || index == 8) {
        waveSpeed = 1000 / 0.3048;
    }
    else{
        waveSpeed = 1000 / 0.3048;
    } // */

    if( diameter <= 0.4 / 0.3048 ) {
        waveSpeed = 350 / 0.3048;
    }
    else{
        waveSpeed = 1050 / 0.3048;
    } // */

    // Compute wave speed (a) during water hammer
    //double waveSpeed = sqrt((K / rho) / (1 + (K / E) * (D / e) * typeCoefficient));

    return waveSpeed;  // Returns wave speed in ft/s
}


double Pipe::getFrictionFactor(double Re) const {
    // For zero or extremely low Reynolds numbers, avoid division by zero
    if (Re < 1.0) {
        return 64.0;  // Limiting case as Re approaches zero
    }

    // For laminar flow (Re < 2000)
    if (Re < 2000.0) {
        return 64.0 / Re;
    }

    // For transitional flow (2000 = Re < 4000)
    else if (Re < 4000.0) {
        // Smooth transition between laminar and turbulent regimes
        double fLaminar = 64.0 / 2000.0;
        double fTurbulent = calcSwameeJain(4000.0);
        double t = (Re - 2000.0) / 2000.0;  // Transition factor
        return fLaminar * (1.0 - t) + fTurbulent * t;
    }

    // For turbulent flow (Re = 4000)
    else {
        return calcSwameeJain(Re);
    }
}

// Helper method using the Swamee-Jain approximation formula
double Pipe::calcSwameeJain(double Re) const {
    // Relative roughness (e/D)
    double relativeRoughness = roughness / diameter;

    // Swamee-Jain approximation formula
    return 0.25 / pow(log10(relativeRoughness / 3.7 + 5.74 / pow(Re, 0.9)), 2);
}

// Getter for pipe length
double Pipe::getLength() const {
    return length;
}






