/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "demand.h"
#include "pattern.h"
#include <random>
#include <iostream>

//-----------------------------------------------------------------------------

//  Demand Constructor

Demand::Demand() :
    baseDemand(0.0),
    fullDemand(0.0),
    timePattern(nullptr)
{
}

//-----------------------------------------------------------------------------

//  Demand Destructor

Demand::~Demand() {}

//-----------------------------------------------------------------------------

//  Demand copy Constructor

Demand::Demand(const Demand& other)
{
	baseDemand = other.baseDemand;
	fullDemand = other.fullDemand;
	timePattern = other.timePattern;
}

Demand& Demand::operator=(const Demand& other)
{
    if (this != &other)
    {
        // Copy all properties
        
        baseDemand = other.baseDemand;
		fullDemand = other.fullDemand;
		timePattern = other.timePattern;
       
        // If you have any other members in the Demand class, copy them here
    }
    return *this;
}

//    Find pattern-adjusted full demand

double Demand::getFullDemand(double multiplier, double patternFactor) const
{
    if (timePattern) patternFactor = timePattern->currentFactor();

    double fullDemand = multiplier * baseDemand * patternFactor;

    return fullDemand;
}
