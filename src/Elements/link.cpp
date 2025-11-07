/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "link.h"
#include "pipe.h"
#include "pump.h"
#include "valve.h"
#include "Core/constants.h"
#include "Core/network.h"
#include "Utilities/mempool.h"
#include "Core/datamanager.h"
#include "epanet3.h"
#include "node.h"
#include "junction.h"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <algorithm>

using namespace std;

//-----------------------------------------------------------------------------

static const string s_From = " status changed from ";
static const string s_To =   " to ";
static const string linkStatusWords[] = {"CLOSED", "OPEN", "ACTIVE", "TEMP_CLOSED"};

//-----------------------------------------------------------------------------

/// Constructor

Link::Link(std::string name_) :
    Element(name_),
    rptFlag(false),
    fromNode(nullptr),
    toNode(nullptr),
    initStatus(LINK_OPEN),
    diameter(0.0),
    lossCoeff(0.0),
    initSetting(1.0),
    status(0),
    previousStatus(0),
    flow(0.0),
    pastFlow(0.0),
    leakage(0.0),
    hLoss(0.0),
    pastHloss(0.0),
    hGrad(0.0),
    setting(0.0),
    quality(0.0),
    inertialTerm(0.0)
{
}

/// Destructor

Link::~Link() = default;

/// Factory Method

Link* Link::factory(int type, std::string name_, MemPool* memPool, Network* network)
{
    Link* link = nullptr;
    switch (type)
    {
    case Link::PIPE:
        link = new Pipe(name_);
        break;
    case Link::PUMP:
        link = new Pump(name_);
        break;
    case Link::VALVE:
        link = new Valve(name_);
        break;
    default:
        throw std::runtime_error("Unknown link type in Link::factory().");
    }
    // Assign pointers if needed
    link->memPool = memPool;
    link->network = network;
    return link;
}

//-----------------------------------------------------------------------------

void Link::initialize(bool reInitFlow)
{
    status = initStatus;
    setting = initSetting;
    if ( reInitFlow )
    {
        if ( status == LINK_CLOSED ) flow = ZERO_FLOW;
        else setInitFlow();
    }
    leakage = 0.0;
	
    inertialTerm = 0; // length / (GRAVITY * area);

}

//-----------------------------------------------------------------------------

double Link::getUnitHeadLoss()
{
    return hLoss;
}

//-----------------------------------------------------------------------------

string Link::writeStatusChange(int oldStatus)
{
    stringstream ss;
    ss << "    " << typeStr() << " " <<
	    name << s_From << linkStatusWords[oldStatus] << s_To <<
	    linkStatusWords[status];
    return ss.str();
}
















