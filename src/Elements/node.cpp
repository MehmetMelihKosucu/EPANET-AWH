/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "node.h"
#include "junction.h"
#include "reservoir.h"
#include "tank.h"
#include "qualsource.h"
#include "Utilities/mempool.h"
#include "Core/network.h"
#include <iostream>
#include "link.h"
using namespace std;

//-----------------------------------------------------------------------------

// Constructor

Node::Node(const std::string& name_) :
    Element(name_),
    name(name_)             // Initialize the name
{
    // Initialize other members in the constructor body
    rptFlag = false;
    elev = 0.0;
    xCoord = -1e20;
    yCoord = -1e20;
    initQual = 0.0;
    qualSource = nullptr;
    fixedGrade = false;
    head = 0.0;
    h1ini = 0.0;
    h2ini = 0.0;
    pastHead = 0.0;
    ph = 0.0;
    qGrad = 0.0;
    fullDemand = 0.0;
    actualDemand = 0.0;
    outflow = 0.0;
    quality = 0.0;
    pastOutflow = 0.0;
}

// Destructor

Node::~Node()
{
    delete qualSource;
}

//-----------------------------------------------------------------------------

// Factory Method

Node* Node::factory(int type_, string name_, MemPool* memPool)
{
    Node* newNode = nullptr;

    switch (type_) {
    case JUNCTION:
        newNode = new(memPool->alloc(sizeof(Junction))) Junction(name_);
        break;
    case TANK:
        newNode = new(memPool->alloc(sizeof(Tank))) Tank(name_); 
        break;
    case RESERVOIR:
        newNode = new(memPool->alloc(sizeof(Reservoir))) Reservoir(name_); 
        break;
    default:
        return nullptr;
    }

    if (!newNode) {
        throw std::runtime_error("Node::factory() failed: Memory allocation returned nullptr");
    }

    return newNode;
}

//-----------------------------------------------------------------------------

void Node::initialize(Network* nw)
{
    head = elev;
    pastHead = elev;   // head on the previous timestep
    ph = elev;         // synonym of the past head
    h1ini = elev;
    h2ini = elev;
    quality = initQual;

    if (qualSource)
        qualSource->quality = quality;

    actualDemand = 0.0;
    outflow = 0.0;

    if (type() == JUNCTION)
        fixedGrade = false;
    else
        fixedGrade = true;
}



