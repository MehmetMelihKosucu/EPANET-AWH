/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

//! \file segpool.h
//! \brief Describes the SegPool class used for water quality transport.

#ifndef SEGPOOL_H_
#define SEGPOOL_H_

#include <vector>
#include "node.h"

class MemPool;

// Segment class to represent individual pipe segments
class Segment {
public:
    double flow;   // Flow through the segment
    double pastFlow; // Flow in the previous timestep
    double head;   // Head at the segment
    double pastHead; // Head in the previous timestep
    double length; // Length of the segment
    Node* nodeA;   // Upstream node
    Node* nodeB;   // Downstream node
    double  v;                //!< volume (ft3)
    double  c;                //!< constituent concentration (mass/ft3)
    Segment* next;    //!< next upstream volume segment

    Segment(double len, Node* upNode, Node* downNode);
};

class SegPool
{
  public:
    SegPool();
    ~SegPool();
    void init();
    Segment* getSegment(double v, double c);
    void     freeSegment(Segment* seg);

  private:
	int        segCount;     // number of volume segments allocated
	Segment*   freeSeg;      // first unused segment
	MemPool*   memPool;      // memory pool for volume segments
};

#endif // SEGPOOL_H_
