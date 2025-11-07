/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "segpool.h"
#include "Utilities/mempool.h"
#include "node.h"

#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------

SegPool::SegPool()
{
    memPool = new MemPool();
    freeSeg = nullptr;
    segCount = 0;
}

//-----------------------------------------------------------------------------

SegPool::~SegPool()
{
//    cout << "\n segCount = " << segCount << "\n";

    delete memPool;
}

//-----------------------------------------------------------------------------

void SegPool::init()
{
    segCount = 0;
    memPool->reset();
    freeSeg = nullptr;
}

//-----------------------------------------------------------------------------

Segment* SegPool::getSegment(double v, double c)
{
    // ... if there's a free segment available then use it
    Segment* seg;
    if ( freeSeg )
    {
       seg = freeSeg;
       freeSeg = seg->next;
    }

    // ... otherwise create a new one from the memory pool
    else
    {
        seg = (Segment *) memPool->alloc(sizeof(Segment));
        segCount++;
    }

    // ... assign segment's volume and quality
    if ( seg )
    {
        seg->v = v;
        seg->c = c;
        seg->next = nullptr;
    }
    return seg;
}

//-----------------------------------------------------------------------------

void SegPool::freeSegment(Segment* seg)
{
    seg->flow = 0.0;
    seg->pastFlow = 0.0;
    seg->head = 0.0;
    seg->pastHead = 0.0;
    seg->v = 0.0;
    seg->c = 0.0;
    seg->nodeA = nullptr;
    seg->nodeB = nullptr;
    seg->next = freeSeg;
    freeSeg = seg;
}


Segment::Segment(double len, Node* upNode, Node* downNode)
    : length(len), nodeA(upNode), nodeB(downNode), flow(0.0), pastFlow(0.0),
    head(0.0), pastHead(0.0), v(0.0), c(0.0), next(nullptr) {
}