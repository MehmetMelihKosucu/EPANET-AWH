/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

//////////////////////////////////////////////
// Implementation of the HydBalance class.  //
//////////////////////////////////////////////

// TO DO:
// - compute and report system wide cumulative flow balance

#include "hydbalance.h"
#include "network.h"
#include "Elements/node.h"
#include "Elements/junction.h"
#include "Elements/link.h"
#include "Elements/valve.h"
#include "Elements/pipe.h"
#include "Elements/tank.h"
#include "Core/hydengine.h"
#include "Solvers/rwcggasolver.h"
#include "Models/headlossmodel.h"
#include "constants.h"
#include <utility>
#include <cmath>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include "dualflowhistory.h"
#include <queue>
//#include "Solvers/cggasolver.h" // Include the CGGASolver header file
using std::vector;
using std::map;
using std::set;
using std::less;
using std::allocator;
using std::pair;
using std::make_pair;

void   findNodeOutflows(double lamda, double dH[], double xQ[], Network* nw);
void   findLeakageFlows(double lamda, double dH[], double xQ[], Network* nw);
double findTotalFlowChange(double lamda, double dQ[], Network* nw);

//-----------------------------------------------------------------------------

//  Evaluate the error in satisfying the conservation of flow and energy
//  equations by an updated set of network heads and flows.
//

double HydBalance::evaluate(
    double lamda,      // step size
    double dH[],       // change in nodal heads
    double dQ[],       // change in link flows
    double xQ[],       // nodal inflow minus outflow
    Network* nw,       // network being analyzed
    int currentTime, 
    double tstep,
    int currentMode)   
{
    // Initialize error tracking metrics
    maxFlowErr = 0.0;
    maxHeadErr = 0.0;
    maxFlowChange = 0.0;
    maxHeadErrLink = -1;
    maxFlowErrNode = -1;
    maxFlowChangeLink = -1;
    
    // Calculate the total number of nodes including all interior junctions
    int regularNodeCount = nw->nodes.size();
    int interiorJunctionCount = 0;
    int totalSegmentCount = 0;
    
    for (int i = 0; i < regularNodeCount; i++) {
            xQ[i] = 0.0;
        }
    
    // Different error evaluation based on flow regime
    double norm = 0.0;
    
    // Find head loss errors and update xQ with internal link flows
    norm = findHeadErrorNorm(lamda, dH, dQ, xQ, nw, currentTime, tstep);
        
        // Update xQ with external outflows (demands, emitters, etc.)
    findNodeOutflows(lamda, dH, xQ, nw);
        
        // Add the flow balance error norm
    norm += findFlowErrorNorm(xQ, nw);
    
    // Calculate total relative flow change across all links
    totalFlowChange = findTotalFlowChange(lamda, dQ, nw);
    
    // Return the root mean square error
    return sqrt(norm);
}

//-----------------------------------------------------------------------------

//  Find the error norm in satisfying the head loss equation across each link.


double HydBalance::findHeadErrorNorm(
        double lamda, double dH[], double dQ[], double xQ[], Network* nw, int currentTime, double tstep)
{
    double norm = 0.0;
    double count = 0.0;
    maxHeadErr = 0.0;
    maxFlowChange = 0.0;
    maxFlowChangeLink = 0;

	

    int linkCount = nw->links.size();
    for (int i = 0; i < linkCount; i++)
    {

        // ... identify link's end nodes

        Link* link = nw->links[i];
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;

        //if (!link) continue;
        //Pipe* pipe = dynamic_cast<Pipe*>(link);

        // ... apply updated flow to end node flow balances

        double flowChange = lamda * dQ[i];
        double flow = link->flow + flowChange;
        xQ[n1] -= flow;
        xQ[n2] += flow;

        // ... update network's max. flow change

		previousMaxFlowChange = maxFlowChange;

        double err = abs(flowChange);
        if ( err > maxFlowChange )
        {
            maxFlowChange = err;
            maxFlowChangeLink = i;
        }

        int currentMode = 0; 

        // ... compute head loss and its gradient (head loss is saved
        // ... to link->hLoss and its gradient to link->hGrad)
//*******************************************************************
		//if ((currentMode == 0 || currentMode == 1) || currentTime == 0 || link->type() != Link::PIPE || nw->option(Options::HYD_SOLVER) != "CGGA")
        link->findHeadLoss(nw, flow);
		//else
            //pipe->findSegmentHeadLoss(nw, flow);
//*******************************************************************
	    
		double unsteadyTerm = 0;

        // ... evaluate head loss error according to Steady and Unsteady Flow Conditions
		

		if (currentTime == 0 || nw->option(Options::HYD_SOLVER) == "GGA" || currentMode == 0 )
		{
			unsteadyTerm = 0;
		}

		else if (currentMode == 1)
		{
			unsteadyTerm = (link->inertialTerm) * (link->flow - link->pastFlow) / tstep;
		}
		else if (currentMode == 2)
		{
			unsteadyTerm = (link->inertialTerm) * (link->flow - link->pastFlow) / tstep;
		}
        h1 = link->fromNode->head + lamda * dH[n1];
        h2 = link->toNode->head + lamda * dH[n2];
        if ( link->hGrad == 0.0 ) link->hLoss = h1 - h2;
        //err = h1 - h2 - link->hLoss;


		err = unsteadyTerm - h1 + h2 + link->hLoss;
       
		if ( abs(err) > maxHeadErr )
        {
            maxHeadErr = abs(err);
            maxHeadErrLink = i;
        }

        // ... update sum of squared errors

        norm += err * err;
        count += 1.0;
    }

    // ... return sum of squared errors normalized by link count

    if ( count == 0.0 ) return 0;
    else return norm / count;

}

//-----------------------------------------------------------------------------

//  Find net external outflow at each network node.

void findNodeOutflows(double lamda, double dH[], double xQ[], Network* nw)
{
    // ... initialize node outflows and their gradients w.r.t. head

    // First validate network pointer
    if (!nw) {
        throw std::runtime_error("Network pointer is null in findNodeOutflows");
    }

    // Initialize outflows and gradients with proper null checks
    for (Node* node : nw->nodes) {
        if (!node) {
            std::cerr << "Warning: Encountered null node during initialization" << std::endl;
            continue;
        }
        node->outflow = 0.0;
        node->qGrad = 0.0;
    }

    // ... find pipe leakage flows & assign them to node outflows

    if ( nw->leakageModel ) findLeakageFlows(lamda, dH, xQ, nw);

    // ... add emitter flows and demands to node outflows

    int nodeCount = nw->nodes.size();
    for (int i = 0; i < nodeCount; i++)
    {
        Node* node = nw->node(i);

        // First, perform thorough null checking
        if (!node) {
            std::cerr << "Error: Null node encountered in findNodeOutflows" << std::endl;
            continue;
        }

        // Now perform the type checking
        if (Node::NodeType() != Node::JUNCTION && Node::NodeType() != Node::TANK && Node::NodeType() != Node::RESERVOIR) {
            std::cerr << "Error: Invalid node type at index " << i << std::endl;
            continue;
        }

        double h = node->head + lamda * dH[i];
        double q = 0.0;
        double dqdh = 0.0;

        // ... for junctions, outflow depends on head

        if (node->type() == Node::JUNCTION) {
            // First check if this is an interior node (segment point)
            bool isInteriorNode = false;
            for (Link* link : nw->links) {
                if (link->type() == Link::PIPE) {
                    Pipe* pipe = dynamic_cast<Pipe*>(link);
                    // Check if this node belongs to pipe segments
                    if (node->name.find(pipe->name + "_seg") != std::string::npos) {
                        isInteriorNode = true;
                        break;
                    }
                }
            }

            if (isInteriorNode) {
                // Special handling for interior nodes
                // Interior nodes don't have emitter or demand flows
                node->qGrad = 0.0;
                node->outflow = 0.0;
                node->actualDemand = 0.0;
                // Interior nodes just pass flow through
                q = 0.0;
                dqdh = 0.0;
            }
            else {
                // Original handling for main junctions
                // Contribution from emitter flow
                q = node->findEmitterFlow(h, dqdh);
                node->qGrad += dqdh;
                node->outflow += q;
                xQ[i] -= q;

                // Contribution from demand flow
                if (node->fixedGrade) {
                    q = xQ[i];
                    xQ[i] -= q;
                }
                else {
                    q = node->findActualDemand(nw, h, dqdh);
                    node->qGrad += dqdh;
                    xQ[i] -= q;
                }
                node->actualDemand = q;
                node->outflow += q;
            }
        }


        // ... for tanks and reservoirs all flow excess becomes outflow

        else
        {
            node->outflow = xQ[i];
            xQ[i] = 0.0;
        }
    }
}

//-----------------------------------------------------------------------------

//  Find the error norm in satisfying flow continuity at each node.

double HydBalance::findFlowErrorNorm(double xQ[], Network* nw)
{
// Note: contributions to the nodal flow imbalance array xQ[] were
//       found previously from findHeadErrorNorm() and findNodeOutflows())

    double norm = 0.0;
    maxFlowErr = 0.0;

    int nodeCount = nw->nodes.size();
    for (int i = 0; i < nodeCount; i++)
    {
        // ... update network's max. flow error

        if ( abs(xQ[i]) > maxFlowErr )
        {
            maxFlowErr = abs(xQ[i]);
            maxFlowErrNode = i;
        }

        // ... update sum of squared errors (flow imbalances)

        norm += xQ[i] * xQ[i];
    }

    // ... return sum of squared errors normalized by number of nodes

    return norm / nodeCount;
}

//-----------------------------------------------------------------------------

//  Assign the leakage flow along each network pipe to its end nodes.

void findLeakageFlows(double lamda, double dH[], double xQ[], Network* nw)
{
    double dqdh = 0.0;  // gradient of leakage outflow w.r.t. pressure head

    for (Link* link : nw->links)
    {
        // ... skip links that don't leak

        link->leakage = 0.0;
        dqdh = 0.0;
        if ( !link->canLeak() ) continue;

        // ... identify link's end nodes and their indexes

        Node* node1 = link->fromNode;
        Node* node2 = link->toNode;
        int n1 = node1->index;
        int n2 = node2->index;

        // ... no leakage if neither end node is not a junction

        bool canLeak1 = (node1->type() == Node::JUNCTION);
        bool canLeak2 = (node2->type() == Node::JUNCTION);
        if ( !canLeak1 && !canLeak2 ) continue;

        // ... find link's average pressure head

        double h1 = node1->head + lamda * dH[n1] - node1->elev;
        double h2 = node2->head + lamda * dH[n2] - node2->elev;
        double h = (h1 + h2) / 2.0;
        if ( h <= 0.0 ) continue;

        // ... find leakage and its gradient

        link->leakage = link->findLeakage(nw, h, dqdh);

        // ... split leakage flow between end nodes, unless one cannot
        //     support leakage or has negative pressure head

        double q = link->leakage / 2.0;
        if ( h1 * h2 <= 0.0 || canLeak1 * canLeak2 == 0 ) q = 2.0 * q;

        // ... add leakage to each node's outflow

        if ( h1 > 0.0 && canLeak1 )
        {
            node1->outflow += q;
            node1->qGrad += dqdh;
            xQ[n1] -= q;
        }
        if ( h2 > 0.0 && canLeak2 )
        {
            node2->outflow += q;
            node2->qGrad += dqdh;
            xQ[n2] -= q;
        }
    }
}

//-----------------------------------------------------------------------------

//  Find the sum of all link flow changes relative to the sum of all link flows.

double findTotalFlowChange(double lamda, double dQ[], Network* nw)
{
    double qSum = 0.0;
    double dqSum = 0.0;
    double dq;

    for ( int i = 0; i < nw->count(Element::LINK); i++ )
    {
        Link* link = nw->links[i];

        dq = lamda * dQ[i];
        dqSum += abs(dq);
        qSum += abs(link->flow + dq);
    }
    if ( qSum > 0.0 ) return dqSum / qSum;
    else return dqSum;
}

//  Evaluate the error in satisfying the conservation of flow and energy
//  equations for unsteady (water hammer) analysis.
double HydBalance::evaluateUnsteady(
            double lamda,      // step size
            double dH[],       // change in nodal heads
            double dQ[],       // change in link flows
            double dQ_start[], // change in link start flows
            double dQ_end[],   // change in link end flows
            double xQ[],       // nodal inflow minus outflow
            Network* nw,       // network being analyzed
            double currentTime,   // current time step
            double tstep)      // time step size
{
    // ... initialize which elements have the maximum errors
    maxFlowErr = 0.0;
    maxHeadErr = 0.0;
    maxCharErr = 0.0;  // Add characteristic equation error tracking
    maxFlowChange = 0.0;
    maxHeadErrLink = -1;
    maxFlowErrNode = -1;
    maxFlowChangeLink = -1;
    
    // ... initialize nodal flow imbalances to 0
    int nodeCount = nw->count(Element::NODE);
    memset(&xQ[0], 0, nodeCount*sizeof(double));
    
    // ... find the error norm in satisfying conservation of energy and momentum
    //     (updating xQ with internal link flows)
    double norm = findUnsteadyHeadErrorNorm(lamda, dH, dQ, dQ_start, dQ_end, xQ, nw, currentTime, tstep);
    
    // ... update xQ with external outflows (including proper tank handling)
    findUnsteadyNodeOutflows(lamda, dH, xQ, nw, tstep);
    
    // ... add the error norm in satisfying conservation of flow
    norm += findUnsteadyFlowErrorNorm(xQ, nw);

    // ... evaluate the total relative flow change
    totalFlowChange = findUnsteadyTotalFlowChange(lamda, dQ, dQ_start, dQ_end, nw);

    // ... return the root mean square error
    return sqrt(norm);
}

//  Find the error norm in satisfying the characteristic equations and momentum
//  conservation for unsteady flow analysis
double HydBalance::findUnsteadyHeadErrorNorm(
    double lamda, double dH[], double dQ[], double dQ_start[], double dQ_end[], 
    double xQ[], Network* nw, double currentTime, double tstep)
{
    double norm = 0.0;
    double count = 0.0;
    maxHeadErr = 0.0;
    maxCharErr = 0.0;
    maxFlowChange = 0.0;
    maxFlowChangeLink = 0;

    int nodeCount = nw->nodes.size();
    int linkCount = nw->count(Element::LINK);

    // Create and initialize the vector tracking nodes connected to incompressible links
    std::vector<bool> nodeHasIncompressibleLink(nodeCount, false);

    std::vector<int> nodeDomain(nodeCount, -1);
    identifyHydraulicDomains(nodeDomain, nw);
    
    // Identify nodes connected to incompressible flow links
    for (int i = 0; i < nodeCount; i++) {
        Node* node = nw->nodes[i];
        if (!node) continue;
        
        // Check all links connected to this node
        for (Link* link : nw->getLinksConnectedToNode(node)) {
            // Verify the link is connected to this node (redundant but safe check)
            if (link->fromNode == node || link->toNode == node) {
                // Check if this is an incompressible flow link
                if (link->flowModel != Link::COMPRESSIBLE) // && link->status != Link::LINK_CLOSED) {
                {    
                    nodeHasIncompressibleLink[i] = true;
                    //break; 
                }
            }
        }
    }

    for (int i = 0; i < linkCount; i++)
    {
        // === IDENTIFY LINK'S END NODES AND GET UPDATED STATES ===
        Link* link = nw->link(i);
        int n1 = link->fromNode->index;
        int n2 = link->toNode->index;
        
        // Get updated heads including the proposed changes from the current iteration
        double h1 = link->fromNode->head + lamda * dH[n1];
        double h2 = link->toNode->head + lamda * dH[n2];

        if (link->flowModel == Link::COMPRESSIBLE) {
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (pipe) {
                // Check if either node is an "intermediate node" connected to incompressible links
                bool node1IsIntermediate = nodeHasIncompressibleLink[n1];
                bool node2IsIntermediate = nodeHasIncompressibleLink[n2];

                double avgFlow = (pipe->startFlow + dQ_start[i] + pipe->endFlow + dQ_end[i]) / 2.0;

                // For node1: Either use startFlow or average flow based on connectivity
                if (node1IsIntermediate) // && abs(avgFlow) < 0.001) 
                {
                    // Use average flow for intermediate nodes
                    xQ[n1] -= pipe->startFlow;
                } else {
                    // Use startFlow for pure compressible nodes
                    xQ[n1] -= pipe->startFlow;
                }
                
                // For node2: Either use endFlow or average flow based on connectivity
                if (node2IsIntermediate) // && abs(avgFlow) < 0.001) 
                {
                    // Use average flow for intermediate nodes
                    xQ[n2] += pipe->endFlow;
                } else {
                    // Use endFlow for pure compressible nodes
                    xQ[n2] += pipe->endFlow;
                } // */
            } else {
                // Non-pipe compressible link (fallback)
                xQ[n1] -= link->flow;
                xQ[n2] += link->flow;
            }
        } else {
            // Incompressible flow links always use single flow value
            xQ[n1] -= link->flow;
            xQ[n2] += link->flow;
        }

        
        // === EVALUATE BASED ON FLOW MODEL TYPE ===
        if (link->flowModel == Link::COMPRESSIBLE)
        {
            // === COMPRESSIBLE FLOW MODEL (WATER HAMMER ANALYSIS) ===
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue; // Skip non-pipe links in compressible flow analysis
            
            // Get updated flows including the proposed changes from the current iteration
            double startFlow = pipe->startFlow + lamda * dQ_start[i];
            double endFlow = pipe->endFlow + lamda * dQ_end[i];

            // Update node flow balances based on flow model and node connectivity
        
            // === CALCULATE CHARACTERISTIC PARAMETERS ===
            
            // Physical constants
            double c = pipe->getWaveSpeed();
            double A = pipe->getArea();
            double L = pipe->length;
            double g = GRAVITY;
            double Bj = c / (g * A);  // Wave impedance term (B')
            double D = pipe->diameter;
            double epsilon = 0;  // Friction integration parameter
            
            // CORRECTED: Calculate proper reach-back time based on discretization
            double reachBackTime;
            
            // Calculate wave travel time and reach-back time
            double physicalWaveTime = L / c;
                    
            // Determine reach-back index based on number of reaches
            int reachBackSteps = pipe->numReaches;
            double waveTravel = reachBackSteps * tstep;  // Wave travel distance

            // Retrieve reach-back flows and heads from history
            double Qa = pipe->pastStartFlow;  // Default fallback
            double Qb = pipe->pastEndFlow;    // Default fallback
            double Ha = pipe->fromNode->pastHead;  // Default fallback
            double Hb = pipe->toNode->pastHead;    // Default fallback
                    
            // Use the FlowHistoryManager to get reach-back values
            FlowHistoryResult historyResult = FlowHistoryManager::getInstance().getReachBackValues(pipe, currentTime, physicalWaveTime, nw);
                    
            // If history was found, use the values from the result
            if (historyResult.found) {
                Qa = historyResult.startFlow;
                Qb = historyResult.endFlow;
                Ha = historyResult.startHead;
                Hb = historyResult.endHead;
            } 

            // Convert gradients to resistance coefficients for characteristic equations
            double KA, KB;

            if (nw->option(Options::HEADLOSS_MODEL) == "D-W")
            {
                double frictionFactor = pipe->getFrictionFactor(pipe->getRe(pipe->flow, nw->option(Options::KIN_VISCOSITY)));
                double K = frictionFactor * pipe->length / (2 * g * D * A * A);

                double frictionFactorA = pipe->getFrictionFactor(pipe->getRe(Qa, nw->option(Options::KIN_VISCOSITY)));
                double frictionFactorB = pipe->getFrictionFactor(pipe->getRe(Qb, nw->option(Options::KIN_VISCOSITY)));
                // Split K for each end as needed (KA for upstream, KB for downstream)
                KA = frictionFactorA * pipe->length / (2 * g * D * A * A);
                KB = frictionFactorB * pipe->length / (2 * g * D * A * A);

            }
            else
            {
                KA = pipe->resistance;
                KB = pipe->resistance;
            }


            // === ENSURE NUMERICAL STABILITY ===
            const double MIN_RESISTANCE = 1e-10;
            KA = std::max(KA, MIN_RESISTANCE);
            KB = std::max(KB, MIN_RESISTANCE);
            
            // === UNSTEADY FRICTION TERMS ===
            double kUFa = pipe->computeShearDecayCoeff(Qa);
            double kUFb = pipe->computeShearDecayCoeff(Qb);
            double Ua = 0.0;
            double Ub = 0.0;
            
            if (kUFa > 0.0 && kUFb > 0.0) {
                auto sign = [](double val) -> double { 
                    return (val > 0.0) ? 1.0 : ((val < 0.0) ? -1.0 : 0.0); 
                };
                
                // CORRECTED: Unsteady friction terms
                Ua = 0.5 * Bj * kUFb * (Qb - 2.0 * sign(Qa) * std::abs(Qb - Qa));
                Ub = 0.5 * Bj * kUFa * (Qa - 2.0 * sign(Qb) * std::abs(Qa - Qb));
            }
            
            // CORRECTED: Characteristic equation coefficients
            double BA, BB, RA, RB;
            
            if (nw->option(Options::HEADLOSS_MODEL) == "D-W") {
                // For boundary A: using C+ from point b
                BA = Bj + epsilon * KB * std::abs(Qb) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * std::abs(Qb)) * Qb - Hb + Ub;
                
                // For boundary B: using C- from point a
                BB = Bj + epsilon * KA * std::abs(Qa) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * std::abs(Qa)) * Qa + Ha + Ua;
            }
            else if (nw->option(Options::HEADLOSS_MODEL) == "H-W") {
                // For boundary A: using C+ from point b
                BA = Bj + epsilon * KB * pow(std::abs(Qb), 0.852) + 0.5 * Bj * kUFb;
                RA = (Bj - (1.0 - epsilon) * KB * pow(std::abs(Qb), 0.852)) * Qb - Hb + Ub;
                
                // For boundary B: using C- from point a
                BB = Bj + epsilon * KA * pow(std::abs(Qa), 0.852) + 0.5 * Bj * kUFa;
                RB = (Bj - (1.0 - epsilon) * KA * pow(std::abs(Qa), 0.852)) * Qa + Ha + Ua;
            }
            
            // Apply wave attenuation scheme if enabled (placeholder for future enhancement)
            if (pipe->impingingWaveCount > 0) {
                // Implement wave attenuation corrections if needed
            }
            
            // === COMPUTE CHARACTERISTIC EQUATION RESIDUALS ===

            double charErrStart = 0.0;  // Error in upstream characteristic equation
            double charErrEnd = 0.0;    // Error in downstream characteristic equation
            
            // === HANDLE BOUNDARY CONDITIONS FOR UPSTREAM NODE (A) ===
            if (link->fromNode->type() == Node::RESERVOIR)
            {
                // For reservoir boundary: head is fixed, so check characteristic equation
                // Positive characteristic: BA * Q - H - RA = 0
                charErrStart = BA * startFlow - h1 - RA;
            }
            else if (link->fromNode->type() == Node::TANK)
            {
                // For tank boundary: head varies with storage, but characteristic still applies
                Tank* tank = static_cast<Tank*>(link->fromNode);
                double tankArea = tank->area;
                
                // The characteristic equation residual
                charErrStart = BA * startFlow - h1 - RA;
            }
            else
            {
                // Normal junction: apply positive characteristic equation
                // Positive characteristic: BA * Q - H - RA = 0
                charErrStart = BA * startFlow - h1 - RA;
            }
            
            // === HANDLE BOUNDARY CONDITIONS FOR DOWNSTREAM NODE (B) ===
            if (link->toNode->type() == Node::RESERVOIR)
            {
                // For reservoir boundary: head is fixed, check negative characteristic
                // Negative characteristic: BB * Q + H - RB = 0
                charErrEnd = BB * endFlow + h2 - RB;
            }
            else if (link->toNode->type() == Node::TANK)
            {
                // For tank boundary: similar to reservoir but with variable head
                Tank* tank = static_cast<Tank*>(link->toNode);
                double tankArea = tank->area;
                
                // The characteristic equation residual
                charErrEnd = BB * endFlow + h2 - RB;
            }
            else
            {
                // Normal junction: apply negative characteristic equation
                // Negative characteristic: BB * Q + H - RB = 0
                charErrEnd = BB * endFlow + h2 - RB;
            }

            // Track flow changes for convergence monitoring
            double flowChangeStart = lamda * dQ_start[i];
            if (abs(flowChangeStart) > maxFlowChange)
            {
                maxFlowChange = abs(flowChangeStart);
                maxFlowChangeLink = i;
            }

            double flowChangeEnd = lamda * dQ_end[i];
            if (abs(flowChangeEnd) > maxFlowChange)
            {
                maxFlowChange = abs(flowChangeEnd);
                maxFlowChangeLink = i;
            }

            //link->hLoss = h1 - h2;  // Update head loss based on current heads
            
            // === TRACK MAXIMUM CHARACTERISTIC ERROR FOR CONVERGENCE MONITORING ===
            double maxCharErrLink = std::max(abs(charErrStart), abs(charErrEnd));
            if (maxCharErrLink > maxCharErr)
            {
                maxCharErr = maxCharErrLink;
                maxCharErrLink = i;  // Store which link has the maximum error
            }
            
            // === ADD TO OVERALL ERROR NORM ===
            // Square the residuals and add them to the norm (least squares approach)
            norm += charErrStart * charErrStart + charErrEnd * charErrEnd;
            count += 2.0; // Two characteristic equations per pipe (positive and negative)
        }
        else
        {
            // === INCOMPRESSIBLE FLOW (RWC OR QUASI-STEADY) ===
            // This section remains unchanged as it handles steady-state-like flow analysis
            double flow = link->flow + lamda * dQ[i];
            
            // Track flow changes for convergence monitoring
            double flowChange = lamda * dQ[i];
            if (abs(flowChange) > maxFlowChange)
            {
                maxFlowChange = abs(flowChange);
                maxFlowChangeLink = i;
            }
            
            // Compute head loss using the link's head loss model
            link->findHeadLoss(nw, flow);
            
            // Include inertial term for Rigid Water Column (RWC) model
            double unsteadyTerm = 0;
            if (link->type() == Link::PIPE)
            {
                // L/gA * (dQ/dt) term for flow inertia effects
                unsteadyTerm = 0; // link->inertialTerm * (flow - link->pastFlow) / tstep;
            }
            
            if (link->hGrad == 0.0) link->hLoss = h1 - h2;

            // Calculate head loss error: energy equation residual
            // unsteadyTerm - h1 + h2 + hLoss = 0 for perfect solution
            double err = unsteadyTerm - h1 + h2 + link->hLoss;

            // Handle closed links (commented out in original, preserved)
            if (link->status == Link::LINK_CLOSED)
            {
                err = 0.0; // No error for closed links
            }
            
            // Track maximum head error for convergence monitoring
            if (abs(err) > maxHeadErr)
            {
                maxHeadErr = abs(err);
                maxHeadErrLink = i;
            }

            // Add to error norm
            norm += err * err;
            count += 1.0; // One energy equation per incompressible link
        }
    }
    
    if (count == 0.0) return 0;
    else return norm / count;  // Root mean square error across all equations
}


//  Find net external outflow at each network node for unsteady analysis
void HydBalance::findUnsteadyNodeOutflows(double lamda, double dH[], double xQ[], Network* nw, double tstep)
{
    // ... initialize node outflows and their gradients w.r.t. head
    for (Node* node : nw->nodes)
    {
        node->outflow = 0.0;
        node->qGrad = 0.0;
    }

    // ... find pipe leakage flows & assign them to node outflows
    if (nw->leakageModel) findLeakageFlows(lamda, dH, xQ, nw);

    // ... add emitter flows and demands to node outflows
    int nodeCount = nw->count(Element::NODE);
    for (int i = 0; i < nodeCount; i++)
    {
        Node* node = nw->node(i);
        double h = node->head + lamda * dH[i];
        double q = 0.0;
        double dqdh = 0.0;

        // ... for junctions, outflow depends on head
        if (node->type() == Node::JUNCTION)
        {
            // ... contribution from emitter flow
            q = node->findEmitterFlow(h, dqdh);
            node->qGrad += dqdh;
            node->outflow += q;
            xQ[i] -= q;
            double HASTIRLAN;
            // ... contribution from demand flow
            // ... for fixed grade junction, demand is remaining flow excess
            if (node->fixedGrade)
            {
                q = xQ[i];
                xQ[i] -= q;
            }
            // ... otherwise junction has pressure-dependent demand
            else
            {
                q = node->findActualDemand(nw, h, dqdh);
                node->qGrad += dqdh;
                xQ[i] -= q;
                HASTIRLAN = xQ[i];
            }
            node->actualDemand = q;
            node->outflow += q;
        }
        // ... for tanks and reservoirs all flow excess becomes outflow
        // (same handling as in the original EPANET function)
        else
        {
            node->outflow = xQ[i];
            xQ[i] = 0.0;

            
        }
    }
}

double HydBalance::findUnsteadyFlowErrorNorm(double xQ[], Network* nw) {
    double norm = 0.0;
    maxFlowErr = 0.0;
    maxFlowErrNode = -1;
    
    // Get basic network information
    int nodeCount = nw->nodes.size();
    
    // Step 1: Identify hydraulic domains created by closed valves
    std::vector<int> nodeDomain(nodeCount, -1);
    identifyHydraulicDomains(nodeDomain, nw);
    
    // Step 2: Identify nodes that are directly connected to closed valves
    // These nodes require special treatment in error assessment
    std::vector<bool> nodeConnectedToClosedValve(nodeCount, false);
    
    for (int i = 0; i < nodeCount; i++) {
        Node* node = nw->nodes[i];

        if (!node || node->fixedGrade) continue;
        
        // Check all links connected to this node
        for (Link* link : nw->getLinksConnectedToNode(node)) {
            if (link->status == Link::LINK_CLOSED) // || abs(link->flow) < 1e-3) 
            {
                // This node is directly connected to a closed valve
                nodeConnectedToClosedValve[i] = true;
                break; // No need to check further links for this node
            } 
        }
    } // */
    
    // Step 3: Calculate domain-specific errors (excluding closed valve nodes)
    std::map<int, double> domainTotalError;
    std::map<int, int> domainNodeCount;
    std::map<int, int> domainActiveNodeCount; // Nodes not connected to closed valves
    
    int totalActiveNodes = 0; // Counter for nodes included in error calculation
    
    for (int i = 0; i < nodeCount; i++) {
        Node* node = nw->nodes[i];
        if (!node || node->fixedGrade) continue;
        //if (!node) continue; // Skip null nodes

        int domain = nodeDomain[i];
        
        // Always track basic domain statistics
        domainTotalError[domain] += xQ[i] * xQ[i];
        domainNodeCount[domain]++;
        
        // Skip nodes connected to closed valves for active error calculation
        if (nodeConnectedToClosedValve[i]) {
            // Log this exclusion for transparency
            if (abs(xQ[i]) > 0.1) {
                nw->msgLog << "\nExcluding node " << node->name 
                          << " from error calculation (connected to closed valve, error=" 
                          << abs(xQ[i]) << ")";
            }
            continue; // Skip this node in active error calculations
        } // */
        
        // Count active nodes for this domain
        domainActiveNodeCount[domain]++;
        totalActiveNodes++;
        
        // Track global maximum error (only among active nodes)
        if (abs(xQ[i]) > maxFlowErr) {
            maxFlowErr = abs(xQ[i]);
            maxFlowErrNode = i;
        }
        
        // Add to overall error norm (only active nodes)
        norm += xQ[i] * xQ[i];
    }
    
    // Step 4: Provide comprehensive diagnostic information
    if (maxFlowErr > 0.1) {
        nw->msgLog << "\n--- Flow Error Analysis by Hydraulic Domain ---";
        
        for (auto& pair : domainTotalError) {
            int domain = pair.first;
            int totalNodes = domainNodeCount[domain];
            int activeNodes = domainActiveNodeCount[domain];
            
            if (activeNodes > 0) {
                double avgDomainErrorActive = pair.second / activeNodes;
                
                nw->msgLog << "\nDomain " << domain << ": " << activeNodes 
                          << " active nodes (out of " << totalNodes << " total)";
                
                // Only report significant errors in active nodes
                if (avgDomainErrorActive > 0.01) {
                    nw->msgLog << "  Average continuity error: " << sqrt(avgDomainErrorActive);
                }
            } else {
                nw->msgLog << "\nDomain " << domain << ": All " << totalNodes 
                          << " nodes connected to closed valves (excluded from analysis)";
            }
        }
        
        nw->msgLog << "\nTotal active nodes in error calculation: " << totalActiveNodes 
                  << " out of " << nodeCount;
    }
    
    // Step 5: Return normalized error based only on hydraulically active nodes
    // This gives a more accurate picture of simulation health
    if (totalActiveNodes > 0) {
        return norm / totalActiveNodes;
    } else {
        // Special case: if all nodes are connected to closed valves
        nw->msgLog << "\nWarning: All nodes connected to closed valves - cannot assess flow error";
        return 0.0;
    }
}

void HydBalance::identifyHydraulicDomains(std::vector<int>& nodeDomain, Network* nw) {
    int domainCount = 0;
    
    // Reset all domain assignments
    for (int i = 0; i < nodeDomain.size(); i++) {
        nodeDomain[i] = -1;
    }
    
    // Use a breadth-first search to identify connected domains
    for (int i = 0; i < nw->nodes.size(); i++) {
        if (nodeDomain[i] >= 0) continue;  // Already assigned
        
        // Start a new domain
        nodeDomain[i] = domainCount++;
        
        // Use BFS to find all connected nodes
        std::queue<int> nodeQueue;
        nodeQueue.push(i);
        
        while (!nodeQueue.empty()) {
            int currentNodeIdx = nodeQueue.front();
            nodeQueue.pop();
            
            Node* currentNode = nw->nodes[currentNodeIdx];
            if (!currentNode) continue;
            
            // Check all connected links
            for (Link* link : nw->getLinksConnectedToNode(currentNode)) {
                // Skip closed links - they create domain boundaries
                if (link->status == Link::LINK_CLOSED) continue;
                
                // Find the node on the other side
                Node* otherNode = (link->fromNode == currentNode) ? 
                                   link->toNode : link->fromNode;
                
                if (!otherNode) continue;
                
                int otherNodeIdx = otherNode->index;
                
                // If not yet assigned, add to this domain
                if (nodeDomain[otherNodeIdx] < 0) {
                    nodeDomain[otherNodeIdx] = nodeDomain[currentNodeIdx];
                    nodeQueue.push(otherNodeIdx);
                }
            }
        }
    }
}

//  Find the sum of all link flow changes relative to the sum of all link flows
//  for unsteady analysis, considering both start and end flows
double HydBalance::findUnsteadyTotalFlowChange(double lamda, double dQ[], double dQ_start[], 
                                             double dQ_end[], Network* nw)
{
    double qSum = 0.0;
    double dqSum = 0.0;
    double dq;

    for (int i = 0; i < nw->count(Element::LINK); i++)
    {
        Link* link = nw->link(i);
        
        if (link->flowModel == Link::COMPRESSIBLE)
        {
            // For compressible flow, consider both start and end flows
            Pipe* pipe = dynamic_cast<Pipe*>(link);
            if (!pipe) continue;
            
            double dq_s = lamda * dQ_start[i];
            double dq_e = lamda * dQ_end[i];
            
            dqSum += abs(dq_s) + abs(dq_e);
            qSum += abs(pipe->startFlow + dq_s) + abs(pipe->endFlow + dq_e);
        }
        else
        {
            // For incompressible flow
            dq = lamda * dQ[i];
            dqSum += abs(dq);
            qSum += abs(link->flow + dq);
        }
    }
    
    if (qSum > 0.0) return dqSum / qSum;
    else return dqSum;


} 


