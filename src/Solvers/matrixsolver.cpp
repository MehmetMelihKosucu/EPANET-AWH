/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

#include "matrixsolver.h"

// Include headers for the different matrix solvers here
#include "sparspaksolver.h"
//#include "cholmodsolver.h"

using namespace std;

MatrixSolver::MatrixSolver() {}

MatrixSolver::~MatrixSolver() {}

MatrixSolver* MatrixSolver::factory(const string name, ostream& logger)
{
    //if (name == "CHOLMOD") return new CholmodSolver();
    if (name == "SPARSPAK") return new SparspakSolver(logger);
    return nullptr;
}

// Function to set coupling with IT1 nodes
void MatrixSolver::setIT1Coupling(int offDiagIndex, Node* it1Node, double coefficient) {
    // For standard matrix solvers, this might just store the coefficient
    // in the off-diagonal matrix
    if (offDiagIndex >= 0 && offDiagIndex < offDiagSize) {
        offDiag[offDiagIndex] = coefficient;
        
        // Ensure it1CouplingNodes is sized properly before accessing
        if (it1CouplingNodes.size() <= offDiagIndex) {
            it1CouplingNodes.resize(offDiagSize, nullptr);
        }
        
        // Store the node for more complex solvers
        it1CouplingNodes[offDiagIndex] = it1Node;
    }
}
