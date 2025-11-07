/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Licensed under the terms of the MIT License (see the LICENSE file for details).
 *
 */

 //! \file network.h
 //! \brief Describes the Network class.

#ifndef NETWORK_H_
#define NETWORK_H_

#include "Core/options.h"
#include "Core/units.h"
#include "Core/qualbalance.h"
#include "Elements/element.h"
#include "Utilities/graph.h"

#include <vector>
#include <ostream>
#include <unordered_map>
#include <string>

class Node;
class Link;
class Pattern;
class Curve;
class Control;
class HeadLossModel;
class DemandModel;
class LeakageModel;
class QualModel;
class MemPool;

//! \class Network
//! \brief Contains the data elements that describe a pipe network.

class Network {
public:
    Network();
    ~Network();

    // Clears all elements from the network
    void clear();

    // Adds an element to the network
    bool addElement(Element::ElementType eType, int subType, std::string name);

    // Finds element counts by type and index by ID name
    int count(Element::ElementType eType);
    int indexOf(Element::ElementType eType, const std::string& name);

    // Gets an analysis option by type
    int           option(Options::IndexOption type);
    double        option(Options::ValueOption type);
    long          option(Options::TimeOption type);
    std::string   option(Options::StringOption type);


    // Gets a network element by ID name
    Node* node(const std::string& name);
    Link* link(const std::string& name);
    Pattern* pattern(const std::string& name);
    Curve* curve(const std::string& name);
    Control* control(const std::string& name);

    // Gets a network element by index
    Node* node(const int index);
    Link* link(const int index);
    Pattern* pattern(const int index);
    Curve* curve(const int index);
    Control* control(const int index);

    // Creates analysis models
    bool createHeadLossModel();
    bool createDemandModel();
    bool createLeakageModel();
    bool createQualModel();

    // Unit conversions
    double ucf(Units::Quantity quantity);
    std::string getUnits(Units::Quantity quantity);
    void convertUnits();

    // Adds/writes network title
    void addTitleLine(const std::string& line);
    void writeTitle(std::ostream& out);

    // Network graph theory operations
    Graph graph;

    // Computational sub-models
    HeadLossModel* headLossModel;
    DemandModel* demandModel;
    LeakageModel* leakageModel;
    QualModel* qualModel;

    // Additional utility methods
    int nodeCount() const;
    int linkCount() const;
    std::vector<Link*> getLinksConnectedToNode(Node* node);

    // Getter for options (inline functions)
    Options& getOptions();

    //void assignNodeIndices();

    // Elements of a network
    std::vector<std::string> title;         //!< Descriptive title for the network
    std::vector<Node*> nodes;               //!< Collection of node objects
    std::vector<Link*> links;               //!< Collection of link objects
    std::vector<Curve*> curves;             //!< Collection of data curve objects
    std::vector<Pattern*> patterns;         //!< Collection of time pattern objects
    std::vector<Control*> controls;         //!< Collection of control rules

    // Conversion factors and options
    Units units;                            //!< Unit conversion factors
    Options options;                        //!< Analysis options
    QualBalance qualBalance;                //!< Water quality mass balance
    std::ostringstream msgLog;              //!< Status message log

    MemPool* memPool;                       //!< Memory pool for network objects

    std::unordered_map<std::string, Element*> nodeTable;
    std::unordered_map<std::string, Element*> linkTable;

private:

    // Hash tables for element indexing
    
    std::unordered_map<std::string, Element*> curveTable;
    std::unordered_map<std::string, Element*> patternTable;
    std::unordered_map<std::string, Element*> controlTable;
};
    //-----------------------------------------------------------------------------
    //    Inline Functions
    //-----------------------------------------------------------------------------

inline int Network::option(Options::IndexOption type) {
    return options.indexOption(type);
}

inline double Network::option(Options::ValueOption type) {
    return options.valueOption(type);
}

inline long Network::option(Options::TimeOption type) {
    return options.timeOption(type);
}

inline std::string Network::option(Options::StringOption type)
{
    return options.stringOption(type);
}

inline double Network::ucf(Units::Quantity quantity) {
    return units.factor(quantity);
}

inline std::string Network::getUnits(Units::Quantity quantity) {
    return units.name(quantity);
}

inline void Network::addTitleLine(const std::string& line) {
    title.push_back(line);
}
#endif // NETWORK_H_    
//         void Network::writeTitle(std::ostream& out) 