/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

///////////////////////////////////////////////
// Implementation of EPANET 3.1's API library  //
///////////////'''''///////////////////////////

// TO DO:
// - finish implementing all of the functions declared in EPANET3.H
// - provide a brief comment on what each function does

#include "epanet3.h"
#include "Core/Project.h"
#include "Core/datamanager.h"
#include "Core/constants.h"
#include "Core/error.h"
#include "Core/network.h"
#include "Core/hydengine.h"
#include "Core/hydbalance.h"
#include "Elements/valve.h"
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Elements/pumpcurve.h"
#include "Elements/control.h"
#include "Elements/junction.h"
#include "Elements/tank.h"
#include "Elements/link.h"
#include "Utilities/utilities.h"
#include "rwcggasolver.h"
#include "matrixsolver.h"
#include "linkparser.h"
#include "options.h"
#include <cstring>
#include <cmath>
#include <limits>
#include <iostream>   //for debugging
#include <iomanip>
#include <algorithm>
#include <string>
#include <time.h>



using namespace Epanet;
using namespace std;

#define project(p) ((Project *)p)

extern "C" {

//-----------------------------------------------------------------------------

int EN_getVersion(int* version)
{
    *version = VERSION;
    return 0;
}

//-----------------------------------------------------------------------------

int EN_runEpanet(const char* inpFile, const char* rptFile, const char* outFile) {
    std::cout << "\n... EPANET with AWH Scheme\n";

    Epanet::Project p;
    int err = 0;

    // Open result file
    //std::ofstream keslan("C:\\EPANET-AWH\\Onizuka_1986\\Results-Onizuka1986.txt");

    std::ofstream keslan("C:\\EPANET-AWH\\Networks\\Results-Nault2018-AWH.txt");
    
    if (!keslan.is_open()) {
        std::cerr << "Error: Could not open Results.txt for writing\n";
        //return ERROR_OPENING_FILE;
    }

    int IndexJ1, IndexJ2, IndexJ3, IndexJ4, IndexJ5, IndexJ6, IndexT1;
    double headJ1, headJ2, headJ3, headJ4, headJ5, headJ6, headT1;

    int IndexJ71, IndexJ70, IndexJ874;
    double headJ71, headJ70, headJ874;

	//keslan << "Time" << "\t\t" << "Head_J1_(m)" << "\t\t" << "Head_J2_(m)" << "\t\t" << "Head_J3_(m)" << "\t\t" << "Head_J4_(m)" << "\t\t" << "Head_J5_(m)" << "\t\t" << "Head_J6_(m)" << "\t\t" << "Head_T1_(m)" << "\n";

    keslan << "Time" << "\t\t" << "Head_J71_(m)" << "\t\t" << "Head_J70_(m)" << "\t\t" << "Head_J874_(m)" << "\n";

    clock_t start_t = clock();
    
    double T_end = p.getTotalSimulationTime(); // Total simulation time

    for (;;)
    {
        // ... open the command line files and load network data
        if ((err = p.openReport(rptFile))) break;
        std::cout << "\n    Reading input file ...";
        if ((err = p.load(inpFile))) break;
        if ((err = p.openOutput(outFile))) break;
        p.writeSummary();

        // ... step through each time period
        double t = 0;
        double tstep = 0;
    
        // ... initialize the solver
        std::cout << "\n    Initializing solver ...";
        if ((err = p.initSolver(false))) break;
        std::cout << "\n    ";

        do {
            std::cout << "\r    Solving network at "                     //r
                << Utilities::getTime(static_cast<double>(t + tstep)) << " hrs ...        ";

            err = p.runSolver(&t);
            if (err) break;

            if (!err) err = p.advanceSolver(&tstep);
            if (err) break;

            /*int ErrorJ1 = EN_getNodeIndex("J1", &IndexJ1, &p);
            int ErrorJ2 = EN_getNodeIndex("J2", &IndexJ2, &p);
            int ErrorJ3 = EN_getNodeIndex("J3", &IndexJ3, &p);
            int ErrorJ4 = EN_getNodeIndex("J4", &IndexJ4, &p);
            int ErrorJ5 = EN_getNodeIndex("J5", &IndexJ5, &p);
            int ErrorJ6 = EN_getNodeIndex("J6", &IndexJ6, &p);
            int ErrorT1 = EN_getNodeIndex("T1", &IndexT1, &p);

            double ErrorValJ1 = EN_getNodeValue(IndexJ1, EN_HEAD, &headJ1, &p);
            double ErrorValJ2 = EN_getNodeValue(IndexJ2, EN_HEAD, &headJ2, &p);
            double ErrorValJ3 = EN_getNodeValue(IndexJ3, EN_HEAD, &headJ3, &p);
            double ErrorValJ4 = EN_getNodeValue(IndexJ4, EN_HEAD, &headJ4, &p);
            double ErrorValJ5 = EN_getNodeValue(IndexJ5, EN_HEAD, &headJ5, &p);
            double ErrorValJ6 = EN_getNodeValue(IndexJ6, EN_HEAD, &headJ6, &p);
            double ErrorValT1 = EN_getNodeValue(IndexT1, EN_HEAD, &headT1, &p); // */

            int ErrorJ71 = EN_getNodeIndex("J71", &IndexJ71, &p);
            int ErrorJ70 = EN_getNodeIndex("J70", &IndexJ70, &p);
            int ErrorJ874 = EN_getNodeIndex("J874", &IndexJ874, &p);

            double ErrorValJ71 = EN_getNodeValue(IndexJ71, EN_HEAD, &headJ71, &p);
            double ErrorValJ70 = EN_getNodeValue(IndexJ70, EN_HEAD, &headJ70, &p);
            double ErrorValJ874 = EN_getNodeValue(IndexJ874, EN_HEAD, &headJ874, &p); // */

            // (G) Output results if desired
            keslan << Utilities::getTime(t) << "\t\t" << headJ71 << "\t\t" << headJ70 << "\t\t" << headJ874 << "\n"; 
            //keslan << Utilities::getTime(t) << "\t\t" << headJ1 << "\t\t" << headJ2 << "\t\t" << headJ3 << "\t\t" << headJ4 << "\t\t" << headJ5 << "\t\t" << headJ6 << "\t\t" << headT1 << "\n";
            keslan.flush();
            
            // Output results
            //myfile << Utilities::getTime(t) << " " << Xm << "\n";

        } while (tstep > 0 && !err);
        
        break;
    }
    

    // ... simulation was successful
    if (!err)
    {
        // ... report execution time
        clock_t end_t = clock();
        double cpu_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
        std::stringstream ss;
        ss << "\n  Simulation completed in ";
        p.writeMsg(ss.str());
        ss.str("");
        if (cpu_t < 0.001) ss << "< 0.001 sec.";
        else ss << std::setprecision(3) << cpu_t << " sec.";
        p.writeMsg(ss.str());

        // ... report simulation results
        std::cout << "\n    Writing report ...                           ";
        err = p.writeReport();
        std::cout << "\n    Simulation completed.                         \n";
        std::cout << "\n... EPANET completed in " << ss.str() << "\n"; //
    }

    if (err)
    {
        p.writeMsgLog();
        std::cout << "\n\n    There were errors. See report file for details.\n";
        return err;
    }
    return 0;
}


//-----------------------------------------------------------------------------

EN_Project EN_createProject()
{
    Epanet::Project* p = new Epanet::Project();
    return (EN_Project *)p;
}

//-----------------------------------------------------------------------------

int EN_deleteProject(EN_Project p)
{
    delete (Epanet::Project *)p;
    return 0;
}

//-----------------------------------------------------------------------------

int EN_loadProject(const char* fname, EN_Project p)
{
    return project(p)->load(fname);
}

//-----------------------------------------------------------------------------

int EN_saveProject(const char* fname, EN_Project p)
{
    return project(p)->save(fname);
}

//-----------------------------------------------------------------------------

int EN_clearProject(EN_Project p)
{
    project(p)->clear();
    return 0;
}

//-----------------------------------------------------------------------------

////////////////////////////////////////////////////////////////
//  NOT SURE IF THIS METHOD WORKS CORRECTLY -- NEEDS TESTING  //
////////////////////////////////////////////////////////////////
int EN_cloneProject(EN_Project pClone, EN_Project pSource)
{
    if ( pSource == nullptr || pClone == nullptr ) return 102;
    int err = 0;
    std::string tmpFile;
    if ( Utilities::getTmpFileName(tmpFile) )
    {
        try
        {
            EN_saveProject(tmpFile.c_str(), pSource);
            EN_loadProject(tmpFile.c_str(), pClone);
        }
        catch (ENerror const& e)
        {
            static_cast<Epanet::Project*>(pSource)->writeMsg(e.msg);
            err = e.code;
  	    }
        catch (...)
        {
            err = 208; //Unspecified error
        }
        if ( err > 0 )
        {
            EN_clearProject(pClone);
        }
        remove(tmpFile.c_str());
        return err;
    }
    return 208;
}

//-----------------------------------------------------------------------------

/*int EN_runProject(EN_Project p)    // <<=============  TO BE COMPLETED
{
    return 0;
}

//-----------------------------------------------------------------------------

int EN_initSolver(int initFlows, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->initSolver(initFlows);
}

//-----------------------------------------------------------------------------

int EN_runSolver(double* t, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->runSolver(t);
}

//-----------------------------------------------------------------------------

int EN_advanceSolver(int *dt, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->advanceSolver(dt);
}

//-----------------------------------------------------------------------------

int EN_openOutputFile(const char* fname, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->openOutput(fname);
}

//-----------------------------------------------------------------------------

int EN_saveOutput(EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->saveOutput();
}

//-----------------------------------------------------------------------------

int EN_openReportFile(const char* fname, EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->openReport(fname);
}

//-----------------------------------------------------------------------------

int EN_writeReport(EN_Project p)
{
    return static_cast<Epanet::Project*>(p)->writeReport();
}

//-----------------------------------------------------------------------------

int EN_writeSummary(EN_Project p)
{
    static_cast<Epanet::Project*>(p)->writeSummary();
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeResults(int t, EN_Project p)
{
    static_cast<Epanet::Project*>(p)->writeResults(t);
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeMsgLog(EN_Project p)
{
    static_cast<Epanet::Project*>(p)->writeMsgLog();
    return 0;
} // */

int EN_runProject(EN_Project p)    // <<=============  TO BE COMPLETED
{
    return 0;
}

//-----------------------------------------------------------------------------

int EN_initSolver(int initFlows, EN_Project p)
{
    return project(p)->initSolver(initFlows);
}

//-----------------------------------------------------------------------------

int EN_runSolver(double* t, EN_Project p)
{
    return project(p)->runSolver(t);
}

//-----------------------------------------------------------------------------

int EN_advanceSolver(double* dt, EN_Project p)
{
    return project(p)->advanceSolver(dt);
}

//-----------------------------------------------------------------------------

int EN_openOutputFile(const char* fname, EN_Project p)
{
    return project(p)->openOutput(fname);
}

//-----------------------------------------------------------------------------

int EN_saveOutput(EN_Project p)
{
    return project(p)->saveOutput();
}

//-----------------------------------------------------------------------------

int EN_openReportFile(const char* fname, EN_Project p)
{
    return project(p)->openReport(fname);
}

//-----------------------------------------------------------------------------

int EN_writeReport(EN_Project p)
{
    return project(p)->writeReport();
}

//-----------------------------------------------------------------------------

int EN_writeSummary(EN_Project p)
{
    project(p)->writeSummary();
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeResults(int t, EN_Project p)
{
    project(p)->writeResults(t);
    return 0;
}

//-----------------------------------------------------------------------------

int EN_writeMsgLog(EN_Project p)
{
    project(p)->writeMsgLog();
    return 0;
}

//-----------------------------------------------------------------------------

int EN_getCount(int element, int* result, EN_Project p)
{
    return DataManager::getCount(element, result, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeIndex(char* name, int* index, EN_Project p)
{
    return DataManager::getNodeIndex(name, index, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeId(int index, char* id, EN_Project p)
{
    return DataManager::getNodeId(index, id, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeType(int index, int* type, EN_Project p)
{
    return DataManager::getNodeType(index, type, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getNodeValue(int index, int param, double* value, EN_Project p)
{
    return DataManager::getNodeValue(index, param, value, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkIndex(char* name, int* index, EN_Project p)
{
    return DataManager::getLinkIndex(name, index, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkId(int index, char* id, EN_Project p)
{
    return DataManager::getLinkId(index, id, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkType(int index, int* type, EN_Project p)
{
    return DataManager::getLinkType(index, type, project(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkNodes(int index, int* fromNode, int* toNode, EN_Project p)
{
    return DataManager::getLinkNodes(index, fromNode, toNode,
        static_cast<Epanet::Project*>(p)->getNetwork());
}

//-----------------------------------------------------------------------------

int EN_getLinkValue(int index, int param, double* value, EN_Project p)
{
   return DataManager::getLinkValue(index, param, value, static_cast<Epanet::Project*>(p)->getNetwork());
}

int EN_setLinkValue(int index, int param, double value, EN_Project p)
{
	return DataManager::setLinkValue(index, param, value, static_cast<Epanet::Project*>(p)->setNetwork());
}

int       EN_getTimeParam(int, long*, EN_Project)
{
    return 0;
}
}  // end of namespace
