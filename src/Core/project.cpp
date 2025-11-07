/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

 /////////////////////////////////////////////
 //  Implementation of the Project class.  //
 ////////////////////////////////////////////

#include "project.h"
#include "epanet3.h"
#include "Core/error.h"
#include "Core/diagnostics.h"
#include "Core/datamanager.h"
#include "Core/constants.h"
#include "Core/network.h"
#include "Core/hydengine.h"
#include "Core/hydbalance.h"
#include "Elements/pipe.h"
#include "Elements/pump.h"
#include "Elements/valve.h"
#include "Elements/pumpcurve.h"
#include "Elements/control.h"
#include "Elements/junction.h"
#include "Elements/tank.h"
#include "Elements/link.h"
#include "Input/inputreader.h"
#include "Output/projectwriter.h"
#include "Output/reportwriter.h"
#include "Utilities/utilities.h"
#include "linkparser.h"
#include "rwcggasolver.h"
#include "matrixsolver.h"
#include "options.h"
#include <cstring>
#include <cmath>
#include <limits>
#include <iostream>   //for debugging
#include <iomanip>
#include <algorithm>
#include <string>
#include <time.h>
#include <fstream>
#include "segpool.h"
using namespace std;

//-----------------------------------------------------------------------------

namespace Epanet
{
	//  Constructor

	Project::Project() :
		inpFileName(""),
		outFileName(""),
		tmpFileName(""),
		rptFileName(""),
		networkEmpty(true),
		hydEngineOpened(false),
		qualEngineOpened(false),
		outputFileOpened(false),
		solverInitialized(false),
		runQuality(false),
		currentMode(SolverMode::QUASI_STEADY)
	{
		Utilities::getTmpFileName(tmpFileName);

		// Initialize hydEngine here
		
		//hydEngine = std::make_unique<HydEngine>(this); // Use new or smart pointer as needed
	}

	//  Destructor

	Project::~Project()
	{
		//cout << "\nDestructing Project.";

		closeReport();
		outputFile.close();
		remove(tmpFileName.c_str());
		
		//cout << "\nProject destructed.\n";
	}

	//-----------------------------------------------------------------------------

	//  Load a project from a file.

	int Project::load(const char* fname)
	{
		try
		{
			// ... clear any current project
			clear();

			// ... check for duplicate file names
			string s = fname;
			if (s.size() == rptFileName.size() && Utilities::match(s, rptFileName))
			{
				throw FileError(FileError::DUPLICATE_FILE_NAMES);
			}
			if (s.size() == outFileName.size() && Utilities::match(s, outFileName))
			{
				throw FileError(FileError::DUPLICATE_FILE_NAMES);
			}

			// ... save name of input file
			inpFileName = fname;

			// ... use an InputReader to read project data from the input file
			InputReader inputReader;
			inputReader.readFile(fname, &network);
			networkEmpty = false;
			runQuality = network.option(Options::QUAL_TYPE) != Options::NOQUAL;

			// ... convert all network data to internal units
			network.convertUnits();
			network.getOptions().adjustOptions();

			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Save the project to a file.

	int Project::save(const char* fname)
	{
		try
		{
			if (networkEmpty) return 0;
			ProjectWriter projectWriter;
			projectWriter.writeFile(fname, &network);
			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Clear the project of all data.

	void Project::clear()
	{
		hydEngine.close();
		hydEngineOpened = false;

		qualEngine.close();
		qualEngineOpened = false;

		network.clear();
		networkEmpty = true;

		solverInitialized = false;
		inpFileName = "";
	}

	//-----------------------------------------------------------------------------

	//  Initialize the project's solvers.

	int Project::initSolver(bool initFlows)
	{
		try
		{
			if (networkEmpty) return 0;
			solverInitialized = false;
			Diagnostics diagnostics;
			diagnostics.validateNetwork(&network);

			// ... open & initialize the hydraulic engine
			if (!hydEngineOpened)
			{
				initFlows = true;
				hydEngine.open(&network);
				hydEngineOpened = true;
			}
			hydEngine.init(initFlows);

			// ... open and initialize the water quality engine
			if (runQuality == true)
			{
				if (!qualEngineOpened)
				{
					qualEngine.open(&network);
					qualEngineOpened = true;
				}
				qualEngine.init();
			}

			// ... mark solvers as being initialized
			solverInitialized = true;

			// ... initialize the binary output file
			outputFile.initWriter();
			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Solve network hydraulics at the current point in time.

	int Project::runSolver(double* t)
	{
		try
		{
			if (!solverInitialized) throw SystemError(SystemError::SOLVER_NOT_INITIALIZED);
			hydEngine.solve(t);
			if (outputFileOpened && std::fmod(*t, network.option(Options::REPORT_STEP)) == 0)
			{
				outputFile.writeNetworkResults();
			}
			return 0; // */
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Advance the hydraulic solver to the next point in time while updating
	//  water quality.

	int Project::advanceSolver(double* dt)
	{
		try
		{
			
			// Advance the hydraulic engine
			hydEngine.advance(dt);

			// If simulation is complete, finalize results
			if (*dt == 0) {
				finalizeSolver();
			}
			// Otherwise, advance the water quality solver
			else if (runQuality) {
				qualEngine.solve(*dt);
			}

			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Open a binary file that saves computed results.

	int Project::openOutput(const char* fname)
	{
		//... close an already opened output file
		if (networkEmpty) return 0;
		outputFile.close();
		outputFileOpened = false;

		// ... save the name of the output file
		outFileName = fname;
		if (strlen(fname) == 0) outFileName = tmpFileName;

		// ... open the file
		try
		{
			outputFile.open(outFileName, &network);
			outputFileOpened = true;
			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Save results for the current time period to the binary output file.

	int Project::saveOutput()
	{
		if (!outputFileOpened) return 0;
		try
		{
			outputFile.writeNetworkResults();
			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	//  Finalize computed quantities at the end of a run

	void Project::finalizeSolver()
	{
		if (!solverInitialized) return;

		// Save energy usage results to the binary output file.
		if (outputFileOpened)
		{
			double totalHrs = hydEngine.getElapsedTime() / 3600.0;
			double peakKwatts = hydEngine.getPeakKwatts();
			outputFile.writeEnergyResults(totalHrs, peakKwatts);
		}

		// Write mass balance results for WQ constituent to message log
		if (runQuality && network.option(Options::REPORT_STATUS))
		{
			network.qualBalance.writeBalance(network.msgLog);
		}
	}

	//-----------------------------------------------------------------------------

	//  Open the project's status/report file.

	int  Project::openReport(const char* fname)
	{
		try
		{
			//... close an already opened report file
			if (rptFile.is_open()) closeReport();

			// ... check that file name is different from input file name
			string s = fname;
			if (s.size() == inpFileName.size() && Utilities::match(s, inpFileName))
			{
				throw FileError(FileError::DUPLICATE_FILE_NAMES);
			}
			if (s.size() == outFileName.size() && Utilities::match(s, outFileName))
			{
				throw FileError(FileError::DUPLICATE_FILE_NAMES);
			}

			// ... open the report file
			rptFile.open(fname);
			if (!rptFile.is_open())
			{
				throw FileError(FileError::CANNOT_OPEN_REPORT_FILE);
			}
			ReportWriter rw(rptFile, &network);
			rw.writeHeading();
			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	//-----------------------------------------------------------------------------

	// Write a message to the project's message log.

	void  Project::writeMsg(const std::string& msg)
	{
		network.msgLog << msg;
	}

	//-----------------------------------------------------------------------------

	//  Write the project's title and option summary to the report file.

	void Project::writeSummary()
	{
		if (!rptFile.is_open()) return;
		ReportWriter reportWriter(rptFile, &network);
		reportWriter.writeSummary(inpFileName);
	}

	//-----------------------------------------------------------------------------

	//  Close the project's report file.

	void Project::closeReport()
	{
		if (rptFile.is_open()) rptFile.close();
	}

	//-----------------------------------------------------------------------------

	//  Write the project's message log to an output stream.

	void Project::writeMsgLog(ostream& out)
	{
		out << network.msgLog.str();
		network.msgLog.str("");
	}

	//-----------------------------------------------------------------------------

	//  Write the project's message log to the report file.

	void Project::writeMsgLog()
	{
		if (rptFile.is_open())
		{
			rptFile << network.msgLog.str();
			network.msgLog.str("");
		}
	}

	//-----------------------------------------------------------------------------

	//  Write results at the current time period to the report file.

	void Project::writeResults(int t)
	{
		if (!rptFile.is_open()) return;
		ReportWriter reportWriter(rptFile, &network);
		reportWriter.writeResults(t);
	}

	//-----------------------------------------------------------------------------

	//  Write all results saved to the binary output file to a report file.

	int Project::writeReport()
	{
		try
		{
			if (!outputFileOpened)
			{
				throw FileError(FileError::NO_RESULTS_SAVED_TO_REPORT);
			}
			ReportWriter reportWriter(rptFile, &network);
			reportWriter.writeReport(inpFileName, &outputFile);
			return 0;
		}
		catch (ENerror const& e)
		{
			writeMsg(e.msg);
			return e.code;
		}
	}

	// write this function properly
	int Project::getTotalSimulationTime() const {
		
		return Options::TOTAL_DURATION;
	}

	double Project::pressureManagement(double Xm, double error, double delta_Xm, double Xm_Last, int t)
	{
		network.links;

		double ref;
		double pToNode;
		double pRemoteNode;

		for (int j = 0; j < network.count(Element::LINK); j++)
		{
		//	Link* link = getNetwork()->link(j);
			Link* link = network.link(j);
			if (link->type() == Link::VALVE)
			{
				Valve* valve = static_cast<Valve*>(link);
				if (valve->valveType == Valve::DPRV)
				{
					pToNode = valve->toNode->head - valve->toNode->elev;
					
					if (valve->presManagType == Valve::FO)
					{
						ref = valve->fixedOutletPressure / network.ucf(Units::LENGTH);

						error = ref - pToNode;
					}

					else if (valve->presManagType == Valve::TM)
					{
						if (0 <= t && t <= 25200)
							ref = valve->nightPressure / network.ucf(Units::LENGTH);
						else if (25200 <= t && t <= 86400)
							ref = valve->dayPressure / network.ucf(Units::LENGTH);
						else if (86400 <= t && t <= 111600)
							ref = valve->nightPressure / network.ucf(Units::LENGTH);
						else if (111600 <= t && t <= 172800)
							ref = valve->dayPressure / network.ucf(Units::LENGTH);

						error = ref - pToNode;
					}

					else if (valve->presManagType == Valve::FM)
					{
						// Flow Modulated Pressure Control
						double q6;
						int Index6;
						int Error6 = EN_getLinkIndex("6", &Index6, getNetwork());
						double get_Flow6 = EN_getLinkValue(Index6, EN_FLOW, &q6, getNetwork());
						q6;

						ref = (valve->a_FM * (valve->flow * network.ucf(Units::FLOW))* (valve->flow * network.ucf(Units::FLOW)) + valve->b_FM * (valve->flow * network.ucf(Units::FLOW)) + valve->c_FM) / network.ucf(Units::LENGTH);

						error = ref - pToNode;
					}

					else if (valve->presManagType == Valve::RNM)
					{
						pRemoteNode = valve->remoteNode->head - valve->remoteNode->elev;
						ref = valve->rnmPressure / network.ucf(Units::LENGTH);
						error = ref - pRemoteNode;
					}
					valve->dprvOutletPressure = ref;
					double Vcontrol = 0.0047; // 0.1659789336 ft3; // 0.0047 m3
					double lift = 0.057; // 0.187007874 ft; // 0.057 m
					double k5 = 1.30;
					double k6 = 0.56;
					double alfaopen = 1.1 * pow(10, -6); // 1.184030146 * pow(10, -5) ft2/s ; // 1.1 * pow(10, -6) m2/s
					double alfaclose = 10 * pow(10, -6); // 1.076391042 * pow(10, -4) ft2/s; // 10 * pow(10, -6) m2/s
					double Acs = (k5 * Xm * Xm + k6) * Vcontrol / lift; // m2

					double q3;

					// Physical Based Control
					if (error >= 0)
					{
						q3 = alfaopen * error;
					}
					else if (error <= 0)
					{
						q3 = alfaclose * error;
					}

					delta_Xm = q3 / Acs;

					Xm = delta_Xm + Xm_Last;

					if (valve->status == Valve::LINK_OPEN)
						Xm = 1.0;

					if (Xm > 1.0) Xm = 1.0;
					return Xm;
				}
				else
					continue;
			}
			else
				continue;
		}
	} 
}