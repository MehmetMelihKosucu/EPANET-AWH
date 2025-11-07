/* EPANET 3 ALGEBRAIC WATER HAMMER EXTENSION
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

 //! \file options.h
 //! \brief Describes the Options class.

#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "Output/reportfields.h"
#include <string>
#include <sstream>
#include <map>
#include <stdexcept>

class Network;

//! \class Options
//! \brief User-supplied options for analyzing a pipe network.

class Options {
public:
    // Enumerated categories for options
    enum UnitSystem { US, SI };
    enum FlowUnits { CFS, GPM, MGD, IMGD, AFD, LPS, LPM, MLD, CMH, CMD };
    enum PressureUnits { PSI, METERS, KPA };
    enum FileMode { SCRATCH, USE, SAVE };
    enum IfUnbalanced { STOP, CONTINUE };
    enum QualType { NOQUAL, AGE, TRACE, CHEM };
    enum QualUnits { NOUNITS, HRS, PCNT, MGL, UGL };
    enum ReportedItems { NONE, ALL, SOME };

    // String options
    enum StringOption {
        HYD_FILE_NAME, OUT_FILE_NAME, RPT_FILE_NAME, MAP_FILE_NAME,
        HEADLOSS_MODEL, DEMAND_MODEL, LEAKAGE_MODEL, HYD_SOLVER, STEP_SIZING,
        VALVE_REP_TYPE, MATRIX_SOLVER, DEMAND_PATTERN_NAME,
        QUAL_MODEL, QUAL_NAME, QUAL_UNITS_NAME, TRACE_NODE_NAME,
        MAX_STRING_OPTIONS
    };

    // Integer or categorical options
    enum IndexOption {
        UNIT_SYSTEM, FLOW_UNITS, PRESSURE_UNITS, MAX_TRIALS, IF_UNBALANCED,
        HYD_FILE_MODE, DEMAND_PATTERN, ENERGY_PRICE_PATTERN,
        QUAL_TYPE, QUAL_UNITS, TRACE_NODE,
        REPORT_SUMMARY, REPORT_ENERGY, REPORT_STATUS, REPORT_TRIALS,
        REPORT_NODES, REPORT_LINKS,
        MAX_INDEX_OPTIONS
    };

    // Numerical value options
    enum ValueOption {
        SPEC_GRAVITY, KIN_VISCOSITY, DEMAND_MULTIPLIER, MINIMUM_PRESSURE,
        SERVICE_PRESSURE, PRESSURE_EXPONENT, EMITTER_EXPONENT, LEAKAGE_COEFF1,
        LEAKAGE_COEFF2, RELATIVE_ACCURACY, HEAD_TOLERANCE, FLOW_TOLERANCE,
        FLOW_CHANGE_LIMIT, TIME_WEIGHT, TEMP_DISC_PARA, MOLEC_DIFFUSIVITY,
        QUAL_TOLERANCE, BULK_ORDER, WALL_ORDER, TANK_ORDER, BULK_COEFF,
        WALL_COEFF, LIMITING_CONCEN, ROUGHNESS_FACTOR, ENERGY_PRICE,
        PEAKING_CHARGE, PUMP_EFFICIENCY,
        MAX_VALUE_OPTIONS
    };

    // Time-based options
    enum TimeOption {
        START_TIME, HYD_STEP, QUAL_STEP, PATTERN_STEP, PATTERN_START,
        REPORT_STEP, REPORT_START, RULE_STEP, TOTAL_DURATION,
        REPORT_STATISTIC,
        MAX_TIME_OPTIONS
    };

    // Constructor and Destructor
    Options();
    ~Options() {}

    // Getter methods for options
    int flowUnits();
    int pressureUnits() ;
    std::string stringOption(StringOption option);
    int         indexOption(IndexOption option);
    double      valueOption(ValueOption option);
    int         timeOption(TimeOption option);

    // Setter methods for options
    void setDefaults();
    void adjustOptions();
    int setOption(StringOption option, const std::string& value);
    int setOption(IndexOption option, const std::string& value, Network* nw);
    void setOption(IndexOption option, int value);
    void setOption(ValueOption option, double value);
    void setOption(TimeOption option, int value);
    void setReportFieldOption(int fieldType, int fieldIndex, int enabled, int precision, double lowerLimit, double upperLimit);

    // Methods for converting options to strings
    std::string Options::hydOptionsToStr();
    std::string Options::qualOptionsToStr();
    std::string Options::demandOptionsToStr();
    std::string Options::timeOptionsToStr();
    std::string Options::reactOptionsToStr();
    std::string energyOptionsToStr(Network* network);
    std::string Options::reportOptionsToStr();

    // Data members for managing options
    //std::map<TimeOption, int> timeOptions; // Time options stored in a map

private:
    std::string  stringOptions[MAX_STRING_OPTIONS];
    int          indexOptions[MAX_INDEX_OPTIONS];
    double       valueOptions[MAX_VALUE_OPTIONS];
    int          timeOptions[MAX_TIME_OPTIONS];
    ReportFields reportFields;
};

//-----------------------------------------------------------------------------
// Inline Functions
//-----------------------------------------------------------------------------

inline int Options::flowUnits()  {
    return indexOptions[FLOW_UNITS];
}

inline int Options::pressureUnits() {
    return indexOptions[PRESSURE_UNITS];
}

inline std::string Options::stringOption(StringOption option) {
    return stringOptions[option];
}

inline int Options::indexOption(IndexOption option)  {
    return indexOptions[option];
} // */

inline double Options::valueOption(ValueOption option) {
    return valueOptions[option];
}

inline int Options::timeOption(TimeOption option)
{
    return timeOptions[option];
}

#endif // OPTIONS_H_
