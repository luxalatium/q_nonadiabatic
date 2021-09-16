// ==============================================================================
//
//  Job_Nonadiabatic5d.cpp
//  QTR
//
//  Created by Albert Lu on 4/2/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 4/2/21
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Nonadiabatic5d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Nonadiabatic5d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobNonadiabatic5d::JobNonadiabatic5d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    nadb = new Nonadiabatic5d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobNonadiabatic5d::~JobNonadiabatic5d()
{
    delete nadb;
}

/* ------------------------------------------------------------------------------- */

void JobNonadiabatic5d::run(class QTR *qtr)
{     
    nadb->Evolve();
    log->log("[Job_Nonadiabatic5d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

