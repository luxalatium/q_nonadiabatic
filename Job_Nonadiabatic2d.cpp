// ==============================================================================
//
//  Job_Nonadiabatic2d.cpp
//  QTR
//
//  Created by Albert Lu on 1/14/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 1/14/21
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Nonadiabatic2d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Nonadiabatic2d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobNonadiabatic2d::JobNonadiabatic2d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    nadb = new Nonadiabatic2d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobNonadiabatic2d::~JobNonadiabatic2d()
{
    delete nadb;
}

/* ------------------------------------------------------------------------------- */

void JobNonadiabatic2d::run(class QTR *qtr)
{     
    nadb->Evolve();
    log->log("[Job_Nonadiabatic2d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

