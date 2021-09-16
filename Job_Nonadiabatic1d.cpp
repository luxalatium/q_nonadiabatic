// ==============================================================================
//
//  Job_Nonadiabatic1d.cpp
//  QTR
//
//  Created by Albert Lu on 11/2/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 11/2/20
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Nonadiabatic1d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Nonadiabatic1d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobNonadiabatic1d::JobNonadiabatic1d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    nadb = new Nonadiabatic1d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobNonadiabatic1d::~JobNonadiabatic1d()
{
    delete nadb;
}

/* ------------------------------------------------------------------------------- */

void JobNonadiabatic1d::run(class QTR *qtr)
{     
    nadb->Evolve();
    log->log("[Job_Nonadiabatic1d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

