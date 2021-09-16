// ==============================================================================
//
//  Job_Nonadiabatic4d.cpp
//  QTR
//
//  Created by Albert Lu on 3/12/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 3/12/21
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Nonadiabatic4d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Nonadiabatic4d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobNonadiabatic4d::JobNonadiabatic4d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    nadb = new Nonadiabatic4d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobNonadiabatic4d::~JobNonadiabatic4d()
{
    delete nadb;
}

/* ------------------------------------------------------------------------------- */

void JobNonadiabatic4d::run(class QTR *qtr)
{     
    nadb->Evolve();
    log->log("[Job_Nonadiabatic4d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

