// ==============================================================================
//
//  Job_Nonadiabatic3d.cpp
//  QTR
//
//  Created by Albert Lu on 1/26/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 1/26/21
//
//  Note:
//
// ==============================================================================

#include <string>

#include "Job_Nonadiabatic3d.h"
#include "Log.h"
#include "Parameters.h"
#include "Qtr.h"
#include "Nonadiabatic3d.h"

using namespace QTR_NS;

/* ------------------------------------------------------------------------------- */

JobNonadiabatic3d::JobNonadiabatic3d(class QTR *q)
{
    qtr = q;
    log = qtr->log;
    parameters = qtr->parameters;
    nadb = new Nonadiabatic3d(qtr);
}
/* ------------------------------------------------------------------------------- */

JobNonadiabatic3d::~JobNonadiabatic3d()
{
    delete nadb;
}

/* ------------------------------------------------------------------------------- */

void JobNonadiabatic3d::run(class QTR *qtr)
{     
    nadb->Evolve();
    log->log("[Job_Nonadiabatic3d] Done! \n"); 
}
/* ------------------------------------------------------------------------------- */

