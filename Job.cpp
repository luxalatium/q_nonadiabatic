// ==============================================================================
//
//  Job.cpp
//  QTR
//
//  Created by Albert Lu on 11/2/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 4/2/21
//
//  Note:
//
// ==============================================================================

# include "Job.h"
# include "Parameters.h"

# include "Job_Nonadiabatic1d.h"
# include "Job_Nonadiabatic2d.h"
# include "Job_Nonadiabatic3d.h"
# include "Job_Nonadiabatic4d.h"
# include "Job_Nonadiabatic5d.h"

using namespace QTR_NS;

const char Job::NONADIABATIC1D[] = "nonadiabatic1d";
const char Job::NONADIABATIC2D[] = "nonadiabatic2d";
const char Job::NONADIABATIC3D[] = "nonadiabatic3d";
const char Job::NONADIABATIC4D[] = "nonadiabatic4d";
const char Job::NONADIABATIC5D[] = "nonadiabatic5d";

/* ------------------------------------------------------------------------------- */

Job *Job::getJob(class QTR *qtr) {
    
    Job *job = NULL;
    
    if (qtr->parameters->job == NONADIABATIC1D)
    {
        job = new JobNonadiabatic1d(qtr);
    }
    else if (qtr->parameters->job == NONADIABATIC2D) 
    {
        job = new JobNonadiabatic2d(qtr);
    }
    else if (qtr->parameters->job == NONADIABATIC3D) 
    {
        job = new JobNonadiabatic3d(qtr);
    }
    else if (qtr->parameters->job == NONADIABATIC4D)
    {
        job = new JobNonadiabatic4d(qtr);
    }
    else if (qtr->parameters->job == NONADIABATIC5D)
    {
        job = new JobNonadiabatic5d(qtr);
    }
    return job;
}
/* ----------------------------------------------------------------- */
