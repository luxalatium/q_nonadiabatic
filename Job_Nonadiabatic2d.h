// ==============================================================================
//
//  Job_Nonadiabatic2d.h
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

#ifndef QTR_JOB_NONADIABATIC2D_H
#define QTR_JOB_NONADIABATIC2D_H

#include "Job.h"
#include "Nonadiabatic2d.h"

namespace QTR_NS  {

    class JobNonadiabatic2d: public Job    {
        
    public:
        JobNonadiabatic2d(class QTR *);
        ~JobNonadiabatic2d(void);
        
        void run(class QTR *);
        
    private:
        
        Log *log;
        Parameters *parameters;
        QTR *qtr;
        Nonadiabatic2d *nadb;
    };
}
#endif /* QTR_JOB_NONADIABATIC2D_H */
