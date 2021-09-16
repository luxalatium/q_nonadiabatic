// ==============================================================================
//
//  Job_Nonadiabatic4d.h
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

#ifndef QTR_JOB_NONADIABATIC4D_H
#define QTR_JOB_NONADIABATIC4D_H

#include "Job.h"
#include "Nonadiabatic4d.h"

namespace QTR_NS  {

    class JobNonadiabatic4d: public Job    {
        
    public:
        JobNonadiabatic4d(class QTR *);
        ~JobNonadiabatic4d(void);
        
        void run(class QTR *);
        
    private:
        
        Log *log;
        Parameters *parameters;
        QTR *qtr;
        Nonadiabatic4d *nadb;
    };
}
#endif /* QTR_JOB_NONADIABATIC4D_H */
