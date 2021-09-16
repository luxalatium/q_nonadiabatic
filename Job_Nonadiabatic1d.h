// ==============================================================================
//
//  Job_Nonadiabatic1d.h
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

#ifndef QTR_JOB_NONADIABATIC1D_H
#define QTR_JOB_NONADIABATIC1D_H

#include "Job.h"
#include "Nonadiabatic1d.h"

namespace QTR_NS  {

    class JobNonadiabatic1d: public Job    {
        
    public:
        JobNonadiabatic1d(class QTR *);
        ~JobNonadiabatic1d(void);
        
        void run(class QTR *);
        
    private:
        
        Log *log;
        Parameters *parameters;
        QTR *qtr;
        Nonadiabatic1d *nadb;
    };
}
#endif /* QTR_JOB_NONADIABATIC1D_H */
