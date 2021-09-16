// ==============================================================================
//
//  Job_Nonadiabatic5d.h
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

#ifndef QTR_JOB_NONADIABATIC5D_H
#define QTR_JOB_NONADIABATIC5D_H

#include "Job.h"
#include "Nonadiabatic5d.h"

namespace QTR_NS  {

    class JobNonadiabatic5d: public Job    {
        
    public:
        JobNonadiabatic5d(class QTR *);
        ~JobNonadiabatic5d(void);
        
        void run(class QTR *);
        
    private:
        
        Log *log;
        Parameters *parameters;
        QTR *qtr;
        Nonadiabatic5d *nadb;
    };
}
#endif /* QTR_JOB_NONADIABATIC5D_H */
