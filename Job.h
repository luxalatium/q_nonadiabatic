// ==============================================================================
//
//  Job.h
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

#ifndef QTR_JOB_H
#define QTR_JOB_H

#include "Qtr.h"

namespace QTR_NS {
    
    class Job {
        
    public:
        
        virtual ~Job() {};
        
        virtual void run(class QTR *) = 0;
        
        static Job *getJob(class QTR *);

        static const char NONADIABATIC1D[];
        static const char NONADIABATIC2D[];
        static const char NONADIABATIC3D[];
        static const char NONADIABATIC4D[];
        static const char NONADIABATIC5D[];
    };
}
#endif /* QTR_JOB_H */
