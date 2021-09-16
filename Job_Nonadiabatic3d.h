// ==============================================================================
//
//  Job_Nonadiabatic3d.h
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

#ifndef QTR_JOB_NONADIABATIC3D_H
#define QTR_JOB_NONADIABATIC3D_H

#include "Job.h"
#include "Nonadiabatic3d.h"

namespace QTR_NS  {

    class JobNonadiabatic3d: public Job    {
        
    public:
        JobNonadiabatic3d(class QTR *);
        ~JobNonadiabatic3d(void);
        
        void run(class QTR *);
        
    private:
        
        Log *log;
        Parameters *parameters;
        QTR *qtr;
        Nonadiabatic3d *nadb;
    };
}
#endif /* QTR_JOB_NONADIABATIC3D_H */
