// ==============================================================================
//
//  Nonadiabatic3d.h
//  QTR
//
//  Created by Albert Lu on 1/26/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 5/5/21
//
//  Note:
//
// ==============================================================================

#ifndef QTR_Nonadiabatic3d_H
#define QTR_Nonadiabatic3d_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class Nonadiabatic3d {
        
    public:
        Nonadiabatic3d(class QTR *q);
        ~Nonadiabatic3d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(int idx);
        inline int                    GridToIdx(int x1, int x2, int x3);

        inline std::complex<double>   WaveA_1(double x1, double x2, double x3);
        inline std::complex<double>   WaveB_1(double x1, double x2, double x3);
        inline double                 VA_1(double x1, double x2, double x3);
        inline double                 VB_1(double x1, double x2, double x3);
        inline double                 VC_1(double x1, double x2, double x3);

    private:

        void            init();
        QTR             *qtr;
        Error           *err;
        Log             *log;
        Parameters      *parameters;

        // General parameters
        std::complex<double>  I;      // sqrt(-1)
        std::complex<double>  xZERO;  // complex zero
        int             DIMENSIONS;
        int             EDGE;
        int             PERIOD;
        int             SORT_PERIOD;
        int             PRINT_PERIOD;
        int             PRINT_WAVEFUNC_PERIOD;
        int             GRIDS_TOT;
        bool            QUIET;
        bool            TIMING;
        double          TIME;   
        double          PI_INV;  // 1/pi

        // Grid size
        double          kk;    // time resolution
        VectorXd        H;     // grid size 
        VectorXd        S;  

        // Domain size
        VectorXd        Box;
        VectorXi        BoxShape;
        int             M1, M2;
        int             W1, W2;
        int             O1;

        // Potential parameters  
        double          hb;
        double          m;

        // Wavefunction
        VectorXd        Wave0;
        VectorXd        A;
        VectorXd        P;

        // Truncate parameters
        bool            isEmpty;
        bool            isFullGrid; 
        double          PTolMin;
        double          TolH, TolH_A, TolH_B;
        double          TolL, TolL_A, TolL_B;
        double          TolHd, TolHd_A, TolHd_B;
        double          TolLd, TolLd_A, TolLd_B;
        double          ExReduce;
        int             ExLimit;

        // Domains
        MeshIndex       TB_A;    // Truncation boundary
        MeshIndex       TB_B;    // Truncation boundary
        MeshIndex       TBL_A;
        MeshIndex       TBL_B;
        MeshIndex       TBL_P_A;
        MeshIndex       TBL_P_B;
        MeshIndex       ExFF_A;
        MeshIndex       ExFF_B;
        std::vector<bool> Check_A, Check_B;

        double          toggle_threshold;
        double          trans_x0;
        bool            isTrans;
        bool            isPrintEdge;
        bool            isPrintDensity;
        bool            isToggle;
    };
}

#endif /* QTR_Nonadiabatic3d_H */
