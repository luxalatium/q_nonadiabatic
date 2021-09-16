// ==============================================================================
//
//  Nonadiabatic4d.h
//  QTR
//
//  Created by Albert Lu on 3/12/21.
//  alu@tacc.utexas.edu
//
//  Last modified on 5/5/21
//
//  Note:
//
// ==============================================================================

#ifndef QTR_Nonadiabatic4d_H
#define QTR_Nonadiabatic4d_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class Nonadiabatic4d {
        
    public:
        Nonadiabatic4d(class QTR *q);
        ~Nonadiabatic4d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(unsigned int idx);
        inline unsigned int           GridToIdx(int x1, int x2, int x3, int x4);

        inline std::complex<double>   WaveA_1(double x1, double x2, double x3, double x4);
        inline std::complex<double>   WaveB_1(double x1, double x2, double x3, double x4);
        inline double                 VA_1(double x1, double x2, double x3, double x4);
        inline double                 VB_1(double x1, double x2, double x3, double x4);
        inline double                 VC_1(double x1, double x2, double x3, double x4);

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
        unsigned int    GRIDS_TOT;
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
        unsigned int    M1, M2, M3;
        unsigned int    W1, W2, W3;
        unsigned int    O1, O2;

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
        MeshIndexU       TB_A;    // Truncation boundary
        MeshIndexU       TB_B;    // Truncation boundary
        MeshIndexU       TBL_A;
        MeshIndexU       TBL_B;
        MeshIndexU       TBL_P_A;
        MeshIndexU       TBL_P_B;
        MeshIndexU       ExFF_A;
        MeshIndexU       ExFF_B;
        std::vector<bool> Check_A, Check_B;

        double          toggle_threshold;
        double          trans_x0;
        bool            isTrans;
        bool            isPrintEdge;
        bool            isPrintDensity;
        bool            isToggle;
    };
}

#endif /* QTR_Nonadiabatic4d_H */
