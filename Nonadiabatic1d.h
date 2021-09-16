// ==============================================================================
//
//  Nonadiabatic1d.h
//  QTR
//
//  Created by Albert Lu on 11/2/20.
//  alu@tacc.utexas.edu
//
//  Last modified on 3/20/21
//
//  Note:
//
// ==============================================================================

#ifndef QTR_Nonadiabatic1d_H
#define QTR_Nonadiabatic1d_H

#include <complex>

#include "Containers.h"
#include "Eigen.h"
#include "Pointers.h"

namespace QTR_NS {
    
    class Nonadiabatic1d {
        
    public:
        Nonadiabatic1d(class QTR *q);
        ~Nonadiabatic1d();
  
        void                          Evolve();
        VectorXi                      IdxToGrid(int idx);
        inline int                    GridToIdx(int x1);

        inline std::complex<double>   WaveA_1(double x1);
        inline std::complex<double>   WaveB_1(double x1);
        inline double                 VA_1(double x1);
        inline double                 VB_1(double x1);
        inline double                 VC_1(double x1);

        inline std::complex<double>   WaveA_2(double x1);
        inline std::complex<double>   WaveB_2(double x1);
        inline double                 VA_2(double x1);
        inline double                 VB_2(double x1);
        inline double                 VC_2(double x1);

        inline std::complex<double>   WaveA_3(double x1);
        inline std::complex<double>   WaveB_3(double x1);
        inline double                 VA_3(double x1);
        inline double                 VB_3(double x1);
        inline double                 VC_3(double x1);

        inline std::complex<double>   WaveA_4(double x1);
        inline std::complex<double>   WaveB_4(double x1);
        inline double                 VA_4(double x1);
        inline double                 VB_4(double x1);
        inline double                 VC_4(double x1);

        inline std::complex<double>   WaveA_5(double x1);
        inline std::complex<double>   WaveB_5(double x1);
        inline double                 VA_5(double x1);
        inline double                 VB_5(double x1);
        inline double                 VC_5(double x1);

        inline std::complex<double>   WaveA_6(double x1);
        inline std::complex<double>   WaveB_6(double x1);
        inline double                 VA_6(double x1);
        inline double                 VB_6(double x1);
        inline double                 VC_6(double x1);

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

        // Output
        double          trans_x0;
        bool            isTrans;
        bool            isPrintEdge;
        bool            isPrintDensity;
    };
}

#endif /* QTR_Nonadiabatic1d_H */
