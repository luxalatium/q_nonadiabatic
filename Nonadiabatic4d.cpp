// ==============================================================================
//
//  Nonadiabatic4d.cpp
//  QTR
//
//  Created by Albert Lu on 3/12/21
//  alu@tacc.utexas.edu
//
//  Last modified on 3/29/21
//
//  Note:
//
// ==============================================================================

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <new>
#include <omp.h>
#include <parallel/algorithm>
#include <vector>

#include "Constants.h"
#include "Containers.h"
#include "Error.h"
#include "Log.h"
#include "Parameters.h"
#include "Nonadiabatic4d.h"

using namespace QTR_NS;
using std::vector;
using std::complex;
using std::real;
using std::imag;
using std::max;
using std::min;
using std::nothrow;

#define BIG_NUMBER 2147483647

/* ------------------------------------------------------------------------------- */

// DEFINE POTENTIAL

#if defined NADBPOT_TN1

#define WAVEFUNCTION_A(x1,x2,x3,x4) WaveA_1(x1,x2,x3,x4)
#define WAVEFUNCTION_B(x1,x2,x3,x4) WaveB_1(x1,x2,x3,x4)
#define POTENTIAL_A(x1,x2,x3,x4) VA_1(x1,x2,x3,x4)
#define POTENTIAL_B(x1,x2,x3,x4) VB_1(x1,x2,x3,x4)
#define POTENTIAL_C(x1,x2,x3,x4) VC_1(x1,x2,x3,x4)
#define POTNAME "Tannor1"

#else

#define WAVEFUNCTION_A(x1,x2,x3,x4) WaveA_1(x1,x2,x3,x4)
#define WAVEFUNCTION_B(x1,x2,x3,x4) WaveB_1(x1,x2,x3,x4)
#define POTENTIAL_A(x1,x2,x3,x4) VA_1(x1,x2,x3,x4)
#define POTENTIAL_B(x1,x2,x3,x4) VB_1(x1,x2,x3,x4)
#define POTENTIAL_C(x1,x2,x3,x4) VC_1(x1,x2,x3,x4)
#define POTNAME "Tannor1"

#endif
/* ------------------------------------------------------------------------------- */

Nonadiabatic4d::Nonadiabatic4d(class QTR *q)
{
    qtr = q;
    err = qtr->error;
    log = qtr->log;
    parameters = qtr->parameters;
    init();
} 
/* ------------------------------------------------------------------------------- */

Nonadiabatic4d::~Nonadiabatic4d()
{     
    return;
}
/* ------------------------------------------------------------------------------- */

void Nonadiabatic4d::init()
{
    log->log("\n\n[Nonadiabatic4d] INIT starts ...\n");
    log->log("\n\n[Nonadiabatic4d] Potential type: %s\n", POTNAME);

    // General parameters
    I = {0,1}; // sqrt(-1)
    PI_INV = 1.0 / PI; // 1/PI
    xZERO = {0,0}; // complex zero
    DIMENSIONS = parameters->scxd_dimensions;
    EDGE = parameters->scxd_edge;
    PERIOD = parameters->scxd_period;
    SORT_PERIOD = parameters->scxd_sortperiod;
    PRINT_PERIOD = parameters->scxd_printperiod;
    TIME = parameters->scxd_Tf;
    QUIET = parameters->quiet;
    TIMING = parameters->timing;
    trans_x0 = parameters->scxd_trans_x0;
    isToggle = parameters->scxd_isToggle;
    isTrans = parameters->scxd_isTrans;
    isPrintEdge = parameters->scxd_isPrintEdge;
    isPrintDensity = parameters->scxd_isPrintDensity;

    log->log("[Nonadiabatic4d] DIMENSIONS: %d\n", DIMENSIONS);
    log->log("[Nonadiabatic4d] EDGE: %d\n", EDGE);

    // Grid size
    H.resize(DIMENSIONS);
    S.resize(DIMENSIONS);  
    kk = parameters->scxd_k;

    H[0] = parameters->scxd_h1;
    H[1] = parameters->scxd_h2;
    H[2] = parameters->scxd_h3;
    H[3] = parameters->scxd_h4;

    for (int i = 0; i < DIMENSIONS; i ++)  {
        S[i] = kk / (H[i] * H[i]);
    }

    // Domain size and # grids
    Box.resize(DIMENSIONS * 2);
    Box[0] = parameters->scxd_xi1;
    Box[1] = parameters->scxd_xf1;
    Box[2] = parameters->scxd_xi2;
    Box[3] = parameters->scxd_xf2;
    Box[4] = parameters->scxd_xi3;
    Box[5] = parameters->scxd_xf3;
    Box[6] = parameters->scxd_xi4;
    Box[7] = parameters->scxd_xf4;

    BoxShape.resize(DIMENSIONS);
    GRIDS_TOT = 1;
    log->log("[Nonadiabatic4d] Number of grids = (");

    for (int i = 0; i < DIMENSIONS; i ++)  {

        BoxShape[i] = (int)std::round((Box[2 * i + 1] - Box[2 * i]) / H[i]) + 1;
        GRIDS_TOT *= BoxShape[i];

        if ( i < DIMENSIONS - 1 )
            log->log("%d, ", BoxShape[i]);
        else
            log->log("%d)\n", BoxShape[i]);
    }

    if (DIMENSIONS==4)  {
        O1 = BoxShape[0];
        O2 = BoxShape[1] * BoxShape[2] * BoxShape[3];
        M1 = BoxShape[1] * BoxShape[2] * BoxShape[3];
        M2 = BoxShape[2] * BoxShape[3];
        M3 = BoxShape[3];
        W1 = BoxShape[1] * BoxShape[2] * BoxShape[3];
        W2 = BoxShape[2] * BoxShape[3];
        W3 = BoxShape[3];
    }

    // Parameters
    hb = parameters->scxd_hb;
    m  = parameters->scxd_m;
    log->log("[Nonadiabatic4d] hb: %lf\n", hb);
    log->log("[Nonadiabatic4d] m: %lf\n", m);

    // Wavefunction parameters
    Wave0.resize(DIMENSIONS);
    Wave0[0] = parameters->scxd_x01;
    Wave0[1] = parameters->scxd_x02;
    Wave0[2] = parameters->scxd_x03;
    Wave0[3] = parameters->scxd_x04;
    log->log("[Nonadiabatic4d] Wave0[0]: %lf\n", Wave0[0]);
    log->log("[Nonadiabatic4d] Wave0[1]: %lf\n", Wave0[1]);
    log->log("[Nonadiabatic4d] Wave0[2]: %lf\n", Wave0[2]);
    log->log("[Nonadiabatic4d] Wave0[3]: %lf\n", Wave0[3]);

    A.resize(DIMENSIONS);
    A[0] = parameters->scxd_a1;
    A[1] = parameters->scxd_a2;
    A[2] = parameters->scxd_a3;
    A[3] = parameters->scxd_a4;
    log->log("[Nonadiabatic4d] A[0]: %lf\n", A[0]);
    log->log("[Nonadiabatic4d] A[1]: %lf\n", A[1]);
    log->log("[Nonadiabatic4d] A[2]: %lf\n", A[2]);
    log->log("[Nonadiabatic4d] A[3]: %lf\n", A[3]);

    P.resize(DIMENSIONS);
    P[0] = parameters->scxd_p1;
    P[1] = parameters->scxd_p2;
    P[2] = parameters->scxd_p3;
    P[3] = parameters->scxd_p4;
    log->log("[Nonadiabatic4d] P[0]: %lf\n", P[0]);
    log->log("[Nonadiabatic4d] P[1]: %lf\n", P[1]);
    log->log("[Nonadiabatic4d] P[2]: %lf\n", P[2]);
    log->log("[Nonadiabatic4d] P[3]: %lf\n", P[3]);

    // Truncate parameters
    isFullGrid = parameters->scxd_isFullGrid;
    TolH = parameters->scxd_TolH;    // Tolerance of probability density for Zero point Cutoff
    TolL = parameters->scxd_TolL;    // Tolerance of probability density for Edge point
    TolHd = parameters->scxd_TolHd;  // Tolerance of probability first diff for Zero point Cutoff
    TolLd = parameters->scxd_TolLd;  // Tolerance of probability density for Edge point
    PTolMin = parameters->scxd_PTolMin; 
    ExReduce = parameters->scxd_ExReduce; //Extrapolation reduce factor
    ExLimit = parameters->scxd_ExLimit;   //Extrapolation counts limit

    log->log("[Nonadiabatic4d] isFullGrid: %d\n", (int)isFullGrid);
    log->log("[Nonadiabatic4d] isToggle: %d\n", (int)isToggle);

    if (isToggle)  {
        toggle_threshold = parameters->scxd_toggle_threshold;
        log->log("[Nonadiabatic4d] Toggle threshold: %lf\n", toggle_threshold);
    }
    log->log("[Nonadiabatic4d] TolH: %e\n", TolH);
    log->log("[Nonadiabatic4d] TolL: %e\n", TolL);
    log->log("[Nonadiabatic4d] TolHd: %e\n", TolHd);
    log->log("[Nonadiabatic4d] TolLd: %e\n", TolLd);
    log->log("[Nonadiabatic4d] PTolMin: %e\n", PTolMin);
    log->log("[Nonadiabatic4d] ExReduce: %lf\n", ExReduce);
    log->log("[Nonadiabatic4d] ExLimit: %d\n", ExLimit);
    log->log("[Nonadiabatic4d] INIT done.\n\n");
}
/* ------------------------------------------------------------------------------- */

void Nonadiabatic4d::Evolve()
{
    #pragma omp declare reduction (merge : MeshIndexU : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))

    log->log("[Nonadiabatic4d] Evolve starts ...\n");

    // Files
    FILE *pfile;

    // Variables 
    int count;
    int n1, n2, n3, n4;
    int nx1, nx2, nx3, nx4;
    int ta_size_A, tb_size_A;
    int ta_size_B, tb_size_B;
    int x1_max_A, x1_min_A, x2_max_A, x2_min_A, x3_max_A, x3_min_A, x4_max_A, x4_min_A;
    int x1_max_B, x1_min_B, x2_max_B, x2_min_B, x3_max_B, x3_min_B, x4_max_B, x4_min_B;  
    int x1_max_C, x1_min_C, x2_max_C, x2_min_C, x3_max_C, x3_min_C, x4_max_C, x4_min_C;  
    int x1_max_T, x1_min_T, x2_max_T, x2_min_T, x3_max_T, x3_min_T, x4_max_T, x4_min_T;  
    int idx_x0;  // transmission x0
    complex<double> sum;
    double norm, norm_A, norm_B;  // normalization factor
    double pftrans;
    bool b1, b2, b3, b4, b5, b6, b7;
    bool isEmpty;

    // Timing variables
    double t_0_begin, t_1_begin, t_2_begin;    
    double t_0_end, t_1_end, t_2_end;
    double t_0_elapsed = 0.0;
    double t_1_elapsed = 0.0;
    double t_2_elapsed = 0.0;

    // Core computation time (RK4, normalization, initialization, etc)
    double t_full = 0.0;
    double t_truncate = 0.0;

    // Overhead time (truncate)
    double t_overhead = 0.0;

    // Constants
    complex<double> Ikh2mh0sq = I * kk * hb / (2.0 * m * H[0] * H[0]);
    complex<double> Ikh2mh1sq = I * kk * hb / (2.0 * m * H[1] * H[1]);
    complex<double> Ikh2mh2sq = I * kk * hb / (2.0 * m * H[2] * H[2]);
    complex<double> Ikh2mh3sq = I * kk * hb / (2.0 * m * H[3] * H[3]);
    complex<double> Ik2h = I * kk / hb;
    double TolHd_A_sq, TolHd_B_sq;
    double TolLd_A_sq, TolLd_B_sq;

    // temporary index container
    MeshIndexU tmpVec; 

    // Boundary layer container for extrapolation loop
    MeshIndexU ExBD, ExBD_A, ExBD_B;     
 
    //  1d Grid vector and indices
    VectorXi grid;
    int g1,g2,g3,g4,g5;
    double xx1,xx2,xx3,xx4;

    complex<double> f0_A, f0_B;

    complex<double> f1p1_A, f1m1_A;
    complex<double> f1p1_B, f1m1_B;
    complex<double> f2p1_A, f2m1_A;
    complex<double> f2p1_B, f2m1_B;
    complex<double> f3p1_A, f3m1_A;
    complex<double> f3p1_B, f3m1_B;
    complex<double> f4p1_A, f4m1_A;
    complex<double> f4p1_B, f4m1_B;

    complex<double> kk0_A, kk0_B;
    complex<double> kk1p1_A, kk1m1_A;
    complex<double> kk1p1_B, kk1m1_B;
    complex<double> kk2p1_A, kk2m1_A;
    complex<double> kk2p1_B, kk2m1_B;
    complex<double> kk3p1_A, kk3m1_A;
    complex<double> kk3p1_B, kk3m1_B;
    complex<double> kk4p1_A, kk4m1_A;
    complex<double> kk4p1_B, kk4m1_B;
    
    complex<double> c_f0_a, c_f0_b;
    double pot_a, pot_b, pot_c;

    // Vector iterater
    vector<unsigned int>::iterator it;

    // Extrapolation 
    int min_dir;
    int Excount_A, Excount_B;
    bool isExtrapolate_A, isExtrapolate_B; 
    bool isFirstExtrp_A, isFirstExtrp_B;
    complex<double> val, val_min;
    double val_min_abs;
    vector<complex<double>> ExTBL;

    // Neighborlist
    int nneigh = 0;
    vector<vector<int>> neighlist;
    vector<int> neighs(DIMENSIONS);

    log->log("[Nonadiabatic4d] Initializing containers ...\n");

    // Initialize containers

    t_0_begin = omp_get_wtime();

    bool **TAMask_A;
    bool **TAMask_B;

    if ( !isFullGrid )  {

        TAMask_A = new bool*[O1];
        TAMask_B = new bool*[O1];

        for (unsigned int i = 0; i < O1; i++)  {
            TAMask_A[i] = new(nothrow) bool[O2];
            if (!TAMask_A[i])  {
                log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "TAMask_A");
                std::exit(EXIT_FAILURE);
            }
        }

        for (unsigned int i = 0; i < O1; i++)  {
            TAMask_B[i] = new(nothrow) bool[O2];
            if (!TAMask_B[i])  {
                log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "TAMask_B");
                std::exit(EXIT_FAILURE);
            }
        }
    }

    double **PF_A = new double*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        PF_A[i] = new(nothrow) double[O2];
        if (!PF_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "PF_A");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **F_A = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        F_A[i] = new(nothrow) complex<double>[O2];
        if (!F_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "F_A");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **FF_A = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        FF_A[i] = new(nothrow) complex<double>[O2];
        if (!FF_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "FF_A");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK1_A = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK1_A[i] = new(nothrow) complex<double>[O2];
        if (!KK1_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK1_A");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK2_A = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK2_A[i] = new(nothrow) complex<double>[O2];
        if (!KK2_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK2_A");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK3_A = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK3_A[i] = new(nothrow) complex<double>[O2];
        if (!KK3_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK3_A");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK4_A = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK4_A[i] = new(nothrow) complex<double>[O2];
        if (!KK4_A[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK4_A");
            std::exit(EXIT_FAILURE);
        }
    }

    double **PF_B = new double*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        PF_B[i] = new(nothrow) double[O2];
        if (!PF_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "PF_B");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **F_B = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        F_B[i] = new(nothrow) complex<double>[O2];
        if (!F_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "F_B");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **FF_B = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        FF_B[i] = new(nothrow) complex<double>[O2];
        if (!FF_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "FF_B");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK1_B = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK1_B[i] = new(nothrow) complex<double>[O2];
        if (!KK1_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK1_B");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK2_B = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK2_B[i] = new(nothrow) complex<double>[O2];
        if (!KK2_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK2_B");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK3_B = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK3_B[i] = new(nothrow) complex<double>[O2];
        if (!KK3_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK3_B");
            std::exit(EXIT_FAILURE);
        }
    }
    complex<double> **KK4_B = new complex<double>*[O1];
    for (unsigned int i = 0; i < O1; i++)  {
        KK4_B[i] = new(nothrow) complex<double>[O2];
        if (!KK4_B[i])  {
            log->log("[Nonadiabatic4d] Memory allocation %s failed\n", "KK4_B");
            std::exit(EXIT_FAILURE);
        }
    }

    #pragma omp parallel for schedule(runtime)
    for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
        for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
            for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {
                for (int i4 = 0; i4 < BoxShape[3]; i4 ++)  {
                    PF_A[i1][i2*W2+i3*W3+i4] = 0.0;
                    PF_B[i1][i2*W2+i3*W3+i4] = 0.0;
                    F_A[i1][i2*W2+i3*W3+i4] = xZERO;
                    F_B[i1][i2*W2+i3*W3+i4] = xZERO;
                    FF_A[i1][i2*W2+i3*W3+i4] = xZERO;
                    FF_B[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK1_A[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK1_B[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK2_A[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK2_B[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK3_A[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK3_B[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK4_A[i1][i2*W2+i3*W3+i4] = xZERO;
                    KK4_B[i1][i2*W2+i3*W3+i4] = xZERO;
                }
            }
        }
    }

    if ( !isFullGrid )  {

        t_1_begin = omp_get_wtime();
        
        #pragma omp parallel for schedule(runtime)
        for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
            for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {
                    for (int i4 = 0; i4 < BoxShape[3]; i4 ++)  {
                        TAMask_A[i1][i2*W2+i3*W3+i4] = 0;
                        TAMask_B[i1][i2*W2+i3*W3+i4] = 0;
                    }
                }
            }
        }

        nneigh = 0;

        for (int d = 1; d <= 4; d ++)  {
            for (n1 = -d; n1 <= d; n1 ++)  {
                for (n2 = -(d-abs(n1)); n2 <= d-abs(n1); n2 ++)  {
                    for (n3 = -(d-abs(n1)-abs(n2)); n3 <= d-abs(n1)-abs(n2); n3 ++)  {

                        n4 = d - abs(n1) - abs(n2) - abs(n3);

                        if (n4 != 0)  {
                            neighs = {n1,n2,n3,n4};
                            neighlist.push_back(neighs);
                            neighs = {n1,n2,n3,-n4};
                            neighlist.push_back(neighs);
                            nneigh += 2;
                        }
                        else  {
                            neighs = {n1,n2,n3,0};
                            neighlist.push_back(neighs);
                            nneigh += 1;
                        }
                    }
                }
            }
        }

        log->log("[Nonadiabatic4d] nneigh = %d\n", nneigh); 
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_overhead += t_1_elapsed;
    }
    t_0_end = omp_get_wtime();
    t_0_elapsed = t_0_end - t_0_begin;
    t_full += t_0_elapsed;

    if ( !isFullGrid )
        t_truncate += t_0_elapsed - t_1_elapsed; // subtract overhead

    if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing containers) = %.4e sec\n\n", t_0_elapsed); 

    // .........................................................................................

    log->log("[Nonadiabatic4d] Initializing wavefunction ...\n");  

    t_1_begin = omp_get_wtime();

    // Initialize wavefunction

    /* ---------------- Surface A ---------------- */

    #pragma omp parallel for schedule(runtime)
    for (int i1 = EDGE; i1 <  BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 <  BoxShape[1]-EDGE; i2 ++)  {
            for (int i3 = EDGE; i3 <  BoxShape[2]-EDGE; i3 ++)  {
                for (int i4 = EDGE; i4 <  BoxShape[3]-EDGE; i4 ++)  {
                    F_A[i1][i2*W2+i3*W3+i4] = WAVEFUNCTION_A(Box[0]+i1*H[0], Box[2]+i2*H[1], Box[4]+i3*H[2], Box[6]+i4*H[3]);
                }
            }
        }
    }
    norm_A = 0.0;

    #pragma omp parallel for reduction (+:norm_A) schedule(runtime)
    for (int i1 = EDGE; i1 <  BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 <  BoxShape[1]-EDGE; i2 ++)  {
            for (int i3 = EDGE; i3 <  BoxShape[2]-EDGE; i3 ++)  {
                for (int i4 = EDGE; i4 <  BoxShape[3]-EDGE; i4 ++)  {
                    norm_A += std::abs(F_A[i1][i2*W2+i3*W3+i4] * std::conj(F_A[i1][i2*W2+i3*W3+i4]));
                }
            }
        }
    }
    norm_A *= H[0] * H[1] * H[2] * H[3];
    log->log("[Nonadiabatic4d] Normalization factor (A) = %.16e\n",norm_A);

    /* ---------------- Surface B ---------------- */
    
    #pragma omp parallel for schedule(runtime)
    for (int i1 = EDGE; i1 <  BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 <  BoxShape[1]-EDGE; i2 ++)  {
            for (int i3 = EDGE; i3 <  BoxShape[2]-EDGE; i3 ++)  {
                for (int i4 = EDGE; i4 <  BoxShape[3]-EDGE; i4 ++)  {
                    F_B[i1][i2*W2+i3*W3+i4] = WAVEFUNCTION_B(Box[0]+i1*H[0], Box[2]+i2*H[1], Box[4]+i3*H[2], Box[6]+i4*H[3]);
                }
            }
        }
    }
    norm_B = 0.0;

    #pragma omp parallel for reduction (+:norm_B) schedule(runtime)
    for (int i1 = EDGE; i1 <  BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 <  BoxShape[1]-EDGE; i2 ++)  {
            for (int i3 = EDGE; i3 <  BoxShape[2]-EDGE; i3 ++)  {
                for (int i4 = EDGE; i4 <  BoxShape[3]-EDGE; i4 ++)  {
                    norm_B += std::abs(F_B[i1][i2*W2+i3*W3+i4] * std::conj(F_B[i1][i2*W2+i3*W3+i4]));
                }
            }
        }
    }
    norm_B *= H[0] * H[1] * H[2] * H[3];
    log->log("[Nonadiabatic4d] Normalization factor (B) = %.16e\n",norm_B);

    /* ------------------------------- */

    norm = 1.0 / sqrt(norm_A + norm_B);
    log->log("[Nonadiabatic4d] Normalization factor 1/sqrt(A+B) = %.16e\n",norm);

    #pragma omp parallel for schedule(runtime)
    for (int i1 = EDGE; i1 <  BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 <  BoxShape[1]-EDGE; i2 ++)  {
            for (int i3 = EDGE; i3 <  BoxShape[2]-EDGE; i3 ++)  {
                for (int i4 = EDGE; i4 <  BoxShape[3]-EDGE; i4 ++)  {
                    F_A[i1][i2*W2+i3*W3+i4] = norm * F_A[i1][i2*W2+i3*W3+i4];
                    PF_A[i1][i2*W2+i3*W3+i4] = std::abs(F_A[i1][i2*W2+i3*W3+i4] * std::conj(F_A[i1][i2*W2+i3*W3+i4]));
                }
            }
        }
    }
    #pragma omp parallel for schedule(runtime)
    for (int i1 = EDGE; i1 <  BoxShape[0]-EDGE; i1 ++)  {
        for (int i2 = EDGE; i2 <  BoxShape[1]-EDGE; i2 ++)  {
            for (int i3 = EDGE; i3 <  BoxShape[2]-EDGE; i3 ++)  {
                for (int i4 = EDGE; i4 <  BoxShape[3]-EDGE; i4 ++)  {
                    F_B[i1][i2*W2+i3*W3+i4] = norm * F_B[i1][i2*W2+i3*W3+i4];
                    PF_B[i1][i2*W2+i3*W3+i4] = std::abs(F_B[i1][i2*W2+i3*W3+i4] * std::conj(F_B[i1][i2*W2+i3*W3+i4]));
                }
            }
        }
    }
    t_1_end = omp_get_wtime();
    t_1_elapsed = t_1_end - t_1_begin;
    t_full += t_1_elapsed;
    t_truncate += t_1_elapsed;
    if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing wavefunction) = %.4e sec\n\n", t_1_elapsed); 

    // .........................................................................................

    // Initial truncation & edge point check

    if ( !isFullGrid )
    {
        t_2_begin = omp_get_wtime();

        log->log("[Nonadiabatic4d] Initial truncation ...\n");

        // Population-weighted TG criteria

        TolH_A = TolH * max(norm_A/(norm_A+norm_B), PTolMin);
        TolL_A = TolL * max(norm_A/(norm_A+norm_B), PTolMin);
        TolH_B = TolH * max(norm_B/(norm_A+norm_B), PTolMin);
        TolL_B = TolL * max(norm_B/(norm_A+norm_B), PTolMin);

        TolHd_A = TolHd;
        TolLd_A = TolLd;
        TolHd_B = TolHd;
        TolLd_B = TolLd;

        TolHd_A_sq = TolHd_A * TolHd_A;
        TolLd_A_sq = TolLd_A * TolLd_A;
        TolHd_B_sq = TolHd_B * TolHd_B;
        TolLd_B_sq = TolLd_B * TolLd_B;

        x1_min_A = BIG_NUMBER;
        x2_min_A = BIG_NUMBER;
        x3_min_A = BIG_NUMBER;
        x4_min_A = BIG_NUMBER;
        x1_max_A = -BIG_NUMBER;
        x2_max_A = -BIG_NUMBER;
        x3_max_A = -BIG_NUMBER;
        x4_max_A = -BIG_NUMBER;

        x1_min_B = BIG_NUMBER;
        x2_min_B = BIG_NUMBER;
        x3_min_B = BIG_NUMBER;
        x4_min_B = BIG_NUMBER;
        x1_max_B = -BIG_NUMBER;
        x2_max_B = -BIG_NUMBER;
        x3_max_B = -BIG_NUMBER;
        x4_max_B = -BIG_NUMBER;

        /* ---------------- Surface A ---------------- */

        // Truncation

        ta_size_A = 0;

        t_1_begin = omp_get_wtime();

        #pragma omp parallel for reduction(+: ta_size_A) private(b1,b2,b3,b4,nx1,nx2,nx3,\
                                                                f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A) schedule(runtime)
        for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
            for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                    for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                        nx1 = int(TAMask_A[i1+1][i2*W2+i3*W3+i4]) + int(TAMask_A[i1-1][i2*W2+i3*W3+i4]);
                        nx2 = int(TAMask_A[i1][(i2+1)*W2+i3*W3+i4]) + int(TAMask_A[i1][(i2-1)*W2+i3*W3+i4]);
                        nx3 = int(TAMask_A[i1][i2*W2+(i3+1)*W3+i4]) + int(TAMask_A[i1][i2*W2+(i3-1)*W3+i4]);
                        nx4 = int(TAMask_A[i1][i2*W2+i3*W3+(i4+1)]) + int(TAMask_A[i1][i2*W2+i3*W3+(i4-1)]);

                        f1p1_A = (TAMask_A[i1+1][i2*W2+i3*W3+i4]) ? F_A[i1+1][i2*W2+i3*W3+i4] : F_A[i1][i2*W2+i3*W3+i4];
                        f1m1_A = (TAMask_A[i1-1][i2*W2+i3*W3+i4]) ? F_A[i1-1][i2*W2+i3*W3+i4] : F_A[i1][i2*W2+i3*W3+i4];
                        f2p1_A = (TAMask_A[i1][(i2+1)*W2+i3*W3+i4]) ? F_A[i1][(i2+1)*W2+i3*W3+i4] : F_A[i1][i2*W2+i3*W3+i4];
                        f2m1_A = (TAMask_A[i1][(i2-1)*W2+i3*W3+i4]) ? F_A[i1][(i2-1)*W2+i3*W3+i4] : F_A[i1][i2*W2+i3*W3+i4];
                        f3p1_A = (TAMask_A[i1][i2*W2+(i3+1)*W3+i4]) ? F_A[i1][i2*W2+(i3+1)*W3+i4] : F_A[i1][i2*W2+i3*W3+i4];
                        f3m1_A = (TAMask_A[i1][i2*W2+(i3-1)*W3+i4]) ? F_A[i1][i2*W2+(i3-1)*W3+i4] : F_A[i1][i2*W2+i3*W3+i4];
                        f4p1_A = (TAMask_A[i1][i2*W2+i3*W3+(i4+1)]) ? F_A[i1][i2*W2+i3*W3+(i4+1)] : F_A[i1][i2*W2+i3*W3+i4];
                        f4m1_A = (TAMask_A[i1][i2*W2+i3*W3+(i4-1)]) ? F_A[i1][i2*W2+i3*W3+(i4-1)] : F_A[i1][i2*W2+i3*W3+i4];

                        b1 = PF_A[i1][i2*W2+i3*W3+i4] < TolH_A;
                        b2 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_A - f1m1_A)/(nx1*H[0]),2)) + \
                             ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_A - f2m1_A)/(nx2*H[1]),2)) + \
                             ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_A - f3m1_A)/(nx3*H[2]),2)) + \
                             ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_A - f4m1_A)/(nx4*H[3]),2)) <  TolHd_A_sq;
                        
                        if (b1 && b2)
                            PF_A[i1][i2*W2+i3*W3+i4] = 0.0;
                        else {
                            TAMask_A[i1][i2*W2+i3*W3+i4] = 1;
                            ta_size_A += 1;
                        }
                    }
                }
            }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-1) = %.4e sec\n\n", t_1_elapsed); 

        // TA box and TB

        if (ta_size_A == 0)  {
            tb_size_A = 0;
            log->log("[Nonadiabatic4d] TA_A is empty\n");
        }
        else  {
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min_A,x2_min_A,x3_min_A,x4_min_A) reduction(max: x1_max_A,x2_max_A,x3_max_A,x4_max_A) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                        for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])  {
                                if (i1 < x1_min_A)  x1_min_A = i1;
                                if (i1 > x1_max_A)  x1_max_A = i1;
                                if (i2 < x2_min_A)  x2_min_A = i2;
                                if (i2 > x2_max_A)  x2_max_A = i2;
                                if (i3 < x3_min_A)  x3_min_A = i3;
                                if (i3 > x3_max_A)  x3_max_A = i3;
                                if (i4 < x4_min_A)  x4_min_A = i4;
                                if (i4 > x4_max_A)  x4_max_A = i4;
                            }
                            else  {
                                F_A[i1][i2*W2+i3*W3+i4] = xZERO;
                            }
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-2) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // TB
            #pragma omp parallel for reduction(merge: tmpVec) schedule(runtime)
            for (int i1 = x1_min_A; i1 <= x1_max_A; i1 ++)  {
                for (int i2 = x2_min_A; i2 <= x2_max_A; i2 ++)  {
                    for (int i3 = x3_min_A; i3 <= x3_max_A; i3 ++)  {
                        for (int i4 = x4_min_A; i4 <= x4_max_A; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4] && ( \
                                !TAMask_A[i1+1][i2*W2+i3*W3+i4] || !TAMask_A[i1-1][i2*W2+i3*W3+i4] || \
                                !TAMask_A[i1][(i2+1)*W2+i3*W3+i4] || !TAMask_A[i1][(i2-1)*W2+i3*W3+i4] || \
                                !TAMask_A[i1][i2*W2+(i3+1)*W3+i4] || !TAMask_A[i1][i2*W2+(i3-1)*W3+i4] || \
                                !TAMask_A[i1][i2*W2+i3*W3+(i4+1)] || !TAMask_A[i1][i2*W2+i3*W3+(i4-1)])) 
                                tmpVec.push_back(GridToIdx(i1,i2,i3,i4));
                        }
                    }
                }
            }
            tmpVec.swap(TB_A);
            tmpVec.clear();
            tb_size_A = TB_A.size();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-3) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // TA expansion

            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < TB_A.size(); i++)  {
                g1 = (int)(TB_A[i] / M1);
                g2 = (int)((TB_A[i]%M1) / M2);
                g3 = (int)((TB_A[i]%M2) / M3);
                g4 = (int)(TB_A[i] % M3);

                if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_A[g1+1][g2*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                if ( g1-1 > EDGE && !TAMask_A[g1-1][g2*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));

                if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_A[g1][(g2+1)*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                if ( g2-1 > EDGE && !TAMask_A[g1][(g2-1)*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));

                if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_A[g1][g2*W2+(g3+1)*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                if ( g3-1 > EDGE && !TAMask_A[g1][g2*W2+(g3-1)*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));

                if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_A[g1][g2*W2+g3*W3+(g4+1)])
                    tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                if ( g4-1 > EDGE && !TAMask_A[g1][g2*W2+g3*W3+(g4-1)])
                    tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-4) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // Find unique elements
            __gnu_parallel::sort(tmpVec.begin(),tmpVec.end());
            it = std::unique (tmpVec.begin(), tmpVec.end()); 
            tmpVec.resize(std::distance(tmpVec.begin(),it));

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-5) = %.4e sec\n\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // Update TA box

            #pragma omp parallel for reduction(min: x1_min_A,x2_min_A,x3_min_A,x4_min_A) \
                                     reduction(max: x1_max_A,x2_max_A,x3_max_A,x4_max_A) \
                                     private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < tmpVec.size(); i ++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)((tmpVec[i]%M1) / M2);
                g3 = (int)((tmpVec[i]%M2) / M3);
                g4 = (int)(tmpVec[i] % M3);

                if (!TAMask_A[g1][g2*W2+g3*W3+g4])
                    TAMask_A[g1][g2*W2+g3*W3+g4] = 1;

                x1_min_A = (g1 < x1_min_A) ? g1 : x1_min_A;
                x2_min_A = (g2 < x2_min_A) ? g2 : x2_min_A;
                x3_min_A = (g3 < x3_min_A) ? g3 : x3_min_A;
                x4_min_A = (g4 < x4_min_A) ? g4 : x4_min_A;
                x1_max_A = (g1 > x1_max_A) ? g1 : x1_max_A;
                x2_max_A = (g2 > x2_max_A) ? g2 : x2_max_A;
                x3_max_A = (g3 > x3_max_A) ? g3 : x3_max_A;
                x4_max_A = (g4 > x4_max_A) ? g4 : x4_max_A;
            }
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-6) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // Update ta_size

            ta_size_A = 0;

            #pragma omp parallel for reduction(+: ta_size_A) schedule(runtime)
            for (int i1 = x1_min_A; i1 <= x1_max_A; i1 ++)  {
                for (int i2 = x2_min_A; i2 <= x2_max_A; i2 ++)  {
                    for (int i3 = x3_min_A; i3 <= x3_max_A; i3 ++)  {
                        for (int i4 = x4_min_A; i4 <= x4_max_A; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])
                                ta_size_A += 1;    
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation A-7) = %.4e sec\n\n", t_1_elapsed); 
        }
        log->log("[Nonadiabatic4d] TA size (A) = %d, TB size (A) = %d\n", ta_size_A, tb_size_A);

        if (ta_size_A != 0)
            log->log("[Nonadiabatic4d] TA Range (A) [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min_A, x1_max_A, x2_min_A, x2_max_A, x3_min_A, x3_max_A, x4_min_A, x4_max_A);

        /* ---------------- Surface B ---------------- */

        // Truncation

        ta_size_B = 0;

        t_1_begin = omp_get_wtime();

        // Truncation and TA
        #pragma omp parallel for reduction(+: ta_size_B) private(b1,b2,b3,b4,nx1,nx2,nx3,\
                                                                f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B) schedule(runtime)
        for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
            for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                    for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                        nx1 = int(TAMask_B[i1+1][i2*W2+i3*W3+i4]) + int(TAMask_B[i1-1][i2*W2+i3*W3+i4]);
                        nx2 = int(TAMask_B[i1][(i2+1)*W2+i3*W3+i4]) + int(TAMask_B[i1][(i2-1)*W2+i3*W3+i4]);
                        nx3 = int(TAMask_B[i1][i2*W2+(i3+1)*W3+i4]) + int(TAMask_B[i1][i2*W2+(i3-1)*W3+i4]);
                        nx4 = int(TAMask_B[i1][i2*W2+i3*W3+(i4+1)]) + int(TAMask_B[i1][i2*W2+i3*W3+(i4-1)]);

                        f1p1_B = (TAMask_B[i1+1][i2*W2+i3*W3+i4]) ? F_B[i1+1][i2*W2+i3*W3+i4] : F_B[i1][i2*W2+i3*W3+i4];
                        f1m1_B = (TAMask_B[i1-1][i2*W2+i3*W3+i4]) ? F_B[i1-1][i2*W2+i3*W3+i4] : F_B[i1][i2*W2+i3*W3+i4];
                        f2p1_B = (TAMask_B[i1][(i2+1)*W2+i3*W3+i4]) ? F_B[i1][(i2+1)*W2+i3*W3+i4] : F_B[i1][i2*W2+i3*W3+i4];
                        f2m1_B = (TAMask_B[i1][(i2-1)*W2+i3*W3+i4]) ? F_B[i1][(i2-1)*W2+i3*W3+i4] : F_B[i1][i2*W2+i3*W3+i4];
                        f3p1_B = (TAMask_B[i1][i2*W2+(i3+1)*W3+i4]) ? F_B[i1][i2*W2+(i3+1)*W3+i4] : F_B[i1][i2*W2+i3*W3+i4];
                        f3m1_B = (TAMask_B[i1][i2*W2+(i3-1)*W3+i4]) ? F_B[i1][i2*W2+(i3-1)*W3+i4] : F_B[i1][i2*W2+i3*W3+i4];
                        f4p1_B = (TAMask_B[i1][i2*W2+i3*W3+(i4+1)]) ? F_B[i1][i2*W2+i3*W3+(i4+1)] : F_B[i1][i2*W2+i3*W3+i4];
                        f4m1_B = (TAMask_B[i1][i2*W2+i3*W3+(i4-1)]) ? F_B[i1][i2*W2+i3*W3+(i4-1)] : F_B[i1][i2*W2+i3*W3+i4];

                        b1 = PF_B[i1][i2*W2+i3*W3+i4] < TolH_B;
                        b2 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_B - f1m1_B)/(nx1*H[0]),2)) + \
                             ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_B - f2m1_B)/(nx2*H[1]),2)) + \
                             ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_B - f3m1_B)/(nx3*H[2]),2)) + \
                             ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_B - f4m1_B)/(nx4*H[3]),2)) <  TolHd_B_sq;
                        
                        if (b1 && b2)
                            PF_B[i1][i2*W2+i3*W3+i4] = 0.0;
                        else {
                            TAMask_B[i1][i2*W2+i3*W3+i4] = 1;
                            ta_size_B += 1;
                        }
                    }
                }
            }
        }
        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-1) = %.4e sec\n\n", t_1_elapsed); 

        // TA box and TB

        if (ta_size_B == 0)  {
            tb_size_B = 0;
            log->log("[Nonadiabatic4d] TA_B is empty\n");
        }
        else  {
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min_B,x2_min_B,x3_min_B,x4_min_B) reduction(max: x1_max_B,x2_max_B,x3_max_B,x4_max_B) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                        for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                if (i1 < x1_min_B)  x1_min_B = i1;
                                if (i1 > x1_max_B)  x1_max_B = i1;
                                if (i2 < x2_min_B)  x2_min_B = i2;
                                if (i2 > x2_max_B)  x2_max_B = i2;
                                if (i3 < x3_min_B)  x3_min_B = i3;
                                if (i3 > x3_max_B)  x3_max_B = i3;
                                if (i4 < x4_min_B)  x4_min_B = i4;
                                if (i4 > x4_max_B)  x4_max_B = i4;
                            }
                            else  {
                                F_B[i1][i2*W2+i3*W3+i4] = xZERO;
                            }
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-2) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // TB
            #pragma omp parallel for reduction(merge: tmpVec) schedule(runtime)
            for (int i1 = x1_min_B; i1 <= x1_max_B; i1 ++)  {
                for (int i2 = x2_min_B; i2 <= x2_max_B; i2 ++)  {
                    for (int i3 = x3_min_B; i3 <= x3_max_B; i3 ++)  {
                        for (int i4 = x4_min_B; i4 <= x4_max_B; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4] && ( \
                                !TAMask_B[i1+1][i2*W2+i3*W3+i4] || !TAMask_B[i1-1][i2*W2+i3*W3+i4] || \
                                !TAMask_B[i1][(i2+1)*W2+i3*W3+i4] || !TAMask_B[i1][(i2-1)*W2+i3*W3+i4] || \
                                !TAMask_B[i1][i2*W2+(i3+1)*W3+i4] || !TAMask_B[i1][i2*W2+(i3-1)*W3+i4] || \
                                !TAMask_B[i1][i2*W2+i3*W3+(i4+1)] || !TAMask_B[i1][i2*W2+i3*W3+(i4-1)])) 
                                tmpVec.push_back(GridToIdx(i1,i2,i3,i4));
                        }
                    }
                }
            }
            tmpVec.swap(TB_B);
            tmpVec.clear();
            tb_size_B = TB_B.size();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-3) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // TA expansion

            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < TB_B.size(); i++)  {
                g1 = (int)(TB_B[i] / M1);
                g2 = (int)((TB_B[i]%M1) / M2);
                g3 = (int)((TB_B[i]%M2) / M3);
                g4 = (int)(TB_B[i] % M3);

                if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_B[g1+1][g2*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                if ( g1-1 > EDGE && !TAMask_B[g1-1][g2*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));

                if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_B[g1][(g2+1)*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                if ( g2-1 > EDGE && !TAMask_B[g1][(g2-1)*W2+g3*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));

                if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_B[g1][g2*W2+(g3+1)*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                if ( g3-1 > EDGE && !TAMask_B[g1][g2*W2+(g3-1)*W3+g4])
                    tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));

                if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_B[g1][g2*W2+g3*W3+(g4+1)])
                    tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                if ( g4-1 > EDGE && !TAMask_B[g1][g2*W2+g3*W3+(g4-1)])
                    tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-4) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // Find unique elements
            __gnu_parallel::sort(tmpVec.begin(),tmpVec.end());
            it = std::unique (tmpVec.begin(), tmpVec.end()); 
            tmpVec.resize(std::distance(tmpVec.begin(),it));

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-5) = %.4e sec\n\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            // Update TA box
            #pragma omp parallel for reduction(min: x1_min_B,x2_min_B,x3_min_B,x4_min_B) \
                                     reduction(max: x1_max_B,x2_max_B,x3_max_B,x4_max_B) \
                                     private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < tmpVec.size(); i ++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)((tmpVec[i]%M1) / M2);
                g3 = (int)((tmpVec[i]%M2) / M3);
                g4 = (int)(tmpVec[i] % M3);

                if (!TAMask_B[g1][g2*W2+g3*W3+g4])
                    TAMask_B[g1][g2*W2+g3*W3+g4] = 1;

                x1_min_B = (g1 < x1_min_B) ? g1 : x1_min_B;
                x2_min_B = (g2 < x2_min_B) ? g2 : x2_min_B;
                x3_min_B = (g3 < x3_min_B) ? g3 : x3_min_B;
                x4_min_B = (g4 < x4_min_B) ? g4 : x4_min_B;
                x1_max_B = (g1 > x1_max_B) ? g1 : x1_max_B;
                x2_max_B = (g2 > x2_max_B) ? g2 : x2_max_B;
                x3_max_B = (g3 > x3_max_B) ? g3 : x3_max_B;
                x4_max_B = (g4 > x4_max_B) ? g4 : x4_max_B;
            }
            tmpVec.clear();

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-6) = %.4e sec\n\n", t_1_elapsed); 
            t_1_begin = omp_get_wtime();

            // Update ta_size

            ta_size_B = 0;

            #pragma omp parallel for reduction(+: ta_size_B) schedule(runtime)
            for (int i1 = x1_min_B; i1 <= x1_max_B; i1 ++)  {
                for (int i2 = x2_min_B; i2 <= x2_max_B; i2 ++)  {
                    for (int i3 = x3_min_B; i3 <= x3_max_B; i3 ++)  {
                        for (int i4 = x4_min_B; i4 <= x4_max_B; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])
                                ta_size_B += 1;    
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            if (!QUIET && TIMING) log->log("[Nonadiabatic4d] Elapsed time (initializing truncation B-7) = %.4e sec\n\n", t_1_elapsed); 
        }
        log->log("[Nonadiabatic4d] TA size (B) = %d, TB size (B) = %d\n", ta_size_B, tb_size_B);

        if (ta_size_B != 0)
            log->log("[Nonadiabatic4d] TA Range (B) [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min_B, x1_max_B, x2_min_B, x2_max_B, x3_min_B, x3_max_B, x4_min_B, x4_max_B);

        if (ta_size_A == 0 && ta_size_B == 0)  {
            log->log("[Nonadiabatic4d] STOP: Both TA_A and TA_B are empty\n");
            std::exit(EXIT_FAILURE);
        }

        // Box range union
        x1_max_C = (x1_max_A > x1_max_B) ? x1_max_A : x1_max_B;
        x1_min_C = (x1_min_A < x1_min_B) ? x1_min_A : x1_min_B;
        x2_max_C = (x2_max_A > x2_max_B) ? x2_max_A : x2_max_B;
        x2_min_C = (x2_min_A < x2_min_B) ? x2_min_A : x2_min_B;
        x3_max_C = (x3_max_A > x3_max_B) ? x3_max_A : x3_max_B;
        x3_min_C = (x3_min_A < x3_min_B) ? x3_min_A : x3_min_B;
        x4_max_C = (x4_max_A > x4_max_B) ? x4_max_A : x4_max_B;
        x4_min_C = (x4_min_A < x4_min_B) ? x4_min_A : x4_min_B;

        log->log("[Nonadiabatic4d] Union Range [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min_C, x1_max_C, x2_min_C, x2_max_C, x3_min_C, x3_max_C, x4_min_C, x4_max_C);

        t_2_end = omp_get_wtime();
        t_2_elapsed = t_2_end - t_2_begin;
        t_overhead += t_2_elapsed;

        if (!QUIET && TIMING)  {
            log->log("[Nonadiabatic4d] Elapsed time (initial truncation) = %.4e sec\n\n", t_1_elapsed);
            log->log("[Nonadiabatic4d] Initialization core computation time: %.4e sec\n", t_truncate);            
            log->log("[Nonadiabatic4d] Initialization overhead: %.4e sec\n", t_overhead); 
        }
    }
    else  // Full grid approach
    {
        log->log("[Nonadiabatic4d] Initialization core computation time: %.4e sec\n", t_full);   
    }
    // .........................................................................................

    // Time iteration 

    log->log("=======================================================\n\n"); 
    log->log("[Nonadiabatic4d] Time interation starts ...\n"); 
    log->log("[Nonadiabatic4d] Number of steps = %d\n\n", (int)(TIME / kk)); 
    log->log("=======================================================\n\n"); 

    for (int tt = 0; tt < (int)(TIME / kk); tt ++)
    {
        t_0_begin = omp_get_wtime(); 

        Excount_A = 0;
        Excount_B = 0;

        if ( tt % PRINT_PERIOD == 0 )
        {
            if ( isPrintEdge  && !isFullGrid )  {

                /* ---------------- Surface A ---------------- */

                pfile = fopen ("edge_A.dat","a");
                fprintf(pfile, "%d %lf %lu\n", tt, tt * kk, TB_A.size());

                for (int i = 0; i < TB_A.size(); i++)
                {
                    g1 = (int)(TB_A[i] / M1);
                    g2 = (int)((TB_A[i]%M1) / M2);
                    g3 = (int)((TB_A[i]%M2) / M3);
                    g4 = (int)(TB_A[i] % M3);

                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    xx3 = Box[4] + g3 * H[2];
                    xx4 = Box[6] + g4 * H[3];
                    fprintf(pfile, "%d %d %d %d %lf %lf %lf %lf\n", g1,g2,g3,g4,xx1,xx2,xx3,xx4);            
                }
                fclose(pfile);

                /* ---------------- Surface B ---------------- */

                pfile = fopen ("edge_B.dat","a");
                fprintf(pfile, "%d %lf %lu\n", tt, tt * kk, TB_B.size());

                for (int i = 0; i < TB_B.size(); i++)
                {
                    g1 = (int)(TB_B[i] / M1);
                    g2 = (int)((TB_B[i]%M1) / M2);
                    g3 = (int)((TB_B[i]%M2) / M3);
                    g4 = (int)(TB_B[i] % M3);

                    xx1 = Box[0] + g1 * H[0];
                    xx2 = Box[2] + g2 * H[1];
                    xx3 = Box[4] + g3 * H[2];
                    xx4 = Box[6] + g4 * H[3];
                    fprintf(pfile, "%d %d %d %d %lf %lf %lf %lf\n", g1,g2,g3,g4,xx1,xx2,xx3,xx4);       
                }
                fclose(pfile);
            }
            if ( isPrintDensity && !isFullGrid )  {

                /* ---------------- Surface A ---------------- */

                pfile = fopen ("density_A.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, ta_size_A);

                for (int i1 = x1_min_A; i1 <= x1_max_A; i1 ++)  {
                    for (int i2 = x2_min_A; i2 <= x2_max_A; i2 ++)  {
                        for (int i3 = x3_min_A; i3 <= x3_max_A; i3 ++)  {
                            for (int i4 = x4_min_A; i4 <= x4_max_A; i4 ++)  {
                                if (TAMask_A[i1][i2*W2+i3*W3+i4])  {
                                    xx1 = Box[0] + i1 * H[0];
                                    xx2 = Box[2] + i2 * H[1];
                                    xx3 = Box[4] + i3 * H[2];
                                    xx4 = Box[6] + i4 * H[3];
                                    fprintf(pfile, "%d %d %d %d %lf %lf %lf %lf %.16e\n", i1, i2, i3, i4, xx1, xx2, xx3, xx4, PF_A[i1][i2*W2+i3*W3+i4]);
                                }
                            }
                        }
                    }
                }
                fclose(pfile);

                /* ---------------- Surface B ---------------- */

                pfile = fopen ("density_B.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, ta_size_B);

                for (int i1 = x1_min_B; i1 <= x1_max_B; i1 ++)  {
                    for (int i2 = x2_min_B; i2 <= x2_max_B; i2 ++)  {
                        for (int i3 = x3_min_B; i3 <= x3_max_B; i3 ++)  {
                            for (int i4 = x4_min_B; i4 <= x4_max_B; i4 ++)  {
                                if (TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                    xx1 = Box[0] + i1 * H[0];
                                    xx2 = Box[2] + i2 * H[1];
                                    xx3 = Box[4] + i3 * H[2];
                                    xx4 = Box[6] + i4 * H[3];
                                    fprintf(pfile, "%d %d %d %d %lf %lf %lf %lf %.16e\n", i1, i2, i3, i4, xx1, xx2, xx3, xx4, PF_B[i1][i2*W2+i3*W3+i4]);
                                }
                            }
                        }
                    }
                }
                fclose(pfile);
            }
            if ( isPrintDensity && isFullGrid )  {

                /* ---------------- Surface A ---------------- */

                pfile = fopen ("density_A.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, O1);

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {
                            for (int i4 = 0; i4 < BoxShape[3]; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];
                                fprintf(pfile, "%d %d %d %d %lf %lf %lf %lf %.16e\n", i1, i2, i3, i4, xx1, xx2, xx3, xx4, PF_A[i1][i2*W2+i3*W3+i4]);
                            }
                        }
                    }
                }
                fclose(pfile);

                /* ---------------- Surface B ---------------- */

                pfile = fopen ("density_B.dat","a");
                fprintf(pfile, "%d %lf %d\n", tt, tt * kk, O1);

                for (int i1 = 0; i1 < BoxShape[0]; i1 ++)  {
                    for (int i2 = 0; i2 < BoxShape[1]; i2 ++)  {
                        for (int i3 = 0; i3 < BoxShape[2]; i3 ++)  {
                            for (int i4 = 0; i4 < BoxShape[3]; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];
                                fprintf(pfile, "%d %d %d %d %lf %lf %lf %lf %.16e\n", i1, i2, i3, i4, xx1, xx2, xx3, xx4, PF_B[i1][i2*W2+i3*W3+i4]);
                            }
                        }
                    }
                }
                fclose(pfile);  
            }
        }

        // Check if TB of f is higher than TolL
        
        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();
            t_truncate = 0.0;
            t_overhead = 0.0;

            // TBL = Index of Extrapolating Edge points.
            // TBL_P = Index history of TBL of this iteration. To Prevent from extrapolating the same points multiple times.

            /* ---------------- Surface A ---------------- */

            TBL_A.clear();

            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4,b1,b2,b3,b4,nx1,nx2,nx3,nx4,\
                                                                      f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A) schedule(runtime)
            for (int i = 0; i < TB_A.size(); i++)
            {
                g1 = (int)(TB_A[i] / M1);
                g2 = (int)((TB_A[i]%M1) / M2);
                g3 = (int)((TB_A[i]%M2) / M3);
                g4 = (int)(TB_A[i] % M3);

                nx1 = int(TAMask_A[g1+1][g2*W2+g3*W3+g4]) + int(TAMask_A[g1-1][g2*W2+g3*W3+g4]);
                nx2 = int(TAMask_A[g1][(g2+1)*W2+g3*W3+g4]) + int(TAMask_A[g1][(g2-1)*W2+g3*W3+g4]);
                nx3 = int(TAMask_A[g1][g2*W2+(g3+1)*W3+g4]) + int(TAMask_A[g1][g2*W2+(g3-1)*W3+g4]);
                nx4 = int(TAMask_A[g1][g2*W2+g3*W3+(g4+1)]) + int(TAMask_A[g1][g2*W2+g3*W3+(g4-1)]);

                f1p1_A = (TAMask_A[g1+1][g2*W2+g3*W3+g4]) ? F_A[g1+1][g2*W2+g3*W3+g4] : F_A[g1][g2*W2+g3*W3+g4];
                f1m1_A = (TAMask_A[g1-1][g2*W2+g3*W3+g4]) ? F_A[g1-1][g2*W2+g3*W3+g4] : F_A[g1][g2*W2+g3*W3+g4];
                f2p1_A = (TAMask_A[g1][(g2+1)*W2+g3*W3+g4]) ? F_A[g1][(g2+1)*W2+g3*W3+g4] : F_A[g1][g2*W2+g3*W3+g4];
                f2m1_A = (TAMask_A[g1][(g2-1)*W2+g3*W3+g4]) ? F_A[g1][(g2-1)*W2+g3*W3+g4] : F_A[g1][g2*W2+g3*W3+g4];
                f3p1_A = (TAMask_A[g1][g2*W2+(g3+1)*W3+g4]) ? F_A[g1][g2*W2+(g3+1)*W3+g4] : F_A[g1][g2*W2+g3*W3+g4];
                f3m1_A = (TAMask_A[g1][g2*W2+(g3-1)*W3+g4]) ? F_A[g1][g2*W2+(g3-1)*W3+g4] : F_A[g1][g2*W2+g3*W3+g4];
                f4p1_A = (TAMask_A[g1][g2*W2+g3*W3+(g4+1)]) ? F_A[g1][g2*W2+g3*W3+(g4+1)] : F_A[g1][g2*W2+g3*W3+g4];
                f4m1_A = (TAMask_A[g1][g2*W2+g3*W3+(g4-1)]) ? F_A[g1][g2*W2+g3*W3+(g4-1)] : F_A[g1][g2*W2+g3*W3+g4];

                b1 = PF_A[g1][g2*W2+g3*W3+g4] >= TolL_A;
                b2 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_A - f1m1_A)/(nx1*H[0]),2)) + \
                     ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_A - f2m1_A)/(nx2*H[1]),2)) + \
                     ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_A - f3m1_A)/(nx3*H[2]),2)) + \
                     ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_A - f4m1_A)/(nx4*H[3]),2)) >= TolLd_A_sq;
                b3 = g1 > EDGE && g2 > EDGE && g3 > EDGE && g4 > EDGE;
                b4 = (g1 < BoxShape[0]-EDGE-1) && (g2 < BoxShape[1]-EDGE-1) && (g3 < BoxShape[2]-EDGE-1) && (g4 < BoxShape[3]-EDGE-1);

                if ((b1||b2) && b3 && b4)
                    tmpVec.push_back(TB_A[i]);
            }
            tmpVec.swap(TBL_A);
            tmpVec.clear();
            TBL_P_A = TBL_A;

            /* ---------------- Surface B ---------------- */

            TBL_B.clear();

            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4,b1,b2,b3,b4,nx1,nx2,nx3,nx4,\
                                                                      f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B) schedule(runtime)
            for (int i = 0; i < TB_B.size(); i++)
            {
                g1 = (int)(TB_B[i] / M1);
                g2 = (int)((TB_B[i]%M1) / M2);
                g3 = (int)((TB_B[i]%M2) / M3);
                g4 = (int)(TB_B[i] % M3);

                nx1 = int(TAMask_B[g1+1][g2*W2+g3*W3+g4]) + int(TAMask_B[g1-1][g2*W2+g3*W3+g4]);
                nx2 = int(TAMask_B[g1][(g2+1)*W2+g3*W3+g4]) + int(TAMask_B[g1][(g2-1)*W2+g3*W3+g4]);
                nx3 = int(TAMask_B[g1][g2*W2+(g3+1)*W3+g4]) + int(TAMask_B[g1][g2*W2+(g3-1)*W3+g4]);
                nx4 = int(TAMask_B[g1][g2*W2+g3*W3+(g4+1)]) + int(TAMask_B[g1][g2*W2+g3*W3+(g4-1)]);

                f1p1_B = (TAMask_B[g1+1][g2*W2+g3*W3+g4]) ? F_B[g1+1][g2*W2+g3*W3+g4] : F_B[g1][g2*W2+g3*W3+g4];
                f1m1_B = (TAMask_B[g1-1][g2*W2+g3*W3+g4]) ? F_B[g1-1][g2*W2+g3*W3+g4] : F_B[g1][g2*W2+g3*W3+g4];
                f2p1_B = (TAMask_B[g1][(g2+1)*W2+g3*W3+g4]) ? F_B[g1][(g2+1)*W2+g3*W3+g4] : F_B[g1][g2*W2+g3*W3+g4];
                f2m1_B = (TAMask_B[g1][(g2-1)*W2+g3*W3+g4]) ? F_B[g1][(g2-1)*W2+g3*W3+g4] : F_B[g1][g2*W2+g3*W3+g4];
                f3p1_B = (TAMask_B[g1][g2*W2+(g3+1)*W3+g4]) ? F_B[g1][g2*W2+(g3+1)*W3+g4] : F_B[g1][g2*W2+g3*W3+g4];
                f3m1_B = (TAMask_B[g1][g2*W2+(g3-1)*W3+g4]) ? F_B[g1][g2*W2+(g3-1)*W3+g4] : F_B[g1][g2*W2+g3*W3+g4];
                f4p1_B = (TAMask_B[g1][g2*W2+g3*W3+(g4+1)]) ? F_B[g1][g2*W2+g3*W3+(g4+1)] : F_B[g1][g2*W2+g3*W3+g4];
                f4m1_B = (TAMask_B[g1][g2*W2+g3*W3+(g4-1)]) ? F_B[g1][g2*W2+g3*W3+(g4-1)] : F_B[g1][g2*W2+g3*W3+g4];

                b1 = PF_B[g1][g2*W2+g3*W3+g4] >= TolL_B;
                b2 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_B - f1m1_B)/(nx1*H[0]),2)) + \
                     ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_B - f2m1_B)/(nx2*H[1]),2)) + \
                     ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_B - f3m1_B)/(nx3*H[2]),2)) + \
                     ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_B - f4m1_B)/(nx4*H[3]),2)) >= TolLd_B_sq;
                b3 = g1 > EDGE && g2 > EDGE && g3 > EDGE && g4 > EDGE;
                b4 = (g1 < BoxShape[0]-EDGE-1) && (g2 < BoxShape[1]-EDGE-1) && (g3 < BoxShape[2]-EDGE-1) && (g4 < BoxShape[3]-EDGE-1);

                if ((b1||b2) && b3 && b4)
                    tmpVec.push_back(TB_B[i]);
            }
            tmpVec.swap(TBL_B);
            tmpVec.clear();
            TBL_P_B = TBL_B;              

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-a-1: TBL) = %.4e sec\n", t_1_elapsed);   
        }
        else  
        {
            t_full = 0.0;
        }
        isExtrapolate_A = false;
        isExtrapolate_B = false;
        isFirstExtrp_A = true;
        isFirstExtrp_B = true;

        // CASE 1: Truncating with extrapolation

        while ( !isFullGrid && (TBL_A.size() != 0 || TBL_B.size() != 0) && (Excount_A < ExLimit && Excount_B < ExLimit) )
        {
            // Extrapolation

            /* ---------------- Surface A ---------------- */

            if ( TBL_A.size() != 0 )  {

                t_1_begin = omp_get_wtime();

                isExtrapolate_A = true;

                // Avoid unexpected arrangement of TBL
                __gnu_parallel::sort(TBL_A.begin(),TBL_A.end());
                it = std::unique (TBL_A.begin(), TBL_A.end()); 
                TBL_A.resize(std::distance(TBL_A.begin(),it)); 

                // Find extrapolation target
                // ExFF: Index of Extrapolated points
                ExFF_A.clear();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
                for (int i = 0; i < TBL_A.size(); i++)  {

                    g1 = (int)(TBL_A[i] / M1);
                    g2 = (int)((TBL_A[i]%M1) / M2);
                    g3 = (int)((TBL_A[i]%M2) / M3);
                    g4 = (int)(TBL_A[i] % M3);

                    if ( g1-1 > EDGE && !TAMask_A[g1-1][g2*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                    if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_A[g1+1][g2*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                    if ( g2-1 > EDGE && !TAMask_A[g1][(g2-1)*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                    if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_A[g1][(g2+1)*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                    if ( g3-1 > EDGE && !TAMask_A[g1][g2*W2+(g3-1)*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                    if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_A[g1][g2*W2+(g3+1)*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                    if ( g4-1 > EDGE && !TAMask_A[g1][g2*W2+g3*W3+(g4-1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                    if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_A[g1][g2*W2+g3*W3+(g4+1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                }
                tmpVec.swap(ExFF_A);
                tmpVec.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-b-1: ExFF A) = %.4e sec\n", t_1_elapsed);  

                if ( ExFF_A.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    // ExFF & TBL set difference
                    tmpVec.resize(ExFF_A.size() + TBL_A.size());
                    __gnu_parallel::sort(TBL_A.begin(), TBL_A.end());
                    __gnu_parallel::sort(ExFF_A.begin(), ExFF_A.end());
                    it = std::set_difference(ExFF_A.begin(), ExFF_A.end(), TBL_A.begin(), TBL_A.end(), tmpVec.begin());
                    tmpVec.resize(it - tmpVec.begin()); 
                    tmpVec.swap(ExFF_A);
                    tmpVec.clear();

                    // Find unique elements
                    __gnu_parallel::sort(ExFF_A.begin(),ExFF_A.end());
                    it = std::unique (ExFF_A.begin(), ExFF_A.end()); 
                    ExFF_A.resize(std::distance(ExFF_A.begin(),it));

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-b-2: ExFF A) = %.4e sec\n", t_1_elapsed);   

                    // Find the direction of Outer to Edge points
                    t_1_begin = omp_get_wtime();

                    Check_A.clear();
                    ExTBL.clear();

                    for (int i = 0; i < ExFF_A.size(); i ++)  {
                        Check_A.push_back(true);
                        ExTBL.push_back(xZERO);
                    }
                    
                    #pragma omp parallel for private(g1,g2,g3,g4,sum,count,isEmpty,val_min_abs,val_min,min_dir) schedule(runtime)
                    for (int i = 0; i < ExFF_A.size(); i ++)
                    {
                        g1 = (int)(ExFF_A[i] / M1);
                        g2 = (int)((ExFF_A[i]%M1) / M2);
                        g3 = (int)((ExFF_A[i]%M2) / M3);
                        g4 = (int)(ExFF_A[i] % M3);
                        sum = xZERO;
                        count = 0;
                        isEmpty = true;
                        val_min_abs = BIG_NUMBER;
                        val_min = BIG_NUMBER;
                        min_dir = -1;

                        if ( F_A[g1-1][g2*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_A[g1-1][g2*W2+g3*W3+g4]) < val_min_abs &&  F_A[g1-2][g2*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1-1][g2*W2+g3*W3+g4]);
                                val_min = F_A[g1-1][g2*W2+g3*W3+g4];
                                min_dir = 0;
                            }
                            if ( F_A[g1-2][g2*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1-1][g2*W2+g3*W3+g4]) - std::log(F_A[g1-2][g2*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1+1][g2*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_A[g1+1][g2*W2+g3*W3+g4]) < val_min_abs && F_A[g1+2][g2*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1+1][g2*W2+g3*W3+g4]);
                                val_min = F_A[g1+1][g2*W2+g3*W3+g4];
                                min_dir = 0;
                            }
                            if ( F_A[g1+2][g2*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1+1][g2*W2+g3*W3+g4]) - std::log(F_A[g1+2][g2*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1][(g2-1)*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_A[g1][(g2-1)*W2+g3*W3+g4]) < val_min_abs && F_A[g1][(g2-2)*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1][(g2-1)*W2+g3*W3+g4]);
                                val_min = F_A[g1][(g2-1)*W2+g3*W3+g4];
                                min_dir = 1;
                            }
                            if ( F_A[g1][(g2-2)*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1][(g2-1)*W2+g3*W3+g4]) - std::log(F_A[g1][(g2-2)*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1][(g2+1)*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_A[g1][(g2+1)*W2+g3*W3+g4]) < val_min_abs && F_A[g1][(g2+2)*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1][(g2+1)*W2+g3*W3+g4]);
                                val_min = F_A[g1][(g2+1)*W2+g3*W3+g4];
                                min_dir = 1;
                            }
                            if ( F_A[g1][(g2+2)*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1][(g2+1)*W2+g3*W3+g4]) - std::log(F_A[g1][(g2+2)*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1][g2*W2+(g3-1)*W3+g4] != xZERO )  {

                            if ( std::abs(F_A[g1][g2*W2+(g3-1)*W3+g4]) < val_min_abs && F_A[g1][g2*W2+(g3-2)*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1][g2*W2+(g3-1)*W3+g4]);
                                val_min = F_A[g1][g2*W2+(g3-1)*W3+g4];
                                min_dir = 2;
                            }
                            if ( F_A[g1][g2*W2+(g3-2)*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1][g2*W2+(g3-1)*W3+g4]) - std::log(F_A[g1][g2*W2+(g3-2)*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1][g2*W2+(g3+1)*W3+g4] != xZERO )  {

                            if ( std::abs(F_A[g1][g2*W2+(g3+1)*W3+g4]) < val_min_abs && F_A[g1][g2*W2+(g3+2)*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1][g2*W2+(g3+1)*W3+g4]);
                                val_min = F_A[g1][g2*W2+(g3+1)*W3+g4];
                                min_dir = 2;
                            }
                            if ( F_A[g1][g2*W2+(g3+2)*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1][g2*W2+(g3+1)*W3+g4]) - std::log(F_A[g1][g2*W2+(g3+2)*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1][g2*W2+g3*W3+(g4-1)] != xZERO )  {

                            if ( std::abs(F_A[g1][g2*W2+g3*W3+(g4-1)]) < val_min_abs && F_A[g1][g2*W2+g3*W3+(g4-2)] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1][g2*W2+g3*W3+(g4-1)]);
                                val_min = F_A[g1][g2*W2+g3*W3+(g4-1)];
                                min_dir = 3;
                            }
                            if ( F_A[g1][g2*W2+g3*W3+(g4-2)] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1][g2*W2+g3*W3+(g4-1)]) - std::log(F_A[g1][g2*W2+g3*W3+(g4-2)]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_A[g1][g2*W2+g3*W3+(g4+1)] != xZERO )  {

                            if ( std::abs(F_A[g1][g2*W2+g3*W3+(g4+1)]) < val_min_abs && F_A[g1][g2*W2+g3*W3+(g4+2)] != xZERO )  {
                                val_min_abs = std::abs(F_A[g1][g2*W2+g3*W3+(g4+1)]);
                                val_min = F_A[g1][g2*W2+g3*W3+(g4+1)];
                                min_dir = 3;
                            }
                            if ( F_A[g1][g2*W2+g3*W3+(g4+2)] != xZERO )  {

                                val = exp( 2.0 * std::log(F_A[g1][g2*W2+g3*W3+(g4+1)]) - std::log(F_A[g1][g2*W2+g3*W3+(g4+2)]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( isEmpty )
                        {
                            Check_A[i] = false; 

                        }  else  {

                            // Assume the probability of outer points always smaller than edge point.
                            // if larger, then set the edge point with smallest P to outer
                            // point instead of using the extrapolation result.
                            if ( std::abs(sum/count) > val_min_abs )  {
                                ExTBL[i] = val_min * exp(-ExReduce * H[min_dir]);
                            }
                            else  {
                                ExTBL[i] = sum / count;
                            }
                        }
                    }

                    count = 0;

                    #pragma omp parallel for reduction (+:count) private(g1,g2) schedule(runtime)
                    for ( int i = 0; i < ExFF_A.size(); i++ )  {
                        if (Check_A[i])  {
                            g1 = (int)(ExFF_A[i] / M1);
                            g2 = (int)(ExFF_A[i] % M1);
                            F_A[g1][g2] = ExTBL[i];
                            count += 1;
                        }
                    }

                    if (count == 0)  {
                        ExFF_A.clear();
                        ExTBL.clear();
                        Check_A.clear();
                    }
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-b-3: ExFF A) = %.4e sec\n", t_1_elapsed);  
                }  // ExFF_A.size() > 0
            }  // if TBL_A.size() != 0

            /* ---------------- Surface B ---------------- */

            if ( TBL_B.size() != 0 )  {

                t_1_begin = omp_get_wtime();

                isExtrapolate_B = true;

                // Avoid unexpected arrangement of TBL
                __gnu_parallel::sort(TBL_B.begin(),TBL_B.end());
                it = std::unique (TBL_B.begin(), TBL_B.end()); 
                TBL_B.resize(std::distance(TBL_B.begin(),it)); 

                // Find extrapolation target
                // ExFF: Index of Extrapolated points
                ExFF_B.clear();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
                for (int i = 0; i < TBL_B.size(); i++)  {

                    g1 = (int)(TBL_B[i] / M1);
                    g2 = (int)((TBL_B[i]%M1) / M2);
                    g3 = (int)((TBL_B[i]%M2) / M3);
                    g4 = (int)(TBL_B[i] % M3);

                    if ( g1-1 > EDGE && !TAMask_B[g1-1][g2*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                    if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_B[g1+1][g2*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                    if ( g2-1 > EDGE && !TAMask_B[g1][(g2-1)*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                    if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_B[g1][(g2+1)*W2+g3*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                    if ( g3-1 > EDGE && !TAMask_B[g1][g2*W2+(g3-1)*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                    if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_B[g1][g2*W2+(g3+1)*W3+g4] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                    if ( g4-1 > EDGE && !TAMask_B[g1][g2*W2+g3*W3+(g4-1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                    if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_B[g1][g2*W2+g3*W3+(g4+1)] )
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                }
                tmpVec.swap(ExFF_B);
                tmpVec.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-b-1: ExFF B) = %.4e sec\n", t_1_elapsed);  

                if ( ExFF_B.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    // ExFF & TBL set difference
                    tmpVec.resize(ExFF_B.size() + TBL_B.size());
                    __gnu_parallel::sort(TBL_B.begin(), TBL_B.end());
                    __gnu_parallel::sort(ExFF_B.begin(), ExFF_B.end());
                    it = std::set_difference( ExFF_B.begin(), ExFF_B.end(), TBL_B.begin(), TBL_B.end(), tmpVec.begin() );
                    tmpVec.resize(it - tmpVec.begin()); 
                    tmpVec.swap(ExFF_B);
                    tmpVec.clear();

                    // Find unique elements
                    __gnu_parallel::sort(ExFF_B.begin(),ExFF_B.end());
                    it = std::unique (ExFF_B.begin(), ExFF_B.end()); 
                    ExFF_B.resize(std::distance(ExFF_B.begin(),it));

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-b-2: ExFF B) = %.4e sec\n", t_1_elapsed);   

                    // Find the direction of "outer point to edge point"

                    t_1_begin = omp_get_wtime();

                    Check_B.clear();
                    ExTBL.clear();

                    for (int i = 0; i < ExFF_B.size(); i ++)  {
                        Check_B.push_back(true);
                        ExTBL.push_back(xZERO);
                    }
                    
                    #pragma omp parallel for private(g1,g2,g3,g4,sum,count,isEmpty,val_min_abs,val_min,min_dir) schedule(runtime)
                    for (int i = 0; i < ExFF_B.size(); i ++)
                    {
                        g1 = (int)(ExFF_B[i] / M1);
                        g2 = (int)((ExFF_B[i]%M1) / M2);
                        g3 = (int)((ExFF_B[i]%M2) / M3);
                        g4 = (int)(ExFF_B[i] % M3);
                        sum = xZERO;
                        count = 0;
                        isEmpty = true;
                        val_min_abs = BIG_NUMBER;
                        val_min = BIG_NUMBER;
                        min_dir = -1;

                        if ( F_B[g1-1][g2*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_B[g1-1][g2*W2+g3*W3+g4]) < val_min_abs &&  F_B[g1-2][g2*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1-1][g2*W2+g3*W3+g4]);
                                val_min = F_B[g1-1][g2*W2+g3*W3+g4];
                                min_dir = 0;
                            }
                            if ( F_B[g1-2][g2*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1-1][g2*W2+g3*W3+g4]) - std::log(F_B[g1-2][g2*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1+1][g2*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_B[g1+1][g2*W2+g3*W3+g4]) < val_min_abs && F_B[g1+2][g2*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1+1][g2*W2+g3*W3+g4]);
                                val_min = F_B[g1+1][g2*W2+g3*W3+g4];
                                min_dir = 0;
                            }
                            if ( F_B[g1+2][g2*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1+1][g2*W2+g3*W3+g4]) - std::log(F_B[g1+2][g2*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1][(g2-1)*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_B[g1][(g2-1)*W2+g3*W3+g4]) < val_min_abs && F_B[g1][(g2-2)*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1][(g2-1)*W2+g3*W3+g4]);
                                val_min = F_B[g1][(g2-1)*W2+g3*W3+g4];
                                min_dir = 1;
                            }
                            if ( F_B[g1][(g2-2)*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1][(g2-1)*W2+g3*W3+g4]) - std::log(F_B[g1][(g2-2)*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1][(g2+1)*W2+g3*W3+g4] != xZERO )  {

                            if ( std::abs(F_B[g1][(g2+1)*W2+g3*W3+g4]) < val_min_abs && F_B[g1][(g2+2)*W2+g3*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1][(g2+1)*W2+g3*W3+g4]);
                                val_min = F_B[g1][(g2+1)*W2+g3*W3+g4];
                                min_dir = 1;
                            }
                            if ( F_B[g1][(g2+2)*W2+g3*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1][(g2+1)*W2+g3*W3+g4]) - std::log(F_B[g1][(g2+2)*W2+g3*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1][g2*W2+(g3-1)*W3+g4] != xZERO )  {

                            if ( std::abs(F_B[g1][g2*W2+(g3-1)*W3+g4]) < val_min_abs && F_B[g1][g2*W2+(g3-2)*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1][g2*W2+(g3-1)*W3+g4]);
                                val_min = F_B[g1][g2*W2+(g3-1)*W3+g4];
                                min_dir = 2;
                            }
                            if ( F_B[g1][g2*W2+(g3-2)*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1][g2*W2+(g3-1)*W3+g4]) - std::log(F_B[g1][g2*W2+(g3-2)*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1][g2*W2+(g3+1)*W3+g4] != xZERO )  {

                            if ( std::abs(F_B[g1][g2*W2+(g3+1)*W3+g4]) < val_min_abs && F_B[g1][g2*W2+(g3+2)*W3+g4] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1][g2*W2+(g3+1)*W3+g4]);
                                val_min = F_B[g1][g2*W2+(g3+1)*W3+g4];
                                min_dir = 2;
                            }
                            if ( F_B[g1][g2*W2+(g3+2)*W3+g4] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1][g2*W2+(g3+1)*W3+g4]) - std::log(F_B[g1][g2*W2+(g3+2)*W3+g4]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1][g2*W2+g3*W3+(g4-1)] != xZERO )  {

                            if ( std::abs(F_B[g1][g2*W2+g3*W3+(g4-1)]) < val_min_abs && F_B[g1][g2*W2+g3*W3+(g4-2)] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1][g2*W2+g3*W3+(g4-1)]);
                                val_min = F_B[g1][g2*W2+g3*W3+(g4-1)];
                                min_dir = 3;
                            }
                            if ( F_B[g1][g2*W2+g3*W3+(g4-2)] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1][g2*W2+g3*W3+(g4-1)]) - std::log(F_B[g1][g2*W2+g3*W3+(g4-2)]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( F_B[g1][g2*W2+g3*W3+(g4+1)] != xZERO )  {

                            if ( std::abs(F_B[g1][g2*W2+g3*W3+(g4+1)]) < val_min_abs && F_B[g1][g2*W2+g3*W3+(g4+2)] != xZERO )  {
                                val_min_abs = std::abs(F_B[g1][g2*W2+g3*W3+(g4+1)]);
                                val_min = F_B[g1][g2*W2+g3*W3+(g4+1)];
                                min_dir = 3;
                            }
                            if ( F_B[g1][g2*W2+g3*W3+(g4+2)] != xZERO )  {

                                val = exp( 2.0 * std::log(F_B[g1][g2*W2+g3*W3+(g4+1)]) - std::log(F_B[g1][g2*W2+g3*W3+(g4+2)]) );

                                if ( !(isnan(real(val)) || isnan(-real(val))) && !(isinf(real(val)) || isinf(-real(val))) && !(isnan(imag(val)) || isnan(-imag(val))) && !(isinf(imag(val)) || isinf(-imag(val))) )  
                                {
                                    sum += val;
                                    count += 1;
                                    isEmpty = false;
                                }
                            }
                        }

                        if ( isEmpty )
                        {
                            Check_B[i] = false;

                        }  else  {

                            if ( std::abs(sum/count) > val_min_abs )  {
                                ExTBL[i] = val_min * exp(-ExReduce * H[min_dir]);
                            }
                            else  {
                                ExTBL[i] = sum / count;
                            }
                        }
                    }

                    count = 0;

                    #pragma omp parallel for reduction (+:count) private(g1,g2) schedule(runtime)
                    for ( int i = 0; i < ExFF_B.size(); i++ )  {
                        if (Check_B[i])  {
                            g1 = (int)(ExFF_B[i] / M1);
                            g2 = (int)(ExFF_B[i] % M1);
                            F_B[g1][g2] = ExTBL[i];
                            count += 1;
                        }
                    }

                    if (count == 0)  {
                        ExFF_B.clear();
                        ExTBL.clear();
                        Check_B.clear();
                    }
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-b-3: ExFF B) = %.4e sec\n", t_1_elapsed); 
                }  // ExFF_B.size() > 0
            }  // if TBL_B.size() != 0

            // ............................................................................................. Extrapolation

            if ( isFirstExtrp_A && isFirstExtrp_B )  {

                // Check Extending nonzero Area

                if ( ExFF_A.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    tmpVec.clear();

                    #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
                    for (int i = 0; i < ExFF_A.size(); i++)  {
                        if (Check_A[i])  
                        {
                            g1 = (int)(ExFF_A[i] / M1);
                            g2 = (int)((ExFF_A[i]%M1) / M2);
                            g3 = (int)((ExFF_A[i]%M2) / M3);
                            g4 = (int)(ExFF_A[i] % M3);

                            if ( !TAMask_A[g1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4));
                            if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_A[g1+1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                            if ( g1-1 > EDGE && !TAMask_A[g1-1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                            if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_A[g1][(g2+1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                            if ( g2-1 > EDGE && !TAMask_A[g1][(g2-1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                            if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_A[g1][g2*W2+(g3+1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                            if ( g3-1 > EDGE && !TAMask_A[g1][g2*W2+(g3-1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                            if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_A[g1][g2*W2+g3*W3+(g4+1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                            if ( g4-1 > EDGE && !TAMask_A[g1][g2*W2+g3*W3+(g4-1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                        }
                    }

                    for (int i = 0; i < tmpVec.size(); i ++)  {
                        g1 = (int)(tmpVec[i] / M1);
                        g2 = (int)((tmpVec[i]%M1) / M2);
                        g3 = (int)((tmpVec[i]%M2) / M3);
                        g4 = (int)(tmpVec[i] % M3);
                        TAMask_A[g1][g2*W2+g3*W3+g4] = 1;
                    
                        x1_min_A = (g1 < x1_min_A) ? g1 : x1_min_A;
                        x1_max_A = (g1 > x1_max_A) ? g1 : x1_max_A;
                        x2_min_A = (g2 < x2_min_A) ? g2 : x2_min_A;
                        x2_max_A = (g2 > x2_max_A) ? g2 : x2_max_A;
                        x3_min_A = (g3 < x3_min_A) ? g3 : x3_min_A;
                        x3_max_A = (g3 > x3_max_A) ? g3 : x3_max_A;
                        x4_min_A = (g4 < x4_min_A) ? g4 : x4_min_A;
                        x4_max_A = (g4 > x4_max_A) ? g4 : x4_max_A;
                    }
                    tmpVec.clear();

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-c-1: CASE 1 TA A) = %.4e sec\n", t_1_elapsed);
                }

                if ( ExFF_B.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    tmpVec.clear();

                    #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
                    for (int i = 0; i < ExFF_B.size(); i++)  {
                        if (Check_B[i])  
                        {
                            g1 = (int)(ExFF_B[i] / M1);
                            g2 = (int)((ExFF_B[i]%M1) / M2);
                            g3 = (int)((ExFF_B[i]%M2) / M3);
                            g4 = (int)(ExFF_B[i] % M3);

                            if ( !TAMask_B[g1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4));
                            if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_B[g1+1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                            if ( g1-1 > EDGE && !TAMask_B[g1-1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                            if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_B[g1][(g2+1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                            if ( g2-1 > EDGE && !TAMask_B[g1][(g2-1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                            if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_B[g1][g2*W2+(g3+1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                            if ( g3-1 > EDGE && !TAMask_B[g1][g2*W2+(g3-1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                            if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_B[g1][g2*W2+g3*W3+(g4+1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                            if ( g4-1 > EDGE && !TAMask_B[g1][g2*W2+g3*W3+(g4-1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                        }
                    }

                    for (int i = 0; i < tmpVec.size(); i ++)  {
                        g1 = (int)(tmpVec[i] / M1);
                        g2 = (int)((tmpVec[i]%M1) / M2);
                        g3 = (int)((tmpVec[i]%M2) / M3);
                        g4 = (int)(tmpVec[i] % M3);
                        TAMask_B[g1][g2*W2+g3*W3+g4] = 1;
                    
                        x1_min_B = (g1 < x1_min_B) ? g1 : x1_min_B;
                        x1_max_B = (g1 > x1_max_B) ? g1 : x1_max_B;
                        x2_min_B = (g2 < x2_min_B) ? g2 : x2_min_B;
                        x2_max_B = (g2 > x2_max_B) ? g2 : x2_max_B;
                        x3_min_B = (g3 < x3_min_B) ? g3 : x3_min_B;
                        x3_max_B = (g3 > x3_max_B) ? g3 : x3_max_B;
                        x4_min_B = (g4 < x4_min_B) ? g4 : x4_min_B;
                        x4_max_B = (g4 > x4_max_B) ? g4 : x4_max_B;
                    }
                    tmpVec.clear();

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-c-1: CASE 1 TA B) = %.4e sec\n", t_1_elapsed); 
                }

                // Box range union
                x1_max_C = (x1_max_A > x1_max_B) ? x1_max_A : x1_max_B;
                x1_min_C = (x1_min_A < x1_min_B) ? x1_min_A : x1_min_B;
                x2_max_C = (x2_max_A > x2_max_B) ? x2_max_A : x2_max_B;
                x2_min_C = (x2_min_A < x2_min_B) ? x2_min_A : x2_min_B;
                x3_max_C = (x3_max_A > x3_max_B) ? x3_max_A : x3_max_B;
                x3_min_C = (x3_min_A < x3_min_B) ? x3_min_A : x3_min_B;
                x4_max_C = (x4_max_A > x4_max_B) ? x4_max_A : x4_max_B;
                x4_min_C = (x4_min_A < x4_min_B) ? x4_min_A : x4_min_B;

                // Runge???Kutta 4
                #pragma omp parallel 
                {
                    #pragma omp single nowait
                    {
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-1
                    #pragma omp for private(xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                        for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                            for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                                for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                    if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                        xx1 = Box[0] + i1 * H[0];
                                        xx2 = Box[2] + i2 * H[1];
                                        xx3 = Box[4] + i3 * H[2];
                                        xx4 = Box[6] + i4 * H[3];

                                        f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                        f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];

                                        f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                        f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];

                                        c_f0_a = 2.0 * f0_A;
                                        c_f0_b = 2.0 * f0_B;
                                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                        KK1_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A + f1p1_A - c_f0_a) + 
                                                                    Ikh2mh1sq * (f2m1_A + f2p1_A - c_f0_a) +
                                                                    Ikh2mh2sq * (f3m1_A + f3p1_A - c_f0_a) + 
                                                                    Ikh2mh3sq * (f4m1_A + f4p1_A - c_f0_a) - 
                                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                        FF_A[i1][i2*W2+i3*W3+i4] = F_A[i1][i2*W2+i3*W3+i4] + KK1_A[i1][i2*W2+i3*W3+i4] / 6.0;

                                        KK1_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B + f1p1_B - c_f0_b) + 
                                                                    Ikh2mh1sq * (f2m1_B + f2p1_B - c_f0_b) +
                                                                    Ikh2mh2sq * (f3m1_B + f3p1_B - c_f0_b) + 
                                                                    Ikh2mh3sq * (f4m1_B + f4p1_B - c_f0_b) - 
                                                                    Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                        FF_B[i1][i2*W2+i3*W3+i4] = F_B[i1][i2*W2+i3*W3+i4] + KK1_B[i1][i2*W2+i3*W3+i4] / 6.0;
                                    }
                                }
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-11: CASE 1 KK1) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-2
                    #pragma omp for private(xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                        for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                            for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                                for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                    if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                        xx1 = Box[0] + i1 * H[0];
                                        xx2 = Box[2] + i2 * H[1];
                                        xx3 = Box[4] + i3 * H[2];
                                        xx4 = Box[6] + i4 * H[3];

                                        f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                        f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                        kk0_A = KK1_A[i1][i2*W2+i3*W3+i4];
                                        kk1p1_A = KK1_A[i1+1][i2*W2+i3*W3+i4];
                                        kk1m1_A = KK1_A[i1-1][i2*W2+i3*W3+i4];
                                        kk2p1_A = KK1_A[i1][(i2+1)*W2+i3*W3+i4];
                                        kk2m1_A = KK1_A[i1][(i2-1)*W2+i3*W3+i4];
                                        kk3p1_A = KK1_A[i1][i2*W2+(i3+1)*W3+i4];
                                        kk3m1_A = KK1_A[i1][i2*W2+(i3-1)*W3+i4];
                                        kk4p1_A = KK1_A[i1][i2*W2+i3*W3+(i4+1)];
                                        kk4m1_A = KK1_A[i1][i2*W2+i3*W3+(i4-1)];

                                        f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                        f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                        kk0_B = KK1_B[i1][i2*W2+i3*W3+i4];
                                        kk1p1_B = KK1_B[i1+1][i2*W2+i3*W3+i4];
                                        kk1m1_B = KK1_B[i1-1][i2*W2+i3*W3+i4];
                                        kk2p1_B = KK1_B[i1][(i2+1)*W2+i3*W3+i4];
                                        kk2m1_B = KK1_B[i1][(i2-1)*W2+i3*W3+i4];
                                        kk3p1_B = KK1_B[i1][i2*W2+(i3+1)*W3+i4];
                                        kk3m1_B = KK1_B[i1][i2*W2+(i3-1)*W3+i4];
                                        kk4p1_B = KK1_B[i1][i2*W2+i3*W3+(i4+1)];
                                        kk4m1_B = KK1_B[i1][i2*W2+i3*W3+(i4-1)];

                                        c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                                        c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                        KK2_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                                    Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) + 
                                                                    Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                                    Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                        FF_A[i1][i2*W2+i3*W3+i4] += KK2_A[i1][i2*W2+i3*W3+i4] / 3.0;

                                        KK2_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                                    Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                                    Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                                    Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                                    Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                        FF_B[i1][i2*W2+i3*W3+i4] += KK2_B[i1][i2*W2+i3*W3+i4] / 3.0;
                                    }
                                }
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-12: CASE 1 KK2) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-3
                    #pragma omp for private(xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                        for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                            for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                                for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                    if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                        xx1 = Box[0] + i1 * H[0];
                                        xx2 = Box[2] + i2 * H[1];
                                        xx3 = Box[4] + i3 * H[2];
                                        xx4 = Box[6] + i4 * H[3];

                                        f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                        f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                        kk0_A = KK2_A[i1][i2*W2+i3*W3+i4];
                                        kk1p1_A = KK2_A[i1+1][i2*W2+i3*W3+i4];
                                        kk1m1_A = KK2_A[i1-1][i2*W2+i3*W3+i4];
                                        kk2p1_A = KK2_A[i1][(i2+1)*W2+i3*W3+i4];
                                        kk2m1_A = KK2_A[i1][(i2-1)*W2+i3*W3+i4];
                                        kk3p1_A = KK2_A[i1][i2*W2+(i3+1)*W3+i4];
                                        kk3m1_A = KK2_A[i1][i2*W2+(i3-1)*W3+i4];
                                        kk4p1_A = KK2_A[i1][i2*W2+i3*W3+(i4+1)];
                                        kk4m1_A = KK2_A[i1][i2*W2+i3*W3+(i4-1)];

                                        f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                        f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                        kk0_B = KK2_B[i1][i2*W2+i3*W3+i4];
                                        kk1p1_B = KK2_B[i1+1][i2*W2+i3*W3+i4];
                                        kk1m1_B = KK2_B[i1-1][i2*W2+i3*W3+i4];
                                        kk2p1_B = KK2_B[i1][(i2+1)*W2+i3*W3+i4];
                                        kk2m1_B = KK2_B[i1][(i2-1)*W2+i3*W3+i4];
                                        kk3p1_B = KK2_B[i1][i2*W2+(i3+1)*W3+i4];
                                        kk3m1_B = KK2_B[i1][i2*W2+(i3-1)*W3+i4];
                                        kk4p1_B = KK2_B[i1][i2*W2+i3*W3+(i4+1)];
                                        kk4m1_B = KK2_B[i1][i2*W2+i3*W3+(i4-1)];

                                        c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                                        c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                        KK3_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                                Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) +
                                                                Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                                Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                                Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                        FF_A[i1][i2*W2+i3*W3+i4] += KK3_A[i1][i2*W2+i3*W3+i4] / 3.0;

                                        KK3_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                                Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                                Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                                Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                                Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                        FF_B[i1][i2*W2+i3*W3+i4] += KK3_B[i1][i2*W2+i3*W3+i4] / 3.0;
                                    }
                                }
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-13: CASE 1 KK3) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }

                    // RK4-4
                    #pragma omp for private(xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                        for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                            for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                                for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                    if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                        xx1 = Box[0] + i1 * H[0];
                                        xx2 = Box[2] + i2 * H[1];
                                        xx3 = Box[4] + i3 * H[2];
                                        xx4 = Box[6] + i4 * H[3];

                                        f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                        f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                        kk0_A = KK3_A[i1][i2*W2+i3*W3+i4];
                                        kk1p1_A = KK3_A[i1+1][i2*W2+i3*W3+i4];
                                        kk1m1_A = KK3_A[i1-1][i2*W2+i3*W3+i4];
                                        kk2p1_A = KK3_A[i1][(i2+1)*W2+i3*W3+i4];
                                        kk2m1_A = KK3_A[i1][(i2-1)*W2+i3*W3+i4];
                                        kk3p1_A = KK3_A[i1][i2*W2+(i3+1)*W3+i4];
                                        kk3m1_A = KK3_A[i1][i2*W2+(i3-1)*W3+i4];
                                        kk4p1_A = KK3_A[i1][i2*W2+i3*W3+(i4+1)];
                                        kk4m1_A = KK3_A[i1][i2*W2+i3*W3+(i4-1)];

                                        f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                        f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                        f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                        f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                        f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                        f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                        f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                        f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                        f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                        kk0_B = KK3_B[i1][i2*W2+i3*W3+i4];
                                        kk1p1_B = KK3_B[i1+1][i2*W2+i3*W3+i4];
                                        kk1m1_B = KK3_B[i1-1][i2*W2+i3*W3+i4];
                                        kk2p1_B = KK3_B[i1][(i2+1)*W2+i3*W3+i4];
                                        kk2m1_B = KK3_B[i1][(i2-1)*W2+i3*W3+i4];
                                        kk3p1_B = KK3_B[i1][i2*W2+(i3+1)*W3+i4];
                                        kk3m1_B = KK3_B[i1][i2*W2+(i3-1)*W3+i4];
                                        kk4p1_B = KK3_B[i1][i2*W2+i3*W3+(i4+1)];
                                        kk4m1_B = KK3_B[i1][i2*W2+i3*W3+(i4-1)];

                                        c_f0_a = 2.0 * (f0_A+kk0_A);
                                        c_f0_b = 2.0 * (f0_B+kk0_B);
                                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                        KK4_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+kk1m1_A + f1p1_A+kk1p1_A - c_f0_a) + 
                                                                    Ikh2mh1sq * (f2m1_A+kk2m1_A + f2p1_A+kk2p1_A - c_f0_a) + 
                                                                    Ikh2mh2sq * (f3m1_A+kk3m1_A + f3p1_A+kk3p1_A - c_f0_a) +
                                                                    Ikh2mh3sq * (f4m1_A+kk4m1_A + f4p1_A+kk4p1_A - c_f0_a) - 
                                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                        FF_A[i1][i2*W2+i3*W3+i4] += KK4_A[i1][i2*W2+i3*W3+i4] / 6.0;

                                        KK4_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+kk1m1_B + f1p1_B+kk1p1_B - c_f0_b) + 
                                                                    Ikh2mh1sq * (f2m1_B+kk2m1_B + f2p1_B+kk2p1_B - c_f0_b) + 
                                                                    Ikh2mh2sq * (f3m1_B+kk3m1_B + f3p1_B+kk3p1_B - c_f0_b) + 
                                                                    Ikh2mh3sq * (f4m1_B+kk4m1_B + f4p1_B+kk4p1_B - c_f0_b) - 
                                                                    Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                        FF_B[i1][i2*W2+i3*W3+i4] += KK4_B[i1][i2*W2+i3*W3+i4] / 6.0;
                                    }
                                }
                            }
                        }
                    }

                    #pragma omp single nowait
                    {
                        t_1_end = omp_get_wtime();
                        t_1_elapsed = t_1_end - t_1_begin;
                        t_truncate += t_1_elapsed;
                        if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-14: CASE 1 KK4) = %lf sec\n", t_1_elapsed);
                        t_1_begin = omp_get_wtime();
                    }
                } // omp parallel

                isFirstExtrp_A = false;
                isFirstExtrp_B = false;

            } // if ( isFirstExtrp_A && isFirstExtrp_B )
            else if (ExFF_A.size() == 0 && ExFF_B.size() == 0)  {

                // In the case no valid ExFF found, reset TBL_A and TBL_B to break the while loop 
                TBL_A.clear();
                TBL_B.clear();
            }
            else
            {
                // Extrapolation loop when multiple expanding occured

                /* ---------------- Surface A ---------------- */

                if ( ExFF_A.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    tmpVec.clear();

                    #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
                    for (int i = 0; i < ExFF_A.size(); i++)  {
                        if (Check_A[i])  
                        {
                            g1 = (int)(ExFF_A[i] / M1);
                            g2 = (int)((ExFF_A[i]%M1) / M2);
                            g3 = (int)((ExFF_A[i]%M2) / M3);
                            g4 = (int)(ExFF_A[i] % M3);

                            if ( !TAMask_A[g1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4));
                            if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_A[g1+1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                            if ( g1-1 > EDGE && !TAMask_A[g1-1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                            if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_A[g1][(g2+1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                            if ( g2-1 > EDGE && !TAMask_A[g1][(g2-1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                            if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_A[g1][g2*W2+(g3+1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                            if ( g3-1 > EDGE && !TAMask_A[g1][g2*W2+(g3-1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                            if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_A[g1][g2*W2+g3*W3+(g4+1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                            if ( g4-1 > EDGE && !TAMask_A[g1][g2*W2+g3*W3+(g4-1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                        }
                    }

                    for (int i = 0; i < tmpVec.size(); i ++)  {

                        g1 = (int)(tmpVec[i] / M1);
                        g2 = (int)((tmpVec[i]%M1) / M2);
                        g3 = (int)((tmpVec[i]%M2) / M3);
                        g4 = (int)(tmpVec[i] % M3);
                        TAMask_A[g1][g2*W2+g3*W3+g4] = 1;

                        x1_min_A = (g1 < x1_min_A) ? g1 : x1_min_A;
                        x1_max_A = (g1 > x1_max_A) ? g1 : x1_max_A;
                        x2_min_A = (g2 < x2_min_A) ? g2 : x2_min_A;
                        x2_max_A = (g2 > x2_max_A) ? g2 : x2_max_A;
                        x3_min_A = (g3 < x3_min_A) ? g3 : x3_min_A;
                        x3_max_A = (g3 > x3_max_A) ? g3 : x3_max_A;
                        x4_min_A = (g4 < x4_min_A) ? g4 : x4_min_A;
                        x4_max_A = (g4 > x4_max_A) ? g4 : x4_max_A;
                    }
                    tmpVec.clear();
                    ExBD_A.clear();

                    #pragma omp parallel for reduction(merge: ExBD_A) private(g1,g2,g3,g4,n1,n2,n3,n4) schedule(runtime)
                    for (int i = 0; i < ExFF_A.size(); i++)
                    {
                        if (Check_A[i])  
                        {
                            g1 = (int)(ExFF_A[i] / M1);
                            g2 = (int)((ExFF_A[i]%M1) / M2);
                            g3 = (int)((ExFF_A[i]%M2) / M3);
                            g4 = (int)(ExFF_A[i] % M3);

                            ExBD_A.push_back(ExFF_A[i]);

                            for (int j = 0; j < nneigh; j ++)  {

                                n1 = neighlist[j][0];
                                n2 = neighlist[j][1];
                                n3 = neighlist[j][2];
                                n4 = neighlist[j][3];

                                if (TAMask_A[g1+n1][(g2+n2)*W2+(g3+n3)*W3+(g4+n4)])
                                    ExBD_A.push_back(GridToIdx(g1+n1,g2+n2,g3+n3,g4+n4));
                            }
                        }
                    }

                    // Find unique elements (ExBD)
                    __gnu_parallel::sort(ExBD_A.begin(),ExBD_A.end());
                    it = std::unique (ExBD_A.begin(), ExBD_A.end()); 
                    ExBD_A.resize(std::distance(ExBD_A.begin(),it));

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-cx-1: CASE 1 ExBD A) = %.4e sec\n", t_1_elapsed); 

                }  // if ExFF_A.size() > 0

                /* ---------------- Surface B ---------------- */

                if ( ExFF_B.size() > 0 )  {

                    t_1_begin = omp_get_wtime();

                    tmpVec.clear();

                    #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
                    for (int i = 0; i < ExFF_B.size(); i++)  {
                        if (Check_B[i])  
                        {
                            g1 = (int)(ExFF_B[i] / M1);
                            g2 = (int)((ExFF_B[i]%M1) / M2);
                            g3 = (int)((ExFF_B[i]%M2) / M3);
                            g4 = (int)(ExFF_B[i] % M3);

                            if ( !TAMask_B[g1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4));
                            if ( g1+1 < BoxShape[0]-EDGE-1 && !TAMask_B[g1+1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                            if ( g1-1 > EDGE && !TAMask_B[g1-1][g2*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                            if ( g2+1 < BoxShape[1]-EDGE-1 && !TAMask_B[g1][(g2+1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                            if ( g2-1 > EDGE && !TAMask_B[g1][(g2-1)*W2+g3*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                            if ( g3+1 < BoxShape[2]-EDGE-1 && !TAMask_B[g1][g2*W2+(g3+1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                            if ( g3-1 > EDGE && !TAMask_B[g1][g2*W2+(g3-1)*W3+g4] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                            if ( g4+1 < BoxShape[3]-EDGE-1 && !TAMask_B[g1][g2*W2+g3*W3+(g4+1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                            if ( g4-1 > EDGE && !TAMask_B[g1][g2*W2+g3*W3+(g4-1)] )
                                tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                        }
                    }

                    for (int i = 0; i < tmpVec.size(); i ++)  {

                        g1 = (int)(tmpVec[i] / M1);
                        g2 = (int)((tmpVec[i]%M1) / M2);
                        g3 = (int)((tmpVec[i]%M2) / M3);
                        g4 = (int)(tmpVec[i] % M3);
                        TAMask_B[g1][g2*W2+g3*W3+g4] = 1;

                        x1_min_B = (g1 < x1_min_B) ? g1 : x1_min_B;
                        x1_max_B = (g1 > x1_max_B) ? g1 : x1_max_B;
                        x2_min_B = (g2 < x2_min_B) ? g2 : x2_min_B;
                        x2_max_B = (g2 > x2_max_B) ? g2 : x2_max_B;
                        x3_min_B = (g3 < x3_min_B) ? g3 : x3_min_B;
                        x3_max_B = (g3 > x3_max_B) ? g3 : x3_max_B;
                        x4_min_B = (g4 < x4_min_B) ? g4 : x4_min_B;
                        x4_max_B = (g4 > x4_max_B) ? g4 : x4_max_B;
                    }
                    tmpVec.clear();
                    ExBD_B.clear();

                    #pragma omp parallel for reduction(merge: ExBD_B) private(g1,g2,g3,g4,n1,n2,n3,n4) schedule(runtime)
                    for (int i = 0; i < ExFF_B.size(); i++)  {
                        if (Check_B[i])  
                        {
                            g1 = (int)(ExFF_B[i] / M1);
                            g2 = (int)((ExFF_B[i]%M1) / M2);
                            g3 = (int)((ExFF_B[i]%M2) / M3);
                            g4 = (int)(ExFF_B[i] % M3);

                            ExBD_B.push_back(ExFF_B[i]);

                            for (int j = 0; j < nneigh; j ++)  {

                                n1 = neighlist[j][0];
                                n2 = neighlist[j][1];
                                n3 = neighlist[j][2];
                                n4 = neighlist[j][3];

                                if (TAMask_B[g1+n1][(g2+n2)*W2+(g3+n3)*W3+(g4+n4)])
                                    ExBD_B.push_back(GridToIdx(g1+n1,g2+n2,g3+n3,g4+n4));
                            }
                        }
                    }

                    // Find unique elements (ExBD)
                    __gnu_parallel::sort(ExBD_B.begin(),ExBD_B.end());
                    it = std::unique (ExBD_B.begin(), ExBD_B.end()); 
                    ExBD_B.resize(std::distance(ExBD_B.begin(),it));

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-cx-1: CASE 1 ExBD B) = %.4e sec\n", t_1_elapsed); 

                }  // if ExFF_B.size() > 0

                if ( ExFF_A.size() > 0 || ExFF_B.size() > 0 )  {

                    // Combine ExBD_A and ExBD_B
                    if (ExBD_A.size() != 0 && ExBD_B.size() != 0)  {
                        tmpVec.resize(ExBD_A.size() + ExBD_B.size());
                        __gnu_parallel::sort (ExBD_A.begin(), ExBD_A.end());
                        __gnu_parallel::sort (ExBD_B.begin(), ExBD_B.end());
                        it=std::set_union( ExBD_A.begin(), ExBD_A.end(), ExBD_B.begin(), ExBD_B.end(), tmpVec.begin() );
                        tmpVec.resize(it - tmpVec.begin()); 
                        tmpVec.swap(ExBD);
                        tmpVec.clear();
                    }
                    else if (ExBD_A.size() != 0 && ExBD_B.size() == 0)  {
                        ExBD_A.swap(ExBD);
                    }
                    else if (ExBD_A.size() == 0 && ExBD_B.size() != 0)  {
                        ExBD_B.swap(ExBD);
                    }

                    // Box range union
                    x1_max_C = (x1_max_A > x1_max_B) ? x1_max_A : x1_max_B;
                    x1_min_C = (x1_min_A < x1_min_B) ? x1_min_A : x1_min_B;
                    x2_max_C = (x2_max_A > x2_max_B) ? x2_max_A : x2_max_B;
                    x2_min_C = (x2_min_A < x2_min_B) ? x2_min_A : x2_min_B;
                    x3_max_C = (x3_max_A > x3_max_B) ? x3_max_A : x3_max_B;
                    x3_min_C = (x3_min_A < x3_min_B) ? x3_min_A : x3_min_B;
                    x4_max_C = (x4_max_A > x4_max_B) ? x4_max_A : x4_max_B;
                    x4_min_C = (x4_min_A < x4_min_B) ? x4_min_A : x4_min_B;

                    // Runge???Kutta 4

                    t_1_begin = omp_get_wtime();

                    // RK4-1
                    #pragma omp parallel for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,\
                                                    f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                                    f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                                    c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i = 0; i < ExBD.size(); i++)  {

                        g1 = (int)(ExBD[i] / M1);
                        g2 = (int)((ExBD[i]%M1) / M2);
                        g3 = (int)((ExBD[i]%M2) / M3);
                        g4 = (int)(ExBD[i] % M3);

                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];
                        xx3 = Box[4] + g3 * H[2];
                        xx4 = Box[6] + g4 * H[3];

                        f0_A = F_A[g1][g2*W2+g3*W3+g4];
                        f1p1_A = F_A[g1+1][g2*W2+g3*W3+g4];
                        f1m1_A = F_A[g1-1][g2*W2+g3*W3+g4];
                        f2p1_A = F_A[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_A = F_A[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_A = F_A[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_A = F_A[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_A = F_A[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_A = F_A[g1][g2*W2+g3*W3+(g4-1)];

                        f0_B = F_B[g1][g2*W2+g3*W3+g4];
                        f1p1_B = F_B[g1+1][g2*W2+g3*W3+g4];
                        f1m1_B = F_B[g1-1][g2*W2+g3*W3+g4];
                        f2p1_B = F_B[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_B = F_B[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_B = F_B[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_B = F_B[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_B = F_B[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_B = F_B[g1][g2*W2+g3*W3+(g4-1)];

                        c_f0_a = 2.0 * f0_A;
                        c_f0_b = 2.0 * f0_B;
                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                        KK1_A[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_A + f1p1_A - c_f0_a) + 
                                                    Ikh2mh1sq * (f2m1_A + f2p1_A - c_f0_a) +
                                                    Ikh2mh2sq * (f3m1_A + f3p1_A - c_f0_a) + 
                                                    Ikh2mh3sq * (f4m1_A + f4p1_A - c_f0_a) - 
                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                        FF_A[g1][g2*W2+g3*W3+g4] = F_A[g1][g2*W2+g3*W3+g4] + KK1_A[g1][g2*W2+g3*W3+g4] / 6.0;

                        KK1_B[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_B + f1p1_B - c_f0_b) + 
                                                     Ikh2mh1sq * (f2m1_B + f2p1_B - c_f0_b) +
                                                     Ikh2mh2sq * (f3m1_B + f3p1_B - c_f0_b) + 
                                                     Ikh2mh3sq * (f4m1_B + f4p1_B - c_f0_b) - 
                                                     Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                        FF_B[g1][g2*W2+g3*W3+g4] = F_B[g1][g2*W2+g3*W3+g4] + KK1_B[g1][g2*W2+g3*W3+g4] / 6.0;
                    }

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-11: CASE 1 KK1) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();

                    // RK4-2
                    #pragma omp for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i = 0; i < ExBD.size(); i++)  {

                        g1 = (int)(ExBD[i] / M1);
                        g2 = (int)((ExBD[i]%M1) / M2);
                        g3 = (int)((ExBD[i]%M2) / M3);
                        g4 = (int)(ExBD[i] % M3);

                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];
                        xx3 = Box[4] + g3 * H[2];
                        xx4 = Box[6] + g4 * H[3];                  

                        f0_A = F_A[g1][g2*W2+g3*W3+g4];
                        f1p1_A = F_A[g1+1][g2*W2+g3*W3+g4];
                        f1m1_A = F_A[g1-1][g2*W2+g3*W3+g4];
                        f2p1_A = F_A[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_A = F_A[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_A = F_A[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_A = F_A[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_A = F_A[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_A = F_A[g1][g2*W2+g3*W3+(g4-1)];
                        kk0_A = KK1_A[g1][g2*W2+g3*W3+g4];
                        kk1p1_A = KK1_A[g1+1][g2*W2+g3*W3+g4];
                        kk1m1_A = KK1_A[g1-1][g2*W2+g3*W3+g4];
                        kk2p1_A = KK1_A[g1][(g2+1)*W2+g3*W3+g4];
                        kk2m1_A = KK1_A[g1][(g2-1)*W2+g3*W3+g4];
                        kk3p1_A = KK1_A[g1][g2*W2+(g3+1)*W3+g4];
                        kk3m1_A = KK1_A[g1][g2*W2+(g3-1)*W3+g4];
                        kk4p1_A = KK1_A[g1][g2*W2+g3*W3+(g4+1)];
                        kk4m1_A = KK1_A[g1][g2*W2+g3*W3+(g4-1)];

                        f0_B = F_B[g1][g2*W2+g3*W3+g4];
                        f1p1_B = F_B[g1+1][g2*W2+g3*W3+g4];
                        f1m1_B = F_B[g1-1][g2*W2+g3*W3+g4];
                        f2p1_B = F_B[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_B = F_B[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_B = F_B[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_B = F_B[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_B = F_B[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_B = F_B[g1][g2*W2+g3*W3+(g4-1)];
                        kk0_B = KK1_B[g1][g2*W2+g3*W3+g4];
                        kk1p1_B = KK1_B[g1+1][g2*W2+g3*W3+g4];
                        kk1m1_B = KK1_B[g1-1][g2*W2+g3*W3+g4];
                        kk2p1_B = KK1_B[g1][(g2+1)*W2+g3*W3+g4];
                        kk2m1_B = KK1_B[g1][(g2-1)*W2+g3*W3+g4];
                        kk3p1_B = KK1_B[g1][g2*W2+(g3+1)*W3+g4];
                        kk3m1_B = KK1_B[g1][g2*W2+(g3-1)*W3+g4];
                        kk4p1_B = KK1_B[g1][g2*W2+g3*W3+(g4+1)];
                        kk4m1_B = KK1_B[g1][g2*W2+g3*W3+(g4-1)];

                        c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                        c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                        KK2_A[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                    Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) + 
                                                    Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                    Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                        FF_A[g1][g2*W2+g3*W3+g4] += KK2_A[g1][g2*W2+g3*W3+g4] / 3.0;

                        KK2_B[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                    Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                    Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                    Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                    Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                        FF_B[g1][g2*W2+g3*W3+g4] += KK2_B[g1][g2*W2+g3*W3+g4] / 3.0;
                    }

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-12: CASE 1 KK2) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();

                    // RK4-3
                    #pragma omp for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i = 0; i < ExBD.size(); i++)  {

                        g1 = (int)(ExBD[i] / M1);
                        g2 = (int)((ExBD[i]%M1) / M2);
                        g3 = (int)((ExBD[i]%M2) / M3);
                        g4 = (int)(ExBD[i] % M3);

                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];
                        xx3 = Box[4] + g3 * H[2];
                        xx4 = Box[6] + g4 * H[3];                  

                        f0_A = F_A[g1][g2*W2+g3*W3+g4];
                        f1p1_A = F_A[g1+1][g2*W2+g3*W3+g4];
                        f1m1_A = F_A[g1-1][g2*W2+g3*W3+g4];
                        f2p1_A = F_A[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_A = F_A[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_A = F_A[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_A = F_A[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_A = F_A[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_A = F_A[g1][g2*W2+g3*W3+(g4-1)];
                        kk0_A = KK2_A[g1][g2*W2+g3*W3+g4];
                        kk1p1_A = KK2_A[g1+1][g2*W2+g3*W3+g4];
                        kk1m1_A = KK2_A[g1-1][g2*W2+g3*W3+g4];
                        kk2p1_A = KK2_A[g1][(g2+1)*W2+g3*W3+g4];
                        kk2m1_A = KK2_A[g1][(g2-1)*W2+g3*W3+g4];
                        kk3p1_A = KK2_A[g1][g2*W2+(g3+1)*W3+g4];
                        kk3m1_A = KK2_A[g1][g2*W2+(g3-1)*W3+g4];
                        kk4p1_A = KK2_A[g1][g2*W2+g3*W3+(g4+1)];
                        kk4m1_A = KK2_A[g1][g2*W2+g3*W3+(g4-1)];

                        f0_B = F_B[g1][g2*W2+g3*W3+g4];
                        f1p1_B = F_B[g1+1][g2*W2+g3*W3+g4];
                        f1m1_B = F_B[g1-1][g2*W2+g3*W3+g4];
                        f2p1_B = F_B[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_B = F_B[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_B = F_B[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_B = F_B[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_B = F_B[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_B = F_B[g1][g2*W2+g3*W3+(g4-1)];
                        kk0_B = KK2_B[g1][g2*W2+g3*W3+g4];
                        kk1p1_B = KK2_B[g1+1][g2*W2+g3*W3+g4];
                        kk1m1_B = KK2_B[g1-1][g2*W2+g3*W3+g4];
                        kk2p1_B = KK2_B[g1][(g2+1)*W2+g3*W3+g4];
                        kk2m1_B = KK2_B[g1][(g2-1)*W2+g3*W3+g4];
                        kk3p1_B = KK2_B[g1][g2*W2+(g3+1)*W3+g4];
                        kk3m1_B = KK2_B[g1][g2*W2+(g3-1)*W3+g4];
                        kk4p1_B = KK2_B[g1][g2*W2+g3*W3+(g4+1)];
                        kk4m1_B = KK2_B[g1][g2*W2+g3*W3+(g4-1)];

                        c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                        c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                        KK3_A[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                    Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) +
                                                    Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                    Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                        FF_A[g1][g2*W2+g3*W3+g4] += KK3_A[g1][g2*W2+g3*W3+g4] / 3.0;

                        KK3_B[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                    Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                    Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                    Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                    Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                        FF_B[g1][g2*W2+g3*W3+g4] += KK3_B[g1][g2*W2+g3*W3+g4] / 3.0;
                    }

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-13: CASE 1 KK3) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();

                    // RK4-4
                    #pragma omp for private(g1,g2,g3,g4,xx1,xx2,xx3,xx4,\
                                            f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                            kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                            f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                            kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                            c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                    for (int i = 0; i < ExBD.size(); i++)  {

                        g1 = (int)(ExBD[i] / M1);
                        g2 = (int)((ExBD[i]%M1) / M2);
                        g3 = (int)((ExBD[i]%M2) / M3);
                        g4 = (int)(ExBD[i] % M3);

                        xx1 = Box[0] + g1 * H[0];
                        xx2 = Box[2] + g2 * H[1];
                        xx3 = Box[4] + g3 * H[2];
                        xx4 = Box[6] + g4 * H[3];                  

                        f0_A = F_A[g1][g2*W2+g3*W3+g4];
                        f1p1_A = F_A[g1+1][g2*W2+g3*W3+g4];
                        f1m1_A = F_A[g1-1][g2*W2+g3*W3+g4];
                        f2p1_A = F_A[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_A = F_A[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_A = F_A[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_A = F_A[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_A = F_A[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_A = F_A[g1][g2*W2+g3*W3+(g4-1)];
                        kk0_A = KK3_A[g1][g2*W2+g3*W3+g4];
                        kk1p1_A = KK3_A[g1+1][g2*W2+g3*W3+g4];
                        kk1m1_A = KK3_A[g1-1][g2*W2+g3*W3+g4];
                        kk2p1_A = KK3_A[g1][(g2+1)*W2+g3*W3+g4];
                        kk2m1_A = KK3_A[g1][(g2-1)*W2+g3*W3+g4];
                        kk3p1_A = KK3_A[g1][g2*W2+(g3+1)*W3+g4];
                        kk3m1_A = KK3_A[g1][g2*W2+(g3-1)*W3+g4];
                        kk4p1_A = KK3_A[g1][g2*W2+g3*W3+(g4+1)];
                        kk4m1_A = KK3_A[g1][g2*W2+g3*W3+(g4-1)];

                        f0_B = F_B[g1][g2*W2+g3*W3+g4];
                        f1p1_B = F_B[g1+1][g2*W2+g3*W3+g4];
                        f1m1_B = F_B[g1-1][g2*W2+g3*W3+g4];
                        f2p1_B = F_B[g1][(g2+1)*W2+g3*W3+g4];
                        f2m1_B = F_B[g1][(g2-1)*W2+g3*W3+g4];
                        f3p1_B = F_B[g1][g2*W2+(g3+1)*W3+g4];
                        f3m1_B = F_B[g1][g2*W2+(g3-1)*W3+g4];
                        f4p1_B = F_B[g1][g2*W2+g3*W3+(g4+1)];
                        f4m1_B = F_B[g1][g2*W2+g3*W3+(g4-1)];
                        kk0_B = KK3_B[g1][g2*W2+g3*W3+g4];
                        kk1p1_B = KK3_B[g1+1][g2*W2+g3*W3+g4];
                        kk1m1_B = KK3_B[g1-1][g2*W2+g3*W3+g4];
                        kk2p1_B = KK3_B[g1][(g2+1)*W2+g3*W3+g4];
                        kk2m1_B = KK3_B[g1][(g2-1)*W2+g3*W3+g4];
                        kk3p1_B = KK3_B[g1][g2*W2+(g3+1)*W3+g4];
                        kk3m1_B = KK3_B[g1][g2*W2+(g3-1)*W3+g4];
                        kk4p1_B = KK3_B[g1][g2*W2+g3*W3+(g4+1)];
                        kk4m1_B = KK3_B[g1][g2*W2+g3*W3+(g4-1)];

                        c_f0_a = 2.0 * (f0_A+kk0_A);
                        c_f0_b = 2.0 * (f0_B+kk0_B);
                        pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                        pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                        pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                        KK4_A[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_A+kk1m1_A + f1p1_A+kk1p1_A - c_f0_a) + 
                                                    Ikh2mh1sq * (f2m1_A+kk2m1_A + f2p1_A+kk2p1_A - c_f0_a) + 
                                                    Ikh2mh2sq * (f3m1_A+kk3m1_A + f3p1_A+kk3p1_A - c_f0_a) +
                                                    Ikh2mh3sq * (f4m1_A+kk4m1_A + f4p1_A+kk4p1_A - c_f0_a) - 
                                                    Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                        FF_A[g1][g2*W2+g3*W3+g4] += KK4_A[g1][g2*W2+g3*W3+g4] / 6.0;

                        KK4_B[g1][g2*W2+g3*W3+g4] = Ikh2mh0sq * (f1m1_B+kk1m1_B + f1p1_B+kk1p1_B - c_f0_b) + 
                                                    Ikh2mh1sq * (f2m1_B+kk2m1_B + f2p1_B+kk2p1_B - c_f0_b) + 
                                                    Ikh2mh2sq * (f3m1_B+kk3m1_B + f3p1_B+kk3p1_B - c_f0_b) + 
                                                    Ikh2mh3sq * (f4m1_B+kk4m1_B + f4p1_B+kk4p1_B - c_f0_b) - 
                                                    Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                        FF_B[g1][g2*W2+g3*W3+g4] += KK4_B[g1][g2*W2+g3*W3+g4] / 6.0;
                    }

                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kkx-14: CASE 1 KK4) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();

                }  // if ( ExFF_A.size() > 0 || ExFF_B.size() > 0 )
            }  // if not ( isFirstExtrp_A && isFirstExtrp_B ) and (ExFF_A.size() == 1 || ExFF_B.size() == 1)

            // Check Multiple Expanding 
            // TBL = index of FF that FF(TBL) is higher than TolL

            /* ---------------- Surface A ---------------- */

            if ( ExFF_A.size() > 0 )  {

                t_1_begin = omp_get_wtime();

                TBL_A.clear();
                tmpVec.clear();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4,b1,b2,b3,b4,nx1,nx2,nx3,nx4,\
                                                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A) schedule(runtime)
                for (int i = 0; i < ExFF_A.size(); i++)
                {
                    if (Check_A[i])  {

                        g1 = (int)(ExFF_A[i] / M1);
                        g2 = (int)((ExFF_A[i]%M1) / M2);
                        g3 = (int)((ExFF_A[i]%M2) / M3);
                        g4 = (int)(ExFF_A[i] % M3);
                        
                        nx1 = int(TAMask_A[g1+1][g2*W2+g3*W3+g4]) + int(TAMask_A[g1-1][g2*W2+g3*W3+g4]);
                        nx2 = int(TAMask_A[g1][(g2+1)*W2+g3*W3+g4]) + int(TAMask_A[g1][(g2-1)*W2+g3*W3+g4]);
                        nx3 = int(TAMask_A[g1][g2*W2+(g3+1)*W3+g4]) + int(TAMask_A[g1][g2*W2+(g3-1)*W3+g4]);
                        nx4 = int(TAMask_A[g1][g2*W2+g3*W3+(g4+1)]) + int(TAMask_A[g1][g2*W2+g3*W3+(g4-1)]);

                        f1p1_A = (TAMask_A[g1+1][g2*W2+g3*W3+g4]) ? FF_A[g1+1][g2*W2+g3*W3+g4] : FF_A[g1][g2*W2+g3*W3+g4];
                        f1m1_A = (TAMask_A[g1-1][g2*W2+g3*W3+g4]) ? FF_A[g1-1][g2*W2+g3*W3+g4] : FF_A[g1][g2*W2+g3*W3+g4];
                        f2p1_A = (TAMask_A[g1][(g2+1)*W2+g3*W3+g4]) ? FF_A[g1][(g2+1)*W2+g3*W3+g4] : FF_A[g1][g2*W2+g3*W3+g4];
                        f2m1_A = (TAMask_A[g1][(g2-1)*W2+g3*W3+g4]) ? FF_A[g1][(g2-1)*W2+g3*W3+g4] : FF_A[g1][g2*W2+g3*W3+g4];
                        f3p1_A = (TAMask_A[g1][g2*W2+(g3+1)*W3+g4]) ? FF_A[g1][g2*W2+(g3+1)*W3+g4] : FF_A[g1][g2*W2+g3*W3+g4];
                        f3m1_A = (TAMask_A[g1][g2*W2+(g3-1)*W3+g4]) ? FF_A[g1][g2*W2+(g3-1)*W3+g4] : FF_A[g1][g2*W2+g3*W3+g4];
                        f4p1_A = (TAMask_A[g1][g2*W2+g3*W3+(g4+1)]) ? FF_A[g1][g2*W2+g3*W3+(g4+1)] : FF_A[g1][g2*W2+g3*W3+g4];
                        f4m1_A = (TAMask_A[g1][g2*W2+g3*W3+(g4-1)]) ? FF_A[g1][g2*W2+g3*W3+(g4-1)] : FF_A[g1][g2*W2+g3*W3+g4];

                        f0_A = FF_A[g1][g2*W2+g3*W3+g4];

                        b1 = std::abs(f0_A * std::conj(f0_A)) >= TolH_A;
                        b2 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_A - f1m1_A)/(nx1*H[0]),2)) + \
                             ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_A - f2m1_A)/(nx2*H[1]),2)) + \
                             ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_A - f3m1_A)/(nx3*H[2]),2)) + \
                             ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_A - f4m1_A)/(nx4*H[3]),2)) >= TolLd_A_sq;

                        b3 = g1 > EDGE && g2 > EDGE && g3 > EDGE && g4 > EDGE;
                        b4 = (g1 < BoxShape[0]-EDGE-1) && (g2 < BoxShape[1]-EDGE-1) && (g3 < BoxShape[2]-EDGE-1) && (g4 < BoxShape[3]-EDGE-1);

                        if ((b1||b2) && b3 && b4)  {
                            tmpVec.push_back(ExFF_A[i]);
                        }
                    }
                }        
                tmpVec.swap(TBL_A);
                tmpVec.clear(); 

                // TBL & TBL_P set difference
                tmpVec.resize(TBL_P_A.size() + TBL_A.size());
                __gnu_parallel::sort (TBL_A.begin(), TBL_A.end());
                __gnu_parallel::sort (TBL_P_A.begin(), TBL_P_A.end());
                it=std::set_difference( TBL_A.begin(), TBL_A.end(), TBL_P_A.begin(), TBL_P_A.end(), tmpVec.begin() );
                tmpVec.resize(it - tmpVec.begin()); 
                tmpVec.swap(TBL_A);
                tmpVec.clear();

                // Combine TBL and TBL_P
                TBL_P_A.reserve(TBL_P_A.size() + TBL_A.size());
                TBL_P_A.insert(TBL_P_A.end(), TBL_A.begin(), TBL_A.end());

                // Find unique elements
                __gnu_parallel::sort(TBL_P_A.begin(),TBL_P_A.end());
                it = std::unique (TBL_P_A.begin(), TBL_P_A.end()); 
                TBL_P_A.resize(std::distance(TBL_P_A.begin(),it));

                // Update isFirstExtrp
                Excount_A += 1;
                if (Excount_A == ExLimit) TBL_A.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-c-3 CASE 1 TBL A) = %.4e sec\n", t_1_elapsed); 
            }

            /* ---------------- Surface B ---------------- */

            if ( ExFF_B.size() > 0 )  {

                t_1_begin = omp_get_wtime();

                TBL_B.clear();
                tmpVec.clear();

                #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4,b1,b2,b3,b4,nx1,nx2,nx3,nx4,\
                                                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B) schedule(runtime)
                for (int i = 0; i < ExFF_B.size(); i++)
                {
                    if (Check_B[i])  {

                        g1 = (int)(ExFF_B[i] / M1);
                        g2 = (int)((ExFF_B[i]%M1) / M2);
                        g3 = (int)((ExFF_B[i]%M2) / M3);
                        g4 = (int)(ExFF_B[i] % M3);
                        
                        nx1 = int(TAMask_B[g1+1][g2*W2+g3*W3+g4]) + int(TAMask_B[g1-1][g2*W2+g3*W3+g4]);
                        nx2 = int(TAMask_B[g1][(g2+1)*W2+g3*W3+g4]) + int(TAMask_B[g1][(g2-1)*W2+g3*W3+g4]);
                        nx3 = int(TAMask_B[g1][g2*W2+(g3+1)*W3+g4]) + int(TAMask_B[g1][g2*W2+(g3-1)*W3+g4]);
                        nx4 = int(TAMask_B[g1][g2*W2+g3*W3+(g4+1)]) + int(TAMask_B[g1][g2*W2+g3*W3+(g4-1)]);

                        f1p1_B = (TAMask_B[g1+1][g2*W2+g3*W3+g4]) ? FF_B[g1+1][g2*W2+g3*W3+g4] : FF_B[g1][g2*W2+g3*W3+g4];
                        f1m1_B = (TAMask_B[g1-1][g2*W2+g3*W3+g4]) ? FF_B[g1-1][g2*W2+g3*W3+g4] : FF_B[g1][g2*W2+g3*W3+g4];
                        f2p1_B = (TAMask_B[g1][(g2+1)*W2+g3*W3+g4]) ? FF_B[g1][(g2+1)*W2+g3*W3+g4] : FF_B[g1][g2*W2+g3*W3+g4];
                        f2m1_B = (TAMask_B[g1][(g2-1)*W2+g3*W3+g4]) ? FF_B[g1][(g2-1)*W2+g3*W3+g4] : FF_B[g1][g2*W2+g3*W3+g4];
                        f3p1_B = (TAMask_B[g1][g2*W2+(g3+1)*W3+g4]) ? FF_B[g1][g2*W2+(g3+1)*W3+g4] : FF_B[g1][g2*W2+g3*W3+g4];
                        f3m1_B = (TAMask_B[g1][g2*W2+(g3-1)*W3+g4]) ? FF_B[g1][g2*W2+(g3-1)*W3+g4] : FF_B[g1][g2*W2+g3*W3+g4];
                        f4p1_B = (TAMask_B[g1][g2*W2+g3*W3+(g4+1)]) ? FF_B[g1][g2*W2+g3*W3+(g4+1)] : FF_B[g1][g2*W2+g3*W3+g4];
                        f4m1_B = (TAMask_B[g1][g2*W2+g3*W3+(g4-1)]) ? FF_B[g1][g2*W2+g3*W3+(g4-1)] : FF_B[g1][g2*W2+g3*W3+g4];

                        f0_B = FF_B[g1][g2*W2+g3*W3+g4];

                        b1 = std::abs(f0_B * std::conj(f0_B)) >= TolH_B;
                        b2 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_B - f1m1_B)/(nx1*H[0]),2)) + \
                             ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_B - f2m1_B)/(nx2*H[1]),2)) + \
                             ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_B - f3m1_B)/(nx3*H[2]),2)) + \
                             ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_B - f4m1_B)/(nx4*H[3]),2)) >= TolLd_B_sq;

                        b3 = g1 > EDGE && g2 > EDGE && g3 > EDGE && g4 > EDGE;
                        b4 = (g1 < BoxShape[0]-EDGE-1) && (g2 < BoxShape[1]-EDGE-1) && (g3 < BoxShape[2]-EDGE-1) && (g4 < BoxShape[3]-EDGE-1);

                        if ((b1||b2) && b3 && b4)  {
                            tmpVec.push_back(ExFF_B[i]);
                        }
                    }
                }         
                tmpVec.swap(TBL_B);
                tmpVec.clear(); 

                // TBL & TBL_P set difference
                tmpVec.resize(TBL_P_B.size() + TBL_B.size());
                __gnu_parallel::sort (TBL_B.begin(), TBL_B.end());
                __gnu_parallel::sort (TBL_P_B.begin(), TBL_P_B.end());
                it=std::set_difference( TBL_B.begin(), TBL_B.end(), TBL_P_B.begin(), TBL_P_B.end(), tmpVec.begin() );
                tmpVec.resize(it - tmpVec.begin()); 
                tmpVec.swap(TBL_B);
                tmpVec.clear();

                // Combine TBL and TBL_P
                TBL_P_B.reserve(TBL_P_B.size() + TBL_B.size());
                TBL_P_B.insert(TBL_P_B.end(), TBL_B.begin(), TBL_B.end());

                // Find unique elements
                __gnu_parallel::sort(TBL_P_B.begin(),TBL_P_B.end());
                it = std::unique (TBL_P_B.begin(), TBL_P_B.end()); 
                TBL_P_B.resize(std::distance(TBL_P_B.begin(),it));

                // Update isFirstExtrp
                Excount_B += 1;
                if (Excount_B == ExLimit) TBL_B.clear();

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin;
                t_overhead += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-c-3 CASE 1 TBL B) = %.4e sec\n", t_1_elapsed); 
            }

        }  // while ( !isFullGrid && (TBL_A.size() != 0 || TBL_B.size() != 0) && (Excount_A < ExLimit && Excount_B < ExLimit) )
        
        // .........................................................................................

        // CASE 2: Truncating without extrapolation

        if ( !isExtrapolate_A && !isExtrapolate_B && !isFullGrid )
        {
            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }

                // RK4-1
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                    for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                        for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                            for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                    xx1 = Box[0] + i1 * H[0];
                                    xx2 = Box[2] + i2 * H[1];
                                    xx3 = Box[4] + i3 * H[2];
                                    xx4 = Box[6] + i4 * H[3];

                                    f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                    f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                    f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                    f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];

                                    f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                    f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                    f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                    f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];

                                    c_f0_a = 2.0 * f0_A;
                                    c_f0_b = 2.0 * f0_B;
                                    pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                    pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                    pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                    KK1_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A + f1p1_A - c_f0_a) + 
                                                                Ikh2mh1sq * (f2m1_A + f2p1_A - c_f0_a) +
                                                                Ikh2mh2sq * (f3m1_A + f3p1_A - c_f0_a) + 
                                                                Ikh2mh3sq * (f4m1_A + f4p1_A - c_f0_a) - 
                                                                Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                    FF_A[i1][i2*W2+i3*W3+i4] = F_A[i1][i2*W2+i3*W3+i4] + KK1_A[i1][i2*W2+i3*W3+i4] / 6.0;

                                    KK1_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B + f1p1_B - c_f0_b) + 
                                                                Ikh2mh1sq * (f2m1_B + f2p1_B - c_f0_b) +
                                                                Ikh2mh2sq * (f3m1_B + f3p1_B - c_f0_b) + 
                                                                Ikh2mh3sq * (f4m1_B + f4p1_B - c_f0_b) - 
                                                                Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                    FF_B[i1][i2*W2+i3*W3+i4] = F_B[i1][i2*W2+i3*W3+i4] + KK1_B[i1][i2*W2+i3*W3+i4] / 6.0;
                                }
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-21: CASE 2 KK1) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-2
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                    for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                        for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                            for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                    xx1 = Box[0] + i1 * H[0];
                                    xx2 = Box[2] + i2 * H[1];
                                    xx3 = Box[4] + i3 * H[2];
                                    xx4 = Box[6] + i4 * H[3];

                                    f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                    f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                    f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                    f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                    kk0_A = KK1_A[i1][i2*W2+i3*W3+i4];
                                    kk1p1_A = KK1_A[i1+1][i2*W2+i3*W3+i4];
                                    kk1m1_A = KK1_A[i1-1][i2*W2+i3*W3+i4];
                                    kk2p1_A = KK1_A[i1][(i2+1)*W2+i3*W3+i4];
                                    kk2m1_A = KK1_A[i1][(i2-1)*W2+i3*W3+i4];
                                    kk3p1_A = KK1_A[i1][i2*W2+(i3+1)*W3+i4];
                                    kk3m1_A = KK1_A[i1][i2*W2+(i3-1)*W3+i4];
                                    kk4p1_A = KK1_A[i1][i2*W2+i3*W3+(i4+1)];
                                    kk4m1_A = KK1_A[i1][i2*W2+i3*W3+(i4-1)];

                                    f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                    f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                    f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                    f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                    kk0_B = KK1_B[i1][i2*W2+i3*W3+i4];
                                    kk1p1_B = KK1_B[i1+1][i2*W2+i3*W3+i4];
                                    kk1m1_B = KK1_B[i1-1][i2*W2+i3*W3+i4];
                                    kk2p1_B = KK1_B[i1][(i2+1)*W2+i3*W3+i4];
                                    kk2m1_B = KK1_B[i1][(i2-1)*W2+i3*W3+i4];
                                    kk3p1_B = KK1_B[i1][i2*W2+(i3+1)*W3+i4];
                                    kk3m1_B = KK1_B[i1][i2*W2+(i3-1)*W3+i4];
                                    kk4p1_B = KK1_B[i1][i2*W2+i3*W3+(i4+1)];
                                    kk4m1_B = KK1_B[i1][i2*W2+i3*W3+(i4-1)];

                                    c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                                    c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                                    pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                    pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                    pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                    KK2_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                                Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) + 
                                                                Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                                Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                                Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                    FF_A[i1][i2*W2+i3*W3+i4] += KK2_A[i1][i2*W2+i3*W3+i4] / 3.0;

                                    KK2_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                                Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                                Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                                Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                                Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                    FF_B[i1][i2*W2+i3*W3+i4] += KK2_B[i1][i2*W2+i3*W3+i4] / 3.0;
                                }
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-22: CASE 2 KK2) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                    for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                        for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                            for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                    xx1 = Box[0] + i1 * H[0];
                                    xx2 = Box[2] + i2 * H[1];
                                    xx3 = Box[4] + i3 * H[2];
                                    xx4 = Box[6] + i4 * H[3];

                                    f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                    f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                    f3m1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                    f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                    kk0_A = KK2_A[i1][i2*W2+i3*W3+i4];
                                    kk1p1_A = KK2_A[i1+1][i2*W2+i3*W3+i4];
                                    kk1m1_A = KK2_A[i1-1][i2*W2+i3*W3+i4];
                                    kk2p1_A = KK2_A[i1][(i2+1)*W2+i3*W3+i4];
                                    kk2m1_A = KK2_A[i1][(i2-1)*W2+i3*W3+i4];
                                    kk3p1_A = KK2_A[i1][i2*W2+(i3+1)*W3+i4];
                                    kk3m1_A = KK2_A[i1][i2*W2+(i3-1)*W3+i4];
                                    kk4p1_A = KK2_A[i1][i2*W2+i3*W3+(i4+1)];
                                    kk4m1_A = KK2_A[i1][i2*W2+i3*W3+(i4-1)];

                                    f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                    f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                    f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                    f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                    kk0_B = KK2_B[i1][i2*W2+i3*W3+i4];
                                    kk1p1_B = KK2_B[i1+1][i2*W2+i3*W3+i4];
                                    kk1m1_B = KK2_B[i1-1][i2*W2+i3*W3+i4];
                                    kk2p1_B = KK2_B[i1][(i2+1)*W2+i3*W3+i4];
                                    kk2m1_B = KK2_B[i1][(i2-1)*W2+i3*W3+i4];
                                    kk3p1_B = KK2_B[i1][i2*W2+(i3+1)*W3+i4];
                                    kk3m1_B = KK2_B[i1][i2*W2+(i3-1)*W3+i4];
                                    kk4p1_B = KK2_B[i1][i2*W2+i3*W3+(i4+1)];
                                    kk4m1_B = KK2_B[i1][i2*W2+i3*W3+(i4-1)];

                                    c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                                    c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                                    pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                    pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                    pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                    KK3_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                            Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) +
                                                            Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                            Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                            Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                    FF_A[i1][i2*W2+i3*W3+i4] += KK3_A[i1][i2*W2+i3*W3+i4] / 3.0;

                                    KK3_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                            Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                            Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                            Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                            Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                    FF_B[i1][i2*W2+i3*W3+i4] += KK3_B[i1][i2*W2+i3*W3+i4] / 3.0;
                                }
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-23: CASE 2 KK3) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                    for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                        for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                            for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                    xx1 = Box[0] + i1 * H[0];
                                    xx2 = Box[2] + i2 * H[1];
                                    xx3 = Box[4] + i3 * H[2];
                                    xx4 = Box[6] + i4 * H[3];

                                    f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                    f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                    f3m1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                    f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                    kk0_A = KK3_A[i1][i2*W2+i3*W3+i4];
                                    kk1p1_A = KK3_A[i1+1][i2*W2+i3*W3+i4];
                                    kk1m1_A = KK3_A[i1-1][i2*W2+i3*W3+i4];
                                    kk2p1_A = KK3_A[i1][(i2+1)*W2+i3*W3+i4];
                                    kk2m1_A = KK3_A[i1][(i2-1)*W2+i3*W3+i4];
                                    kk3p1_A = KK3_A[i1][i2*W2+(i3+1)*W3+i4];
                                    kk3m1_A = KK3_A[i1][i2*W2+(i3-1)*W3+i4];
                                    kk4p1_A = KK3_A[i1][i2*W2+i3*W3+(i4+1)];
                                    kk4m1_A = KK3_A[i1][i2*W2+i3*W3+(i4-1)];

                                    f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                    f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                    f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                    f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                    f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                    f3p1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                    f3m1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                    f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                    f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                    kk0_B = KK3_B[i1][i2*W2+i3*W3+i4];
                                    kk1p1_B = KK3_B[i1+1][i2*W2+i3*W3+i4];
                                    kk1m1_B = KK3_B[i1-1][i2*W2+i3*W3+i4];
                                    kk2p1_B = KK3_B[i1][(i2+1)*W2+i3*W3+i4];
                                    kk2m1_B = KK3_B[i1][(i2-1)*W2+i3*W3+i4];
                                    kk3p1_B = KK3_B[i1][i2*W2+(i3+1)*W3+i4];
                                    kk3m1_B = KK3_B[i1][i2*W2+(i3-1)*W3+i4];
                                    kk4p1_B = KK3_B[i1][i2*W2+i3*W3+(i4+1)];
                                    kk4m1_B = KK3_B[i1][i2*W2+i3*W3+(i4-1)];

                                    c_f0_a = 2.0 * (f0_A+kk0_A);
                                    c_f0_b = 2.0 * (f0_B+kk0_B);
                                    pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                    pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                    pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                    KK4_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+kk1m1_A + f1p1_A+kk1p1_A - c_f0_a) + 
                                                                Ikh2mh1sq * (f2m1_A+kk2m1_A + f2p1_A+kk2p1_A - c_f0_a) + 
                                                                Ikh2mh2sq * (f3m1_A+kk3m1_A + f3p1_A+kk3p1_A - c_f0_a) +
                                                                Ikh2mh3sq * (f4m1_A+kk4m1_A + f4p1_A+kk4p1_A - c_f0_a) - 
                                                                Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                    FF_A[i1][i2*W2+i3*W3+i4] += KK4_A[i1][i2*W2+i3*W3+i4] / 6.0;

                                    KK4_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+kk1m1_B + f1p1_B+kk1p1_B - c_f0_b) + 
                                                                Ikh2mh1sq * (f2m1_B+kk2m1_B + f2p1_B+kk2p1_B - c_f0_b) + 
                                                                Ikh2mh2sq * (f3m1_B+kk3m1_B + f3p1_B+kk3p1_B - c_f0_b) + 
                                                                Ikh2mh3sq * (f4m1_B+kk4m1_B + f4p1_B+kk4p1_B - c_f0_b) - 
                                                                Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                    FF_B[i1][i2*W2+i3*W3+i4] += KK4_B[i1][i2*W2+i3*W3+i4] / 6.0;
                                }
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_overhead += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-24: CASE 2 KK4) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            } // OMP PARALLEL
        } 
        else if ( !isExtrapolate_A && !isExtrapolate_B && isFullGrid )
        {
            // .........................................................................................

            // CASE 3: Full grid

            #pragma omp parallel
            {
                #pragma omp single nowait
                {
                    t_1_begin = omp_get_wtime();
                }
                // RK4-1
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                        for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                            for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];

                                f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                f3m1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];

                                f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                f3m1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];

                                c_f0_a = 2.0 * f0_A;
                                c_f0_b = 2.0 * f0_B;
                                pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                KK1_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A + f1p1_A - c_f0_a) + 
                                                            Ikh2mh1sq * (f2m1_A + f2p1_A - c_f0_a) +
                                                            Ikh2mh2sq * (f3m1_A + f3p1_A - c_f0_a) + 
                                                            Ikh2mh3sq * (f4m1_A + f4p1_A - c_f0_a) - 
                                                            Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                FF_A[i1][i2*W2+i3*W3+i4] = F_A[i1][i2*W2+i3*W3+i4] + KK1_A[i1][i2*W2+i3*W3+i4] / 6.0;

                                KK1_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B + f1p1_B - c_f0_b) + 
                                                            Ikh2mh1sq * (f2m1_B + f2p1_B - c_f0_b) +
                                                            Ikh2mh2sq * (f3m1_B + f3p1_B - c_f0_b) + 
                                                            Ikh2mh3sq * (f4m1_B + f4p1_B - c_f0_b) - 
                                                            Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                FF_B[i1][i2*W2+i3*W3+i4] = F_B[i1][i2*W2+i3*W3+i4] + KK1_B[i1][i2*W2+i3*W3+i4] / 6.0;
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-31: CASE 3 KK1) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-2
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                        for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                            for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];

                                f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                f3m1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                kk0_A = KK1_A[i1][i2*W2+i3*W3+i4];
                                kk1p1_A = KK1_A[i1+1][i2*W2+i3*W3+i4];
                                kk1m1_A = KK1_A[i1-1][i2*W2+i3*W3+i4];
                                kk2p1_A = KK1_A[i1][(i2+1)*W2+i3*W3+i4];
                                kk2m1_A = KK1_A[i1][(i2-1)*W2+i3*W3+i4];
                                kk3p1_A = KK1_A[i1][i2*W2+(i3+1)*W3+i4];
                                kk3m1_A = KK1_A[i1][i2*W2+(i3-1)*W3+i4];
                                kk4p1_A = KK1_A[i1][i2*W2+i3*W3+(i4+1)];
                                kk4m1_A = KK1_A[i1][i2*W2+i3*W3+(i4-1)];

                                f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                f3m1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                kk0_B = KK1_B[i1][i2*W2+i3*W3+i4];
                                kk1p1_B = KK1_B[i1+1][i2*W2+i3*W3+i4];
                                kk1m1_B = KK1_B[i1-1][i2*W2+i3*W3+i4];
                                kk2p1_B = KK1_B[i1][(i2+1)*W2+i3*W3+i4];
                                kk2m1_B = KK1_B[i1][(i2-1)*W2+i3*W3+i4];
                                kk3p1_B = KK1_B[i1][i2*W2+(i3+1)*W3+i4];
                                kk3m1_B = KK1_B[i1][i2*W2+(i3-1)*W3+i4];
                                kk4p1_B = KK1_B[i1][i2*W2+i3*W3+(i4+1)];
                                kk4m1_B = KK1_B[i1][i2*W2+i3*W3+(i4-1)];

                                c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                                c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                                pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                KK2_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                            Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) + 
                                                            Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                            Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                            Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                FF_A[i1][i2*W2+i3*W3+i4] += KK2_A[i1][i2*W2+i3*W3+i4] / 3.0;

                                KK2_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                            Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                            Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                            Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                            Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                FF_B[i1][i2*W2+i3*W3+i4] += KK2_B[i1][i2*W2+i3*W3+i4] / 3.0;
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-32: CASE 3 KK2) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-3
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                        for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                            for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];

                                f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                f3m1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                kk0_A = KK2_A[i1][i2*W2+i3*W3+i4];
                                kk1p1_A = KK2_A[i1+1][i2*W2+i3*W3+i4];
                                kk1m1_A = KK2_A[i1-1][i2*W2+i3*W3+i4];
                                kk2p1_A = KK2_A[i1][(i2+1)*W2+i3*W3+i4];
                                kk2m1_A = KK2_A[i1][(i2-1)*W2+i3*W3+i4];
                                kk3p1_A = KK2_A[i1][i2*W2+(i3+1)*W3+i4];
                                kk3m1_A = KK2_A[i1][i2*W2+(i3-1)*W3+i4];
                                kk4p1_A = KK2_A[i1][i2*W2+i3*W3+(i4+1)];
                                kk4m1_A = KK2_A[i1][i2*W2+i3*W3+(i4-1)];

                                f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                f3m1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                kk0_B = KK2_B[i1][i2*W2+i3*W3+i4];
                                kk1p1_B = KK2_B[i1+1][i2*W2+i3*W3+i4];
                                kk1m1_B = KK2_B[i1-1][i2*W2+i3*W3+i4];
                                kk2p1_B = KK2_B[i1][(i2+1)*W2+i3*W3+i4];
                                kk2m1_B = KK2_B[i1][(i2-1)*W2+i3*W3+i4];
                                kk3p1_B = KK2_B[i1][i2*W2+(i3+1)*W3+i4];
                                kk3m1_B = KK2_B[i1][i2*W2+(i3-1)*W3+i4];
                                kk4p1_B = KK2_B[i1][i2*W2+i3*W3+(i4+1)];
                                kk4m1_B = KK2_B[i1][i2*W2+i3*W3+(i4-1)];

                                c_f0_a = 2.0 * (f0_A+0.5*kk0_A);
                                c_f0_b = 2.0 * (f0_B+0.5*kk0_B);
                                pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                KK3_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+0.5*kk1m1_A + f1p1_A+0.5*kk1p1_A - c_f0_a) + 
                                                        Ikh2mh1sq * (f2m1_A+0.5*kk2m1_A + f2p1_A+0.5*kk2p1_A - c_f0_a) +
                                                        Ikh2mh2sq * (f3m1_A+0.5*kk3m1_A + f3p1_A+0.5*kk3p1_A - c_f0_a) +
                                                        Ikh2mh3sq * (f4m1_A+0.5*kk4m1_A + f4p1_A+0.5*kk4p1_A - c_f0_a) - 
                                                        Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                FF_A[i1][i2*W2+i3*W3+i4] += KK3_A[i1][i2*W2+i3*W3+i4] / 3.0;

                                KK3_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+0.5*kk1m1_B + f1p1_B+0.5*kk1p1_B - c_f0_b) + 
                                                        Ikh2mh1sq * (f2m1_B+0.5*kk2m1_B + f2p1_B+0.5*kk2p1_B - c_f0_b) +
                                                        Ikh2mh2sq * (f3m1_B+0.5*kk3m1_B + f3p1_B+0.5*kk3p1_B - c_f0_b) +
                                                        Ikh2mh3sq * (f4m1_B+0.5*kk4m1_B + f4p1_B+0.5*kk4p1_B - c_f0_b) - 
                                                        Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                FF_B[i1][i2*W2+i3*W3+i4] += KK3_B[i1][i2*W2+i3*W3+i4] / 3.0;
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-33: CASE 3 KK3) = %lf sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }

                // RK4-4
                #pragma omp for private(xx1,xx2,xx3,xx4,\
                                        f0_A,f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A,\
                                        kk0_A,kk1p1_A,kk1m1_A,kk2p1_A,kk2m1_A,kk3p1_A,kk3m1_A,kk4p1_A,kk4m1_A,\
                                        f0_B,f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B,\
                                        kk0_B,kk1p1_B,kk1m1_B,kk2p1_B,kk2m1_B,kk3p1_B,kk3m1_B,kk4p1_B,kk4m1_B,\
                                        c_f0_a,c_f0_b,pot_a,pot_b,pot_c) schedule(runtime)
                for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                    for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                        for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                            for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                                xx1 = Box[0] + i1 * H[0];
                                xx2 = Box[2] + i2 * H[1];
                                xx3 = Box[4] + i3 * H[2];
                                xx4 = Box[6] + i4 * H[3];

                                f0_A = F_A[i1][i2*W2+i3*W3+i4];
                                f1p1_A = F_A[i1+1][i2*W2+i3*W3+i4];
                                f1m1_A = F_A[i1-1][i2*W2+i3*W3+i4];
                                f2p1_A = F_A[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_A = F_A[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_A = F_A[i1][i2*W2+(i3+1)*W3+i4];
                                f3m1_A = F_A[i1][i2*W2+(i3-1)*W3+i4];
                                f4p1_A = F_A[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_A = F_A[i1][i2*W2+i3*W3+(i4-1)];
                                kk0_A = KK3_A[i1][i2*W2+i3*W3+i4];
                                kk1p1_A = KK3_A[i1+1][i2*W2+i3*W3+i4];
                                kk1m1_A = KK3_A[i1-1][i2*W2+i3*W3+i4];
                                kk2p1_A = KK3_A[i1][(i2+1)*W2+i3*W3+i4];
                                kk2m1_A = KK3_A[i1][(i2-1)*W2+i3*W3+i4];
                                kk3p1_A = KK3_A[i1][i2*W2+(i3+1)*W3+i4];
                                kk3m1_A = KK3_A[i1][i2*W2+(i3-1)*W3+i4];
                                kk4p1_A = KK3_A[i1][i2*W2+i3*W3+(i4+1)];
                                kk4m1_A = KK3_A[i1][i2*W2+i3*W3+(i4-1)];

                                f0_B = F_B[i1][i2*W2+i3*W3+i4];
                                f1p1_B = F_B[i1+1][i2*W2+i3*W3+i4];
                                f1m1_B = F_B[i1-1][i2*W2+i3*W3+i4];
                                f2p1_B = F_B[i1][(i2+1)*W2+i3*W3+i4];
                                f2m1_B = F_B[i1][(i2-1)*W2+i3*W3+i4];
                                f3p1_B = F_B[i1][i2*W2+(i3+1)*W3+i4];
                                f3m1_B = F_B[i1][i2*W2+(i3-1)*W3+i4];
                                f4p1_B = F_B[i1][i2*W2+i3*W3+(i4+1)];
                                f4m1_B = F_B[i1][i2*W2+i3*W3+(i4-1)];
                                kk0_B = KK3_B[i1][i2*W2+i3*W3+i4];
                                kk1p1_B = KK3_B[i1+1][i2*W2+i3*W3+i4];
                                kk1m1_B = KK3_B[i1-1][i2*W2+i3*W3+i4];
                                kk2p1_B = KK3_B[i1][(i2+1)*W2+i3*W3+i4];
                                kk2m1_B = KK3_B[i1][(i2-1)*W2+i3*W3+i4];
                                kk3p1_B = KK3_B[i1][i2*W2+(i3+1)*W3+i4];
                                kk3m1_B = KK3_B[i1][i2*W2+(i3-1)*W3+i4];
                                kk4p1_B = KK3_B[i1][i2*W2+i3*W3+(i4+1)];
                                kk4m1_B = KK3_B[i1][i2*W2+i3*W3+(i4-1)];

                                c_f0_a = 2.0 * (f0_A+kk0_A);
                                c_f0_b = 2.0 * (f0_B+kk0_B);
                                pot_a = POTENTIAL_A(xx1,xx2,xx3,xx4);
                                pot_b = POTENTIAL_B(xx1,xx2,xx3,xx4);
                                pot_c = POTENTIAL_C(xx1,xx2,xx3,xx4);

                                KK4_A[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_A+kk1m1_A + f1p1_A+kk1p1_A - c_f0_a) + 
                                                            Ikh2mh1sq * (f2m1_A+kk2m1_A + f2p1_A+kk2p1_A - c_f0_a) + 
                                                            Ikh2mh2sq * (f3m1_A+kk3m1_A + f3p1_A+kk3p1_A - c_f0_a) +
                                                            Ikh2mh3sq * (f4m1_A+kk4m1_A + f4p1_A+kk4p1_A - c_f0_a) - 
                                                            Ik2h * 0.5 * (pot_a * c_f0_a + pot_c * c_f0_b); 

                                FF_A[i1][i2*W2+i3*W3+i4] += KK4_A[i1][i2*W2+i3*W3+i4] / 6.0;

                                KK4_B[i1][i2*W2+i3*W3+i4] = Ikh2mh0sq * (f1m1_B+kk1m1_B + f1p1_B+kk1p1_B - c_f0_b) + 
                                                            Ikh2mh1sq * (f2m1_B+kk2m1_B + f2p1_B+kk2p1_B - c_f0_b) + 
                                                            Ikh2mh2sq * (f3m1_B+kk3m1_B + f3p1_B+kk3p1_B - c_f0_b) + 
                                                            Ikh2mh3sq * (f4m1_B+kk4m1_B + f4p1_B+kk4p1_B - c_f0_b) - 
                                                            Ik2h * 0.5 * (pot_b * c_f0_b + pot_c * c_f0_a); 

                                FF_B[i1][i2*W2+i3*W3+i4] += KK4_B[i1][i2*W2+i3*W3+i4] / 6.0;
                            }
                        }
                    }
                }

                #pragma omp single nowait
                {
                    t_1_end = omp_get_wtime();
                    t_1_elapsed = t_1_end - t_1_begin;
                    t_full += t_1_elapsed;
                    if (!QUIET && TIMING) log->log("Elapsed time (omp-kk-34: CASE 3 KK4) = %.4e sec\n", t_1_elapsed);
                    t_1_begin = omp_get_wtime();
                }
            }
        }
        
        // .........................................................................................

        // NORMALIZATION AND TRUNCATION

        t_1_begin = omp_get_wtime();

        // Combine TAMask_A and TAMask_B

        if (!isFullGrid)  {
            #pragma omp parallel for schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4] || TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                TAMask_A[i1][i2*W2+i3*W3+i4] = 1;
                                TAMask_B[i1][i2*W2+i3*W3+i4] = 1;
                            }
                        }
                    }
                }
            }
        } 

        // Normalization

        /* ---------------- Surface A ---------------- */

        norm_A = 0.0;

        if (!isFullGrid)  {
            #pragma omp parallel for reduction (+:norm_A) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])
                                norm_A += std::abs(FF_A[i1][i2*W2+i3*W3+i4]*std::conj(FF_A[i1][i2*W2+i3*W3+i4]));
                        }
                    }
                }
            }
        }  
        else  {
            #pragma omp parallel for reduction (+:norm_A) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                        for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                            norm_A += std::abs(FF_A[i1][i2*W2+i3*W3+i4]*std::conj(FF_A[i1][i2*W2+i3*W3+i4]));
                        }
                    }
                }
            }
        }
        norm_A *= H[0] * H[1] * H[2] * H[3];

        if ( (tt + 1) % PERIOD == 0 )
            log->log("[Nonadiabatic4d] Normalization factor (A) = %.16e\n",norm_A);

        /* ---------------- Surface B ---------------- */

        norm_B = 0.0;

        if (!isFullGrid)  {
            #pragma omp parallel for reduction (+:norm_B) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])
                                norm_B += std::abs(FF_B[i1][i2*W2+i3*W3+i4]*std::conj(FF_B[i1][i2*W2+i3*W3+i4]));
                        }
                    }
                }
            }
        }  
        else  {
            #pragma omp parallel for reduction (+:norm_B) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                        for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                            norm_B += std::abs(FF_B[i1][i2*W2+i3*W3+i4]*std::conj(FF_B[i1][i2*W2+i3*W3+i4]));
                        }
                    }
                }
            }
        }
        norm_B *= H[0] * H[1] * H[2] * H[3];

        if ( (tt + 1) % PERIOD == 0 )
            log->log("[Nonadiabatic4d] Normalization factor (B) = %.16e\n",norm_B);

        norm = 1.0 / sqrt(norm_A + norm_B);

        /* ------------------------------ */

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-1 Norm) = %.4e sec\n", t_1_elapsed);
        t_1_begin = omp_get_wtime();

        /* ---------------- Surface A ---------------- */

        if (!isFullGrid)  {
            #pragma omp parallel for private(val) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])  {
                                val = norm * FF_A[i1][i2*W2+i3*W3+i4];
                                FF_A[i1][i2*W2+i3*W3+i4] = val;
                                F_A[i1][i2*W2+i3*W3+i4] = val;
                                PF_A[i1][i2*W2+i3*W3+i4] = std::abs(F_A[i1][i2*W2+i3*W3+i4] * std::conj(F_A[i1][i2*W2+i3*W3+i4]));
                            }
                        }
                    }
                }
            }
        }  
        else  {
            #pragma omp parallel for private(val) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                        for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                            val = norm * FF_A[i1][i2*W2+i3*W3+i4];
                            FF_A[i1][i2*W2+i3*W3+i4] = val;
                            F_A[i1][i2*W2+i3*W3+i4] = val;
                            PF_A[i1][i2*W2+i3*W3+i4] = std::abs(F_A[i1][i2*W2+i3*W3+i4] * std::conj(F_A[i1][i2*W2+i3*W3+i4]));
                        }
                    }
                }
            }
        }

        /* ---------------- Surface B ---------------- */

        if (!isFullGrid)  {
            #pragma omp parallel for private(val) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                val = norm * FF_B[i1][i2*W2+i3*W3+i4];
                                FF_B[i1][i2*W2+i3*W3+i4] = val;
                                F_B[i1][i2*W2+i3*W3+i4] = val;
                                PF_B[i1][i2*W2+i3*W3+i4] = std::abs(F_B[i1][i2*W2+i3*W3+i4] * std::conj(F_B[i1][i2*W2+i3*W3+i4]));
                            }
                        }
                    }
                }
            }
        }  
        else  {
            #pragma omp parallel for private(val) schedule(runtime)
            for (int i1 = EDGE; i1 < BoxShape[0]-EDGE; i1 ++)  {
                for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                    for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                        for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                            val = norm * FF_B[i1][i2*W2+i3*W3+i4];
                            FF_B[i1][i2*W2+i3*W3+i4] = val;
                            F_B[i1][i2*W2+i3*W3+i4] = val;
                            PF_B[i1][i2*W2+i3*W3+i4] = std::abs(F_B[i1][i2*W2+i3*W3+i4] * std::conj(F_B[i1][i2*W2+i3*W3+i4]));
                        }
                    }
                }
            }
        }

        /* -------------------------------------------- */

        t_1_end = omp_get_wtime();
        t_1_elapsed = t_1_end - t_1_begin;
        t_full += t_1_elapsed;
        t_truncate += t_1_elapsed;
        if (!QUIET && TIMING) log->log("Elapsed time (omp-e-1-2 FF) = %.4e sec\n", t_1_elapsed); 

        if ( (tt + 1) % PERIOD == 0 )
        {
            // REPORT MEASUREMENTS
            // ----------------------------------------------------------------------------
            if (isTrans)  {
            
                t_1_begin = omp_get_wtime();

                // position for x >= trans_x0
                idx_x0 = (int)std::round(std::round((trans_x0-Box[0])/H[0]));
                log->log("[Nonadiabatic4d] idx_x0 = %d\n", idx_x0 );
                log->log("[Nonadiabatic4d] trans_x0 = %lf\n", idx_x0 * H[0] + Box[0]);

                pftrans = 0.0;

                if (!isFullGrid) {
                    #pragma omp parallel for reduction (+:pftrans)
                    for (int i1 = idx_x0; i1 < x1_max_C; i1 ++)  {
                        for (int i2 = x2_min_C; i2 < x2_max_C; i2 ++)  {
                            for (int i3 = x3_min_C; i3 < x3_max_C; i3 ++)  {
                                for (int i4 = x4_min_C; i4 < x4_max_C; i4 ++)  {
                                    if (TAMask_A[i1][i2*W2+i3*W3+i4])
                                        pftrans += PF_A[i1][i2*W2+i3*W3+i4];
                                }
                            }
                        }
                    }
                }
                else  {
                    #pragma omp parallel for reduction (+:pftrans)
                    for (int i1 = idx_x0; i1 < BoxShape[0]-EDGE; i1 ++)  {
                        for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                            for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                                for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                                    pftrans += PF_A[i1][i2*W2+i3*W3+i4];
                                }
                            }
                        }
                    }
                }
                pftrans *= H[0] * H[1] * H[2] * H[3];
                log->log("[Nonadiabatic4d] Time %lf, Trans_A = %.16e\n", ( tt + 1 ) * kk, pftrans);

                pftrans = 0.0;

                if (!isFullGrid) {
                    #pragma omp parallel for reduction (+:pftrans)
                    for (int i1 = idx_x0; i1 < x1_max_C; i1 ++)  {
                        for (int i2 = x2_min_C; i2 < x2_max_C; i2 ++)  {
                            for (int i3 = x3_min_C; i3 < x3_max_C; i3 ++)  {
                                for (int i4 = x4_min_C; i4 < x4_max_C; i4 ++)  {
                                    if (TAMask_B[i1][i2*W2+i3*W3+i4])
                                        pftrans += PF_B[i1][i2*W2+i3*W3+i4];
                                }
                            }
                        }
                    }
                }
                else  {
                    #pragma omp parallel for reduction (+:pftrans)
                    for (int i1 = idx_x0; i1 < BoxShape[0]-EDGE; i1 ++)  {
                        for (int i2 = EDGE; i2 < BoxShape[1]-EDGE; i2 ++)  {
                            for (int i3 = EDGE; i3 < BoxShape[2]-EDGE; i3 ++)  {
                                for (int i4 = EDGE; i4 < BoxShape[3]-EDGE; i4 ++)  {
                                    pftrans += PF_B[i1][i2*W2+i3*W3+i4];
                                }
                            }
                        }
                    }
                }
                pftrans *= H[0] * H[1] * H[2] * H[3];
                log->log("[Nonadiabatic4d] Time %lf, Trans_B = %.16e\n", ( tt + 1 ) * kk, pftrans);

                t_1_end = omp_get_wtime();
                t_1_elapsed = t_1_end - t_1_begin; 
                t_full += t_1_elapsed;
                t_truncate += t_1_elapsed;
                if (!QUIET && TIMING) log->log("Elapsed time (omp-x-2 trans) = %.4e sec\n", t_1_elapsed); 
            }
            // ----------------------------------------------------------------------------
        }

        // Truncation and TA

        if ( !isFullGrid )
        {
            t_1_begin = omp_get_wtime();

            // Population-weighted TG criteria

            TolH_A = TolH * max(norm_A/(norm_A+norm_B), PTolMin);
            TolL_A = TolL * max(norm_A/(norm_A+norm_B), PTolMin);
            TolH_B = TolH * max(norm_B/(norm_A+norm_B), PTolMin);
            TolL_B = TolL * max(norm_B/(norm_A+norm_B), PTolMin);

            TolHd_A = TolHd;
            TolLd_A = TolLd;
            TolHd_B = TolHd;
            TolLd_B = TolLd;

            TolHd_A_sq = TolHd_A * TolHd_A;
            TolLd_A_sq = TolLd_A * TolLd_A;
            TolHd_B_sq = TolHd_B * TolHd_B;
            TolLd_B_sq = TolLd_B * TolLd_B;

            /* ---------------- Surface A ---------------- */

            #pragma omp parallel for private(b1,nx1,nx2,nx3,nx4,\
                                            f1p1_A,f1m1_A,f2p1_A,f2m1_A,f3p1_A,f3m1_A,f4p1_A,f4m1_A) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])  {
                                if (PF_A[i1][i2*W2+i3*W3+i4] < TolH_A)  {

                                    nx1 = int(TAMask_A[i1+1][i2*W2+i3*W3+i4]) + int(TAMask_A[i1-1][i2*W2+i3*W3+i4]);
                                    nx2 = int(TAMask_A[i1][(i2+1)*W2+i3*W3+i4]) + int(TAMask_A[i1][(i2-1)*W2+i3*W3+i4]);
                                    nx3 = int(TAMask_A[i1][i2*W2+(i3+1)*W3+i4]) + int(TAMask_A[i1][i2*W2+(i3-1)*W3+i4]);
                                    nx4 = int(TAMask_A[i1][i2*W2+i3*W3+(i4+1)]) + int(TAMask_A[i1][i2*W2+i3*W3+(i4-1)]);

                                    f1p1_A = (TAMask_A[i1+1][i2*W2+i3*W3+i4]) ? FF_A[i1+1][i2*W2+i3*W3+i4] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f1m1_A = (TAMask_A[i1-1][i2*W2+i3*W3+i4]) ? FF_A[i1-1][i2*W2+i3*W3+i4] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f2p1_A = (TAMask_A[i1][(i2+1)*W2+i3*W3+i4]) ? FF_A[i1][(i2+1)*W2+i3*W3+i4] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f2m1_A = (TAMask_A[i1][(i2-1)*W2+i3*W3+i4]) ? FF_A[i1][(i2-1)*W2+i3*W3+i4] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f3p1_A = (TAMask_A[i1][i2*W2+(i3+1)*W3+i4]) ? FF_A[i1][i2*W2+(i3+1)*W3+i4] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f3m1_A = (TAMask_A[i1][i2*W2+(i3-1)*W3+i4]) ? FF_A[i1][i2*W2+(i3-1)*W3+i4] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f4p1_A = (TAMask_A[i1][i2*W2+i3*W3+(i4+1)]) ? FF_A[i1][i2*W2+i3*W3+(i4+1)] : FF_A[i1][i2*W2+i3*W3+i4];
                                    f4m1_A = (TAMask_A[i1][i2*W2+i3*W3+(i4-1)]) ? FF_A[i1][i2*W2+i3*W3+(i4-1)] : FF_A[i1][i2*W2+i3*W3+i4];

                                    b1 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_A - f1m1_A)/(nx1*H[0]),2)) + \
                                         ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_A - f2m1_A)/(nx2*H[1]),2)) + \
                                         ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_A - f3m1_A)/(nx3*H[2]),2)) + \
                                         ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_A - f4m1_A)/(nx4*H[3]),2)) < TolHd_A_sq;

                                    if (b1)
                                        PF_A[i1][i2*W2+i3*W3+i4] = 0.0;
                                }
                            }
                        }
                    }                    
                }
            } 
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-1 TA) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();


            #pragma omp parallel for schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (PF_A[i1][i2*W2+i3*W3+i4] == 0.0)  {
                                TAMask_A[i1][i2*W2+i3*W3+i4] = 0;
                                F_A[i1][i2*W2+i3*W3+i4] = xZERO;
                            }  else  {
                                if (!TAMask_A[i1][i2*W2+i3*W3+i4])
                                    TAMask_A[i1][i2*W2+i3*W3+i4] = 1;
			                }
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-2 TA) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            /* ---------------- Surface B ---------------- */

            #pragma omp parallel for private(b1,nx1,nx2,nx3,nx4,\
                                            f1p1_B,f1m1_B,f2p1_B,f2m1_B,f3p1_B,f3m1_B,f4p1_B,f4m1_B) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                if (PF_B[i1][i2*W2+i3*W3+i4] < TolH_B)  {

                                    nx1 = int(TAMask_B[i1+1][i2*W2+i3*W3+i4]) + int(TAMask_B[i1-1][i2*W2+i3*W3+i4]);
                                    nx2 = int(TAMask_B[i1][(i2+1)*W2+i3*W3+i4]) + int(TAMask_B[i1][(i2-1)*W2+i3*W3+i4]);
                                    nx3 = int(TAMask_B[i1][i2*W2+(i3+1)*W3+i4]) + int(TAMask_B[i1][i2*W2+(i3-1)*W3+i4]);
                                    nx4 = int(TAMask_B[i1][i2*W2+i3*W3+(i4+1)]) + int(TAMask_B[i1][i2*W2+i3*W3+(i4-1)]);

                                    f1p1_B = (TAMask_B[i1+1][i2*W2+i3*W3+i4]) ? FF_B[i1+1][i2*W2+i3*W3+i4] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f1m1_B = (TAMask_B[i1-1][i2*W2+i3*W3+i4]) ? FF_B[i1-1][i2*W2+i3*W3+i4] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f2p1_B = (TAMask_B[i1][(i2+1)*W2+i3*W3+i4]) ? FF_B[i1][(i2+1)*W2+i3*W3+i4] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f2m1_B = (TAMask_B[i1][(i2-1)*W2+i3*W3+i4]) ? FF_B[i1][(i2-1)*W2+i3*W3+i4] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f3p1_B = (TAMask_B[i1][i2*W2+(i3+1)*W3+i4]) ? FF_B[i1][i2*W2+(i3+1)*W3+i4] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f3m1_B = (TAMask_B[i1][i2*W2+(i3-1)*W3+i4]) ? FF_B[i1][i2*W2+(i3-1)*W3+i4] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f4p1_B = (TAMask_B[i1][i2*W2+i3*W3+(i4+1)]) ? FF_B[i1][i2*W2+i3*W3+(i4+1)] : FF_B[i1][i2*W2+i3*W3+i4];
                                    f4m1_B = (TAMask_B[i1][i2*W2+i3*W3+(i4-1)]) ? FF_B[i1][i2*W2+i3*W3+(i4-1)] : FF_B[i1][i2*W2+i3*W3+i4];

                                    b1 = ((nx1 == 0) ? 0.0 : pow(std::abs(f1p1_B - f1m1_B)/(nx1*H[0]),2)) + \
                                         ((nx2 == 0) ? 0.0 : pow(std::abs(f2p1_B - f2m1_B)/(nx2*H[1]),2)) + \
                                         ((nx3 == 0) ? 0.0 : pow(std::abs(f3p1_B - f3m1_B)/(nx3*H[2]),2)) + \
                                         ((nx4 == 0) ? 0.0 : pow(std::abs(f4p1_B - f4m1_B)/(nx4*H[3]),2)) < TolHd_B_sq;

                                    if (b1)
                                        PF_B[i1][i2*W2+i3*W3+i4] = 0.0;
                                }
                            }
                        }
                    }                    
                }
            } 
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-3 TA) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (PF_B[i1][i2*W2+i3*W3+i4] == 0.0)  {
                                TAMask_B[i1][i2*W2+i3*W3+i4] = 0;
                                F_B[i1][i2*W2+i3*W3+i4] = xZERO;
                            }
                            else  {
                                if (!TAMask_B[i1][i2*W2+i3*W3+i4])
                                    TAMask_B[i1][i2*W2+i3*W3+i4] = 1;
                            }
                        }
                    }
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-3-4 TA) = %.4e sec\n", t_1_elapsed);

            // Rebuild TA box

            t_1_begin = omp_get_wtime();

            x1_min_A = BIG_NUMBER;
            x1_max_A = -BIG_NUMBER;
            x2_min_A = BIG_NUMBER;
            x2_max_A = -BIG_NUMBER;
            x3_min_A = BIG_NUMBER;
            x3_max_A = -BIG_NUMBER;
            x4_min_A = BIG_NUMBER;
            x4_max_A = -BIG_NUMBER;
            x1_min_B = BIG_NUMBER;
            x1_max_B = -BIG_NUMBER;
            x2_min_B = BIG_NUMBER;
            x2_max_B = -BIG_NUMBER;
            x3_min_B = BIG_NUMBER;
            x3_max_B = -BIG_NUMBER;
            x4_min_B = BIG_NUMBER;
            x4_max_B = -BIG_NUMBER;

            ta_size_A = 0;
            ta_size_B = 0;

            #pragma omp parallel for reduction(min: x1_min_A,x2_min_A,x3_min_A,x4_min_A,x1_min_B,x2_min_B,x3_min_B,x4_min_B) \
                                     reduction(max: x1_max_A,x2_max_A,x3_max_A,x4_max_A,x1_max_B,x2_max_B,x3_max_B,x4_max_B) \
                                     reduction(+: ta_size_A,ta_size_B) schedule(runtime)
            for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                    for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                        for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])  {
                                if (i1 < x1_min_A)  x1_min_A = i1;
                                if (i1 > x1_max_A)  x1_max_A = i1;
                                if (i2 < x2_min_A)  x2_min_A = i2;
                                if (i2 > x2_max_A)  x2_max_A = i2;
                                if (i3 < x3_min_A)  x3_min_A = i3;
                                if (i3 > x3_max_A)  x3_max_A = i3;
                                if (i4 < x4_min_A)  x4_min_A = i4;
                                if (i4 > x4_max_A)  x4_max_A = i4;
                                ta_size_A += 1;
                            }
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                if (i1 < x1_min_B)  x1_min_B = i1;
                                if (i1 > x1_max_B)  x1_max_B = i1;
                                if (i2 < x2_min_B)  x2_min_B = i2;
                                if (i2 > x2_max_B)  x2_max_B = i2;
                                if (i3 < x3_min_B)  x3_min_B = i3;
                                if (i3 > x3_max_B)  x3_max_B = i3;
                                if (i4 < x4_min_B)  x4_min_B = i4;
                                if (i4 > x4_max_B)  x4_max_B = i4;
                                ta_size_B += 1;
                            }
                        }
                    }
                }
            }

            x1_max_T = max(x1_max_A,x1_max_B);
            x2_max_T = max(x2_max_A,x2_max_B);
            x3_max_T = max(x3_max_A,x3_max_B);
            x4_max_T = max(x4_max_A,x4_max_B);
            x1_min_T = min(x1_min_A,x1_min_B);
            x2_min_T = min(x2_min_A,x2_min_B);
            x3_min_T = min(x3_min_A,x3_min_B);
            x4_min_T = min(x4_min_A,x4_min_B);

            x1_max_C = x1_max_T;
            x2_max_C = x2_max_T;
            x3_max_C = x3_max_T;
            x4_max_C = x4_max_T;
            x1_min_C = x1_min_T;
            x2_min_C = x2_min_T;
            x3_min_C = x3_min_T;
            x4_min_C = x4_min_T;

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-4 TA rebuild) = %.4e sec\n", t_1_elapsed);

            // TB

            /* ---------------- Surface A ---------------- */
            
            t_1_begin = omp_get_wtime();

            if (ta_size_A == 0)
                tb_size_A = 0;
            else  {
                #pragma omp parallel for reduction(merge: tmpVec) schedule(runtime)
                for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                    for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                        for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                            for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                if (TAMask_A[i1][i2*W2+i3*W3+i4])  {
                                    if (!TAMask_A[i1-1][i2*W2+i3*W3+i4] || !TAMask_A[i1+1][i2*W2+i3*W3+i4] || \
                                        !TAMask_A[i1][(i2-1)*W2+i3*W3+i4] || !TAMask_A[i1][(i2+1)*W2+i3*W3+i4] || \
                                        !TAMask_A[i1][i2*W2+(i3-1)*W3+i4] || !TAMask_A[i1][i2*W2+(i3+1)*W3+i4] || \
                                        !TAMask_A[i1][i2*W2+i3*W3+(i4-1)] || !TAMask_A[i1][i2*W2+i3*W3+(i4+1)])
                                        tmpVec.push_back(GridToIdx(i1,i2,i3,i4));
                                }
                            }
                        }
                    }
                }
                tmpVec.swap(TB_A);
                tmpVec.clear();
                tb_size_A = TB_A.size();
            }

            /* ---------------- Surface B ---------------- */

            if (ta_size_B == 0)
                tb_size_B = 0;
            else  {
                #pragma omp parallel for reduction(merge: tmpVec) schedule(runtime)
                for (int i1 = x1_min_C; i1 <= x1_max_C; i1 ++)  {
                    for (int i2 = x2_min_C; i2 <= x2_max_C; i2 ++)  {
                        for (int i3 = x3_min_C; i3 <= x3_max_C; i3 ++)  {
                            for (int i4 = x4_min_C; i4 <= x4_max_C; i4 ++)  {
                                if (TAMask_B[i1][i2*W2+i3*W3+i4])  {
                                    if (!TAMask_B[i1-1][i2*W2+i3*W3+i4] || !TAMask_B[i1+1][i2*W2+i3*W3+i4] || \
                                        !TAMask_B[i1][(i2-1)*W2+i3*W3+i4] || !TAMask_B[i1][(i2+1)*W2+i3*W3+i4] || \
                                        !TAMask_B[i1][i2*W2+(i3-1)*W3+i4] || !TAMask_B[i1][i2*W2+(i3+1)*W3+i4] || \
                                        !TAMask_B[i1][i2*W2+i3*W3+(i4-1)] || !TAMask_B[i1][i2*W2+i3*W3+(i4+1)])
                                        tmpVec.push_back(GridToIdx(i1,i2,i3,i4));
                                }
                            }
                        }
                    }
                }
                tmpVec.swap(TB_B);
                tmpVec.clear();
                tb_size_B = TB_B.size();
            }

            /* -------------------------------------------- */

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-5 TB) = %.4e sec\n", t_1_elapsed);

            // TA expansion

            /* ---------------- Surface A ---------------- */

            t_1_begin = omp_get_wtime();
            
            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < TB_A.size(); i++)
            {
                g1 = (int)(TB_A[i] / M1);
                g2 = (int)((TB_A[i]%M1) / M2);
                g3 = (int)((TB_A[i]%M2) / M3);
                g4 = (int)(TB_A[i] % M3);

                if (!TAMask_A[g1+1][g2*W2+g3*W3+g4])  {
                    if (g1+1 < BoxShape[0]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                }
                if (!TAMask_A[g1-1][g2*W2+g3*W3+g4])  {
                    if (g1-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                }
                if (!TAMask_A[g1][(g2+1)*W2+g3*W3+g4])  {
                    if (g2+1 < BoxShape[1]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                }
                if (!TAMask_A[g1][(g2-1)*W2+g3*W3+g4])  {
                    if (g2-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                }
                if (!TAMask_A[g1][g2*W2+(g3+1)*W3+g4])  {
                    if (g3+1 < BoxShape[2]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                }
                if (!TAMask_A[g1][g2*W2+(g3-1)*W3+g4])  {
                    if (g3-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                }
                if (!TAMask_A[g1][g2*W2+g3*W3+(g4+1)])  {
                    if (g4+1 < BoxShape[3]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                }
                if (!TAMask_A[g1][g2*W2+g3*W3+(g4-1)])  {
                    if (g4-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-6 TAEX-A) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min_A,x2_min_A,x3_min_A,x4_min_A) \
                                     reduction(max: x1_max_A,x2_max_A,x3_max_A,x4_max_A) \
                                     private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < tmpVec.size(); i ++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)((tmpVec[i]%M1) / M2);
                g3 = (int)((tmpVec[i]%M2) / M3);
                g4 = (int)(tmpVec[i] % M3);
                TAMask_A[g1][g2*W2+g3*W3+g4] = 1;
                x1_min_A = (g1 < x1_min_A) ? g1 : x1_min_A;
                x1_max_A = (g1 > x1_max_A) ? g1 : x1_max_A;
                x2_min_A = (g2 < x2_min_A) ? g2 : x2_min_A;
                x2_max_A = (g2 > x2_max_A) ? g2 : x2_max_A;
                x3_min_A = (g3 < x3_min_A) ? g3 : x3_min_A;
                x3_max_A = (g3 > x3_max_A) ? g3 : x3_max_A;
                x4_min_A = (g4 < x4_min_A) ? g4 : x4_min_A;
                x4_max_A = (g4 > x4_max_A) ? g4 : x4_max_A;
            }
            tmpVec.clear();

            ta_size_A = 0;
            #pragma omp parallel for reduction(+: ta_size_A) schedule(runtime)
            for (int i1 = x1_min_A; i1 <= x1_max_A; i1 ++)  {
                for (int i2 = x2_min_A; i2 <= x2_max_A; i2 ++)  {
                    for (int i3 = x3_min_A; i3 <= x3_max_A; i3 ++)  {
                        for (int i4 = x4_min_A; i4 <= x4_max_A; i4 ++)  {
                            if (TAMask_A[i1][i2*W2+i3*W3+i4])
                                ta_size_A += 1;
                        }
                    }    
                }
            }

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-7 TARB-A) = %.4e sec\n", t_1_elapsed);

            /* ---------------- Surface B ---------------- */

            t_1_begin = omp_get_wtime();
            
            #pragma omp parallel for reduction(merge: tmpVec) private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < TB_B.size(); i++)
            {
                g1 = (int)(TB_B[i] / M1);
                g2 = (int)((TB_B[i]%M1) / M2);
                g3 = (int)((TB_B[i]%M2) / M3);
                g4 = (int)(TB_B[i] % M3);

                if (!TAMask_B[g1+1][g2*W2+g3*W3+g4])  {
                    if (g1+1 < BoxShape[0]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1+1,g2,g3,g4));
                }
                if (!TAMask_B[g1-1][g2*W2+g3*W3+g4])  {
                    if (g1-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1-1,g2,g3,g4));
                }
                if (!TAMask_B[g1][(g2+1)*W2+g3*W3+g4])  {
                    if (g2+1 < BoxShape[1]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1,g2+1,g3,g4));
                }
                if (!TAMask_B[g1][(g2-1)*W2+g3*W3+g4])  {
                    if (g2-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1,g2-1,g3,g4));
                }
                if (!TAMask_B[g1][g2*W2+(g3+1)*W3+g4])  {
                    if (g3+1 < BoxShape[2]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1,g2,g3+1,g4));
                }
                if (!TAMask_B[g1][g2*W2+(g3-1)*W3+g4])  {
                    if (g3-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1,g2,g3-1,g4));
                }
                if (!TAMask_B[g1][g2*W2+g3*W3+(g4+1)])  {
                    if (g4+1 < BoxShape[3]-EDGE-1)
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4+1));
                }
                if (!TAMask_B[g1][g2*W2+g3*W3+(g4-1)])  {
                    if (g4-1 > EDGE)
                        tmpVec.push_back(GridToIdx(g1,g2,g3,g4-1));
                }
            }
            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-6 TAEX-B) = %.4e sec\n", t_1_elapsed);
            t_1_begin = omp_get_wtime();

            #pragma omp parallel for reduction(min: x1_min_B,x2_min_B,x3_min_B,x4_min_B) \
                                     reduction(max: x1_max_B,x2_max_B,x3_max_B,x4_max_B) \
                                     private(g1,g2,g3,g4) schedule(runtime)
            for (int i = 0; i < tmpVec.size(); i ++)  {
                g1 = (int)(tmpVec[i] / M1);
                g2 = (int)((tmpVec[i]%M1) / M2);
                g3 = (int)((tmpVec[i]%M2) / M3);
                g4 = (int)(tmpVec[i] % M3);
                TAMask_B[g1][g2*W2+g3*W3+g4] = 1;
                x1_min_B = (g1 < x1_min_B) ? g1 : x1_min_B;
                x1_max_B = (g1 > x1_max_B) ? g1 : x1_max_B;
                x2_min_B = (g2 < x2_min_B) ? g2 : x2_min_B;
                x2_max_B = (g2 > x2_max_B) ? g2 : x2_max_B;
                x3_min_B = (g3 < x3_min_B) ? g3 : x3_min_B;
                x3_max_B = (g3 > x3_max_B) ? g3 : x3_max_B;
                x4_min_B = (g4 < x4_min_B) ? g4 : x4_min_B;
                x4_max_B = (g4 > x4_max_B) ? g4 : x4_max_B;
            }
            tmpVec.clear();

            ta_size_B = 0;
            #pragma omp parallel for reduction(+: ta_size_B) schedule(runtime)
            for (int i1 = x1_min_B; i1 <= x1_max_B; i1 ++)  {
                for (int i2 = x2_min_B; i2 <= x2_max_B; i2 ++)  {
                    for (int i3 = x3_min_B; i3 <= x3_max_B; i3 ++)  {
                        for (int i4 = x4_min_B; i4 <= x4_max_B; i4 ++)  {
                            if (TAMask_B[i1][i2*W2+i3*W3+i4])
                                ta_size_B += 1;   
                        }
                    } 
                }
            }
            /* ------------------------------------------ */

            // Box range union
            x1_max_C = (x1_max_A > x1_max_B) ? x1_max_A : x1_max_B;
            x1_min_C = (x1_min_A < x1_min_B) ? x1_min_A : x1_min_B;
            x2_max_C = (x2_max_A > x2_max_B) ? x2_max_A : x2_max_B;
            x2_min_C = (x2_min_A < x2_min_B) ? x2_min_A : x2_min_B;
            x3_max_C = (x3_max_A > x3_max_B) ? x3_max_A : x3_max_B;
            x3_min_C = (x3_min_A < x3_min_B) ? x3_min_A : x3_min_B;
            x4_max_C = (x4_max_A > x4_max_B) ? x4_max_A : x4_max_B;
            x4_min_C = (x4_min_A < x4_min_B) ? x4_min_A : x4_min_B;

            t_1_end = omp_get_wtime();
            t_1_elapsed = t_1_end - t_1_begin;
            t_overhead += t_1_elapsed;
            if (!QUIET && TIMING) log->log("Elapsed time (omp-e-7 TARB-B) = %.4e sec\n", t_1_elapsed);

            if (isToggle && ((ta_size_A * 1.0 / GRIDS_TOT) > toggle_threshold || (ta_size_B * 1.0 / GRIDS_TOT) > toggle_threshold ))  {
                isFullGrid = true;
                log->log("[Nonadiabatic4d] threshold = %lf, ta_A/ntotal = %lf, ta_B/ntotal = %lf\n", toggle_threshold, ta_size_A*1.0/GRIDS_TOT, ta_size_B*1.0/GRIDS_TOT);
                log->log("[Nonadiabatic4d] Switching to full-grid mode at step: %d\n", tt + 1);
            }
        }

        if ( (tt + 1) % PERIOD == 0 )
        {   
            t_0_end = omp_get_wtime();
            t_0_elapsed = t_0_end - t_0_begin;

            if ( !QUIET ) log->log("[Nonadiabatic4d] Step: %d, Elapsed time: %.4e sec\n", tt + 1, t_0_elapsed);

            if ( !isFullGrid && !QUIET )  {
                log->log("[Nonadiabatic4d] norm (A) = %.2e, norm (B) = %.2e\n", norm_A, norm_B);
                log->log("[Nonadiabatic4d] TolH (A) = %.2e, TolH (B) = %.2e\n", TolH_A, TolH_B);
                log->log("[Nonadiabatic4d] TolL (A) = %.2e, TolL (B) = %.2e\n", TolL_A, TolL_B);
                log->log("[Nonadiabatic4d] TolHd (A) = %.2e, TolHd (B) = %.2e\n", TolHd_A, TolHd_B);
                log->log("[Nonadiabatic4d] TolLd (A) = %.2e, TolLd (B) = %.2e\n", TolLd_A, TolLd_B);
                log->log("[Nonadiabatic4d] TA size (A) = %d, TB size (A) = %d\n", ta_size_A, tb_size_A);
                log->log("[Nonadiabatic4d] TA Range (A) [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min_A, x1_max_A, x2_min_A, x2_max_A, x3_min_A, x3_max_A, x4_min_A, x4_max_A);
                log->log("[Nonadiabatic4d] TA / total grids (A) = %lf\n", ( ta_size_A * 1.0 ) / GRIDS_TOT);
                log->log("[Nonadiabatic4d] ExCount (A) = %d , ExLimit =  %d\n", Excount_A, ExLimit);
                log->log("[Nonadiabatic4d] TA size (B) = %d, TB size (B) = %d\n", ta_size_B, tb_size_B);
                log->log("[Nonadiabatic4d] TA Range (B) [%d, %d][%d, %d][%d, %d][%d, %d]\n", x1_min_B, x1_max_B, x2_min_B, x2_max_B, x3_min_B, x3_max_B, x4_min_B, x4_max_B);
                log->log("[Nonadiabatic4d] TA / total grids (B) = %lf\n", ( ta_size_B * 1.0 ) / GRIDS_TOT);
                log->log("[Nonadiabatic4d] ExCount (B) = %d , ExLimit = %d\n", Excount_B, ExLimit);
                log->log("[Nonadiabatic4d] Core computation time = %lf\n", t_truncate);
                log->log("[Nonadiabatic4d] Overhead time = %lf\n\n", t_overhead);
            }
            else if ( isFullGrid && !QUIET )  {

                log->log("[Nonadiabatic4d] Core computation time = %lf\n", t_full);
            }
            if ( !QUIET ) log->log("\n........................................................\n\n");
        }         
    } // Time iteration 

    delete F_A, F_B;
    delete FF_A, FF_B;
    delete PF_A, PF_B;
    delete KK1_A, KK1_B;
    delete KK2_A, KK2_B;
    delete KK3_A, KK3_B;
    delete KK4_A, KK4_B;

    if ( !isFullGrid )
        delete TAMask_A, TAMask_B;

    log->log("[Nonadiabatic4d] Evolve done.\n");
}
/* =============================================================================== */

/* Potential */

inline complex<double> Nonadiabatic4d::WaveA_1(double x1, double x2, double x3, double x4)
{
    // beta_0 = 8
    // E : 6000 cm^-1 = 6000/219474.63068 = 0.0273380116025718 hartree
    // P[0] = sqrt(2*m*E) = sqrt(2*2000*E) = 10.4571528826104 
    return exp(-8.0*((x1-Wave0[0])*(x1-Wave0[0])+(x2-Wave0[1])*(x2-Wave0[1])+(x3-Wave0[2])*(x3-Wave0[2])+(x4-Wave0[3])*(x4-Wave0[3])) + I/hb * 10.4571528826104 * x1);
}
/* ------------------------------------------------------------------------------- */

inline complex<double> Nonadiabatic4d::WaveB_1(double x1, double x2, double x3, double x4)
{
    return xZERO;
}
/* ------------------------------------------------------------------------------- */

inline double Nonadiabatic4d::VA_1(double x1, double x2, double x3, double x4)
{
    // Va = -a/2*(1+tanhx) + 0.5*ky^2*y^2
    // a : 6000 cm^-1 = 6000/219474.63068 = 0.0273380116025718 hartree
    return  -0.0136690058012859 * (1.0 + tanh(x1)) + 0.1025 * x2 * x2 + 0.1025 * x3 * x3 + 0.1025 * x4 * x4;
}
/* ------------------------------------------------------------------------------- */

inline double Nonadiabatic4d::VB_1(double x1, double x2, double x3, double x4)
{
    // Vb = 0.5 * ky * y^2
    // ky : 0.205
    return -0.00683450290064295 * (1.0 + tanh(x1)) + 0.1025 * (x2 * x2 + x3 * x3 + x4 * x4);
}
/* ------------------------------------------------------------------------------- */

inline double Nonadiabatic4d::VC_1(double x1, double x2, double x3, double x4)
{
    // Vc = c * exp[-b[(x-xc)^2+y^2]]
    // c: 250 cm^-1 = 250/219474.63068 = 0.00113908381677382 hartree
    // b: 1
    // xc: -2
    return  0.00113908381677382 * exp(-(x1*x1+4*x1+4 + x2*x2 + x3*x3 + x4*x4));
    //return  0.0;
}
/* =============================================================================== */

VectorXi Nonadiabatic4d::IdxToGrid(unsigned int idx)
{
    int x1 = (int)(idx / M1);
    int x2 = (int)((idx % M1) / M2);
    int x3 = (int)((idx % M2) / M3);
    int x4 = (int)(idx % M3);

    VectorXi grid;
    grid.resize(DIMENSIONS);
    grid << x1, x2, x3, x4;

    return grid;
}
/* ------------------------------------------------------------------------------- */

inline unsigned int Nonadiabatic4d::GridToIdx(int x1, int x2, int x3, int x4)
{
    return (unsigned int)(x1 * W1 + x2 * W2 + x3 * W3 + x4);
}
/* ------------------------------------------------------------------------------- */

