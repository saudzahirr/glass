#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H

#include "types.h"
#include "mangling.h"

/* =========================
 * FORTRAN API DECLARATIONS
 * ========================= */

fortran API_cgssv(INTEGER* N, COMPLEX* A, COMPLEX* B, COMPLEX* X,
                  INTEGER* MAX_ITER, REAL* TOL, INTEGER* STATUS);

fortran API_dgssv(INTEGER* N, DOUBLE* A, DOUBLE* B, DOUBLE* X,
                  INTEGER* MAX_ITER, DOUBLE* TOL, INTEGER* STATUS);

fortran API_sgssv(INTEGER* N, REAL* A, REAL* B, REAL* X,
                  INTEGER* MAX_ITER, REAL* TOL, INTEGER* STATUS);

fortran API_zgssv(INTEGER* N, DOUBLE_COMPLEX* A, DOUBLE_COMPLEX* B, DOUBLE_COMPLEX* X,
                  INTEGER* MAX_ITER, DOUBLE* TOL, INTEGER* STATUS);


#ifdef __cplusplus

    /* ==============
     * C++ INTERFACE
     * ============== */

    SUBROUTINE GAUSS_SEIDEL(INTEGER* N, REAL* A, REAL* B, REAL* X,
                            INTEGER* MAX_ITER, REAL* TOL, INTEGER* STATUS) {
        API_sgssv(N, A, B, X, MAX_ITER, TOL, STATUS);
    }

    SUBROUTINE GAUSS_SEIDEL(INTEGER* N, DOUBLE* A, DOUBLE* B, DOUBLE* X,
                            INTEGER* MAX_ITER, DOUBLE* TOL, INTEGER* STATUS) {
        API_dgssv(N, A, B, X, MAX_ITER, TOL, STATUS);
    }

    SUBROUTINE GAUSS_SEIDEL(INTEGER* N, COMPLEX* A, COMPLEX* B, COMPLEX* X,
                            INTEGER* MAX_ITER, REAL* TOL, INTEGER* STATUS) {
        API_cgssv(N, A, B, X, MAX_ITER, TOL, STATUS);
    }

    SUBROUTINE GAUSS_SEIDEL(INTEGER* N, DOUBLE_COMPLEX* A, DOUBLE_COMPLEX* B, DOUBLE_COMPLEX* X,
                            INTEGER* MAX_ITER, DOUBLE* TOL, INTEGER* STATUS) {
        API_zgssv(N, A, B, X, MAX_ITER, TOL, STATUS);
    }

#else  // C-only fallback

    /* ===========
     * C INTERFACE
     * ============ */

    #define GAUSS_SEIDEL(N, A, B, X, MAX_ITER, TOL, STATUS)          \
        _Generic((A),                                                \
            REAL*:            API_sgssv,                             \
            DOUBLE*:          API_dgssv,                             \
            COMPLEX*:         API_cgssv,                             \
            DOUBLE_COMPLEX*:  API_zgssv                              \
        )(N, A, B, X, MAX_ITER, TOL, STATUS)

#endif

#endif
