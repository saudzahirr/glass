#include <stdio.h>
#include <complex.h>
#include "gauss-seidel.h"

void test_float() {
    int N = 4;
    float A[16] = {
        10.0f, -1.0f, 2.0f, 0.0f,
        -1.0f, 11.0f, -1.0f, 3.0f,
        2.0f, -1.0f, 10.0f, -1.0f,
        0.0f, 3.0f, -1.0f, 8.0f
    };
    float B[4] = {6.0f, 25.0f, -11.0f, 15.0f};
    float X[4] = {0};
    int max_iter = 1000;
    float tol = 1e-4f;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    printf("\n--- SGSSV (float) ---\nStatus: %d\n", status);
    for (int i = 0; i < N; i++) printf("X[%d] = %f\n", i, X[i]);
}

void test_double() {
    int N = 4;
    double A[16] = {
        10.0, -1.0, 2.0, 0.0,
        -1.0, 11.0, -1.0, 3.0,
        2.0, -1.0, 10.0, -1.0,
        0.0, 3.0, -1.0, 8.0
    };
    double B[4] = {6.0, 25.0, -11.0, 15.0};
    double X[4] = {0};
    int max_iter = 1000;
    double tol = 1e-10;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    printf("\n--- DGSSV (double) ---\nStatus: %d\n", status);
    for (int i = 0; i < N; i++) printf("X[%d] = %lf\n", i, X[i]);
}

void test_cfloat() {
    int N = 4;
    float _Complex A[16] = {
        10.0f + 1.0f*I, -1.0f + 0.0f*I, 2.0f + 0.0f*I, 0.0f + 0.0f*I,
        -1.0f + 0.0f*I, 11.0f + 1.0f*I, -1.0f + 0.0f*I, 3.0f + 0.0f*I,
        2.0f + 0.0f*I, -1.0f + 0.0f*I, 10.0f + 1.0f*I, -1.0f + 0.0f*I,
        0.0f + 0.0f*I, 3.0f + 0.0f*I, -1.0f + 0.0f*I, 8.0f + 1.0f*I
    };
    float _Complex B[4] = {
        6.0f + 1.0f*I, 25.0f + 2.0f*I, -11.0f + 1.0f*I, 15.0f - 1.0f*I
    };
    float _Complex X[4] = {0};
    int max_iter = 1000;
    float tol = 1e-4f;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    printf("\n--- CGSSV (complex float) ---\nStatus: %d\n", status);
    for (int i = 0; i < N; i++) {
        printf("X[%d] = (%f, %f)\n", i, crealf(X[i]), cimagf(X[i]));
    }
}

void test_cdouble() {
    int N = 4;
    double _Complex A[16] = {
        10.0 + 1.0*I, -1.0 + 0.0*I, 2.0 + 0.0*I, 0.0 + 0.0*I,
        -1.0 + 0.0*I, 11.0 + 1.0*I, -1.0 + 0.0*I, 3.0 + 0.0*I,
        2.0 + 0.0*I, -1.0 + 0.0*I, 10.0 + 1.0*I, -1.0 + 0.0*I,
        0.0 + 0.0*I, 3.0 + 0.0*I, -1.0 + 0.0*I, 8.0 + 1.0*I
    };
    double _Complex B[4] = {
        6.0 + 1.0*I, 25.0 + 2.0*I, -11.0 + 1.0*I, 15.0 - 1.0*I
    };
    double _Complex X[4] = {0};
    int max_iter = 1000;
    double tol = 1e-12;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    printf("\n--- ZGSSV (complex double) ---\nStatus: %d\n", status);
    for (int i = 0; i < N; i++) {
        printf("X[%d] = (%lf, %lf)\n", i, creal(X[i]), cimag(X[i]));
    }
}

int main() {
    test_float();
    test_double();
    test_cfloat();
    test_cdouble();
    return 0;
}
