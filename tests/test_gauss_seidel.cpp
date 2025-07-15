#include <iostream>
#include <complex>
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
    float X[4] = {};
    int max_iter = 1000;
    float tol = 1e-4f;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    std::cout << "\n--- SGSSV (float) ---\nStatus: " << status << "\n";
    for (int i = 0; i < N; i++) std::cout << "X[" << i << "] = " << X[i] << "\n";
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
    double X[4] = {};
    int max_iter = 1000;
    double tol = 1e-10;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    std::cout << "\n--- DGSSV (double) ---\nStatus: " << status << "\n";
    for (int i = 0; i < N; i++) std::cout << "X[" << i << "] = " << X[i] << "\n";
}

void test_cfloat() {
    using cf = std::complex<float>;
    int N = 4;
    cf A[16] = {
        {10.0f, 1.0f}, {-1.0f, 0.0f}, {2.0f, 0.0f}, {0.0f, 0.0f},
        {-1.0f, 0.0f}, {11.0f, 1.0f}, {-1.0f, 0.0f}, {3.0f, 0.0f},
        {2.0f, 0.0f}, {-1.0f, 0.0f}, {10.0f, 1.0f}, {-1.0f, 0.0f},
        {0.0f, 0.0f}, {3.0f, 0.0f}, {-1.0f, 0.0f}, {8.0f, 1.0f}
    };
    cf B[4] = { {6.0f, 1.0f}, {25.0f, 2.0f}, {-11.0f, 1.0f}, {15.0f, -1.0f} };
    cf X[4] = {};
    int max_iter = 1000;
    float tol = 1e-4f;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    std::cout << "\n--- CGSSV (complex float) ---\nStatus: " << status << "\n";
    for (int i = 0; i < N; i++) std::cout << "X[" << i << "] = " << X[i] << "\n";
}

void test_cdouble() {
    using cd = std::complex<double>;
    int N = 4;
    cd A[16] = {
        {10.0, 1.0}, {-1.0, 0.0}, {2.0, 0.0}, {0.0, 0.0},
        {-1.0, 0.0}, {11.0, 1.0}, {-1.0, 0.0}, {3.0, 0.0},
        {2.0, 0.0}, {-1.0, 0.0}, {10.0, 1.0}, {-1.0, 0.0},
        {0.0, 0.0}, {3.0, 0.0}, {-1.0, 0.0}, {8.0, 1.0}
    };
    cd B[4] = { {6.0, 1.0}, {25.0, 2.0}, {-11.0, 1.0}, {15.0, -1.0} };
    cd X[4] = {};
    int max_iter = 1000;
    double tol = 1e-12;
    int status;

    GAUSS_SEIDEL(&N, A, B, X, &max_iter, &tol, &status);

    std::cout << "\n--- ZGSSV (complex double) ---\nStatus: " << status << "\n";
    for (int i = 0; i < N; i++) std::cout << "X[" << i << "] = " << X[i] << "\n";
}

int main() {
    test_float();
    test_double();
    test_cfloat();
    test_cdouble();
    return 0;
}
