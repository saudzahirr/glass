#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

extern "C" {
    void sgssv_(int* M, int* N, float* A, int* LDA, float* B, float* X,
                int* MAX_ITER, float* TOL, int* INFO);

    void dgssv_(int* M, int* N, double* A, int* LDA, double* B, double* X,
                int* MAX_ITER, double* TOL, int* INFO);

    void cgssv_(int* M, int* N, std::complex<float>* A, int* LDA,
                std::complex<float>* B, std::complex<float>* X,
                int* MAX_ITER, float* TOL, int* INFO);

    void zgssv_(int* M, int* N, std::complex<double>* A, int* LDA,
                std::complex<double>* B, std::complex<double>* X,
                int* MAX_ITER, double* TOL, int* INFO);
}

template <typename T>
void print_vector(const std::vector<T>& v) {
    for (const auto& x : v)
        std::cout << x << " ";
    std::cout << "\n";
}

void test_sgssv() {
    std::cout << "\n--- Testing SGSSV (float, 4x4) ---\n";

    int N = 4, LDA = 4, MAX_ITER = 1000, INFO;
    float TOL = 1e-4f;

    float A[] = {
        10.0f, -1.0f,  2.0f,  0.0f,
        -1.0f, 11.0f, -1.0f,  3.0f,
         2.0f, -1.0f, 10.0f, -1.0f,
         0.0f,  3.0f, -1.0f,  8.0f
    };

    float B[] = { 6.0f, 25.0f, -11.0f, 15.0f };
    float X[] = { 0.0f, 0.0f, 0.0f, 0.0f };

    sgssv_(&N, &N, A, &LDA, B, X, &MAX_ITER, &TOL, &INFO);

    if (INFO == 0) {
        std::cout << "Converged solution: ";
        print_vector(std::vector<float>(X, X + N));
    } else {
        std::cout << "SGSSV failed with INFO = " << INFO << "\n";
    }
}

void test_dgssv() {
    std::cout << "\n--- Testing DGSSV (double, 4x4) ---\n";

    int N = 4, LDA = 4, MAX_ITER = 1000, INFO;
    double TOL = 1e-10;

    double A[] = {
        10.0, -1.0,  2.0,  0.0,
        -1.0, 11.0, -1.0,  3.0,
         2.0, -1.0, 10.0, -1.0,
         0.0,  3.0, -1.0,  8.0
    };

    double B[] = { 6.0, 25.0, -11.0, 15.0 };
    double X[] = { 0.0, 0.0, 0.0, 0.0 };

    dgssv_(&N, &N, A, &LDA, B, X, &MAX_ITER, &TOL, &INFO);

    if (INFO == 0) {
        std::cout << "Converged solution: ";
        print_vector(std::vector<double>(X, X + N));
    } else {
        std::cout << "DGSSV failed with INFO = " << INFO << "\n";
    }
}

void test_cgssv() {
    std::cout << "\n--- Testing CGSSV (complex<float>, 4x4) ---\n";

    using cf = std::complex<float>;

    int N = 4, LDA = 4, MAX_ITER = 1000, INFO;
    float TOL = 1e-4f;

    cf A[] = {
        {10.0f, 1.0f}, {-1.0f, 0.0f}, { 2.0f, 0.0f}, { 0.0f, 0.0f},
        {-1.0f, 0.0f}, {11.0f, 1.0f}, {-1.0f, 0.0f}, { 3.0f, 0.0f},
        { 2.0f, 0.0f}, {-1.0f, 0.0f}, {10.0f, 1.0f}, {-1.0f, 0.0f},
        { 0.0f, 0.0f}, { 3.0f, 0.0f}, {-1.0f, 0.0f}, { 8.0f, 1.0f}
    };

    cf B[] = { {6.0f, 1.0f}, {25.0f, 2.0f}, {-11.0f, 1.0f}, {15.0f, -1.0f} };
    cf X[] = { {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f} };

    cgssv_(&N, &N, A, &LDA, B, X, &MAX_ITER, &TOL, &INFO);

    if (INFO == 0) {
        std::cout << "Converged solution:\n";
        for (int i = 0; i < N; ++i)
            std::cout << X[i] << "\n";
    } else {
        std::cout << "CGSSV failed with INFO = " << INFO << "\n";
    }
}

void test_zgssv() {
    std::cout << "\n--- Testing ZGSSV (complex<double>, 4x4) ---\n";

    using cd = std::complex<double>;

    int N = 4, LDA = 4, MAX_ITER = 1000, INFO;
    double TOL = 1e-12;

    cd A[] = {
        {10.0, 1.0}, {-1.0, 0.0}, { 2.0, 0.0}, { 0.0, 0.0},
        {-1.0, 0.0}, {11.0, 1.0}, {-1.0, 0.0}, { 3.0, 0.0},
        { 2.0, 0.0}, {-1.0, 0.0}, {10.0, 1.0}, {-1.0, 0.0},
        { 0.0, 0.0}, { 3.0, 0.0}, {-1.0, 0.0}, { 8.0, 1.0}
    };

    cd B[] = { {6.0, 1.0}, {25.0, 2.0}, {-11.0, 1.0}, {15.0, -1.0} };
    cd X[] = { {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0} };

    zgssv_(&N, &N, A, &LDA, B, X, &MAX_ITER, &TOL, &INFO);

    if (INFO == 0) {
        std::cout << "Converged solution:\n";
        for (int i = 0; i < N; ++i)
            std::cout << X[i] << "\n";
    } else {
        std::cout << "ZGSSV failed with INFO = " << INFO << "\n";
    }
}

int main() {
    test_sgssv();
    test_dgssv();
    test_cgssv();
    test_zgssv();
    return 0;
}
