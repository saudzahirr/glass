      SUBROUTINE DGSSV(M, N, A, LDA, B, X, MAX_ITER, TOL, INFO)
C     ----------------------------------------------------------------
C     Gauss-Seidel Solver for A * X = B using flat, row-major arrays
C     Double Precision version (LAPACK-style)
C
C     M        - INTEGER: number of rows in matrix A
C     N        - INTEGER: number of columns in matrix A (should be M for square)
C     A(*)     - DOUBLE PRECISION: matrix stored in row-major order, size >= M*N
C     LDA      - INTEGER: leading dimension, must be >= N
C     B(M)     - DOUBLE PRECISION: right-hand side vector
C     X(M)     - DOUBLE PRECISION: initial guess on input, solution on output
C     MAX_ITER - INTEGER: maximum number of iterations
C     TOL      - DOUBLE PRECISION: convergence tolerance
C     INFO     - INTEGER: status output
C                   0 = success (converged)
C                  >0 = did not converge within MAX_ITER
C                 <0 = illegal or zero diagonal detected
C     ----------------------------------------------------------------

      INTEGER M, N, LDA, MAX_ITER, I, J, K, INDEX
      DOUBLE PRECISION A(*), B(M), X(M), X_NEW(M)
      DOUBLE PRECISION S1, S2, DIFF, MAX_DIFF, TOL
      INTEGER INFO

      IF (M .NE. N) THEN
         INFO = -1
         RETURN
      END IF

      INFO = 1  ! Default: did not converge

      DO K = 1, MAX_ITER
         DO I = 1, M
            S1 = 0.0D0
            S2 = 0.0D0

C           S1 = sum over j = 1 to i - 1 of A(i,j) * X_NEW(j)
            DO J = 1, I - 1
               INDEX = (I - 1) * LDA + J
               S1 = S1 + A(INDEX) * X_NEW(J)
            END DO

C           S2 = sum over j = i + 1 to N of A(i,j) * X(j)
            DO J = I + 1, N
               INDEX = (I - 1) * LDA + J
               S2 = S2 + A(INDEX) * X(J)
            END DO

C           Diagonal A(i,i)
            INDEX = (I - 1) * LDA + I
            IF (A(INDEX) .EQ. 0.0D0) THEN
               INFO = -I
               RETURN
            END IF

            X_NEW(I) = (B(I) - S1 - S2) / A(INDEX)
         END DO

C        Compute max difference between new and old values
         MAX_DIFF = 0.0D0
         DO I = 1, M
            DIFF = ABS(X_NEW(I) - X(I))
            IF (DIFF .GT. MAX_DIFF) MAX_DIFF = DIFF
            X(I) = X_NEW(I)
         END DO

         IF (MAX_DIFF .LT. TOL) THEN
            INFO = 0  ! Success
            RETURN
         END IF
      END DO

C     If we reached here, did not converge in MAX_ITER
      RETURN
      END
