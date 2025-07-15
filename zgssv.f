      SUBROUTINE ZGSSV(M, N, A, LDA, B, X, MAX_ITER, TOL, INFO)
C     ----------------------------------------------------------------
C     Gauss-Seidel solver for A * X = B (complex double, row-major)
C     ----------------------------------------------------------------
      INTEGER M, N, LDA, MAX_ITER, I, J, K, INDEX
      DOUBLE COMPLEX A(*), B(M), X(M), X_NEW(M)
      DOUBLE COMPLEX S1, S2
      DOUBLE PRECISION TOL, DIFF, MAX_DIFF
      INTEGER INFO

      IF (M .NE. N) THEN
         INFO = -1
         RETURN
      END IF

      INFO = 1

      DO K = 1, MAX_ITER
         DO I = 1, M
            S1 = (0.0D0, 0.0D0)
            S2 = (0.0D0, 0.0D0)

            DO J = 1, I - 1
               INDEX = (I - 1) * LDA + J
               S1 = S1 + A(INDEX) * X_NEW(J)
            END DO

            DO J = I + 1, N
               INDEX = (I - 1) * LDA + J
               S2 = S2 + A(INDEX) * X(J)
            END DO

            INDEX = (I - 1) * LDA + I
            IF (A(INDEX) .EQ. (0.0D0, 0.0D0)) THEN
               INFO = -I
               RETURN
            END IF

            X_NEW(I) = (B(I) - S1 - S2) / A(INDEX)
         END DO

         MAX_DIFF = 0.0D0
         DO I = 1, M
            DIFF = ABS(X_NEW(I) - X(I))
            IF (DIFF .GT. MAX_DIFF) MAX_DIFF = DIFF
            X(I) = X_NEW(I)
         END DO

         IF (MAX_DIFF .LT. TOL) THEN
            INFO = 0
            RETURN
         END IF
      END DO

      RETURN
      END
