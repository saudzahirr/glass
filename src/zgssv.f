C     ====================================================================
C     This file is part of glass - General Linear Algebra Subroutines
C
C     Copyright (C) 2025  Saud Zahir
C
C     glass is free software: you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published by
C     the Free Software Foundation, either version 3 of the License, or
C     (at your option) any later version.
C
C     glass is distributed in the hope that it will be useful,
C     but WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with glass.  If not, see <https://www.gnu.org/licenses/>.
C
C     ====================================================================
C       ZGSSV   -   Gauss-Seidel Solver for A * X = B
C     ====================================================================
C       Description:
C       ------------------------------------------------------------------
C         Iterative Gauss-Seidel solver for solving linear systems of
C         equations A * X = B, where A is a square N×N matrix in
C         row-major flat array format. Complex double precision version.
C
C         On input:  X contains initial guess
C         On output: X contains solution
C
C         Convergence is based on maximum absolute difference per iteration.
C     ====================================================================
C       Arguments:
C       ------------------------------------------------------------------
C         N         : INTEGER              -> size of the matrix (NxN)
C         A(*)      : DOUBLE COMPLEX       -> flat array, row-major matrix A
C         B(N)      : DOUBLE COMPLEX       -> right-hand side vector
C         X(N)      : DOUBLE COMPLEX       -> input: initial guess, output: solution
C         MAX_ITER  : INTEGER              -> max number of iterations
C         TOL       : DOUBLE PRECISION     -> convergence tolerance (L∞ norm)
C         INFO      : INTEGER              -> return code:
C                                                 0 = success
C                                                >0 = did not converge
C                                                <0 = illegal or zero diagonal
C     ====================================================================
      SUBROUTINE ZGSSV(N, A, B, X, MAX_ITER, TOL, INFO)

C   I m p l i c i t   T y p e s
C   ------------------------------------------------------------------
      IMPLICIT NONE

C   D u m m y   A r g u m e n t s
C   ------------------------------------------------------------------
      INTEGER              :: N, MAX_ITER, INFO
      DOUBLE COMPLEX       :: A(*), B(N), X(N)
      DOUBLE PRECISION     :: TOL

C   L o c a l   V a r i a b l e s
C   ------------------------------------------------------------------
      INTEGER              :: I, J, K, INDEX
      DOUBLE COMPLEX       :: X_NEW(N)
      DOUBLE COMPLEX       :: S1, S2
      DOUBLE PRECISION     :: DIFF, MAX_DIFF

C   I n i t i a l   S t a t u s
C   ------------------------------------------------------------------
      INFO = 1   ! Default: did not converge

C   M a i n   I t e r a t i o n   L o o p
C   ------------------------------------------------------------------
      DO K = 1, MAX_ITER
         DO I = 1, N
            S1 = (0.0D0, 0.0D0)
            S2 = (0.0D0, 0.0D0)

C           Compute sum: S1 = sum_{j=1}^{i-1} A(i,j) * X_NEW(j)
            DO J = 1, I - 1
               INDEX = (I - 1) * N + J
               S1 = S1 + A(INDEX) * X_NEW(J)
            END DO

C           Compute sum: S2 = sum_{j=i+1}^{N} A(i,j) * X(j)
            DO J = I + 1, N
               INDEX = (I - 1) * N + J
               S2 = S2 + A(INDEX) * X(J)
            END DO

C           Check diagonal element A(i,i)
            INDEX = (I - 1) * N + I
            IF (A(INDEX) .EQ. (0.0D0, 0.0D0)) THEN
               INFO = -I
               RETURN
            END IF

C           Update X_NEW(i)
            X_NEW(I) = (B(I) - S1 - S2) / A(INDEX)
         END DO

C        C o n v e r g e n c e   C h e c k
C        ----------------------------------------------------------------
         MAX_DIFF = 0.0D0
         DO I = 1, N
            DIFF = ABS(X_NEW(I) - X(I))
            IF (DIFF .GT. MAX_DIFF) MAX_DIFF = DIFF
            X(I) = X_NEW(I)
         END DO

         IF (MAX_DIFF .LT. TOL) THEN
            INFO = 0  ! Success
            RETURN
         END IF
      END DO

C   N o n - C o n v e r g e n c e   E x i t
C   ------------------------------------------------------------------
      RETURN
      END
