! ##############################################################
! Module containing 4th and 5th order Legendre polynomials 
! and their derivatives used in ppmwrap.f90 as part of test_advection_slskam_2d.f90
! By : Devin Light 11.29.2012
! ##############################################################

MODULE DGmod
	IMPLICIT NONE
    INTEGER, PARAMETER :: DOUBLE = KIND(1D0)
    

    PRIVATE :: DOUBLE

    ! #######################################################################
    ! node4 is vector of zeros of Jacobi polynomial
    ! associated with the problem, and are also the location of the GLL nodes
    ! #######################################################################

    REAL(KIND = DOUBLE), DIMENSION(0:4) :: node4 = (/ & 
            -1D0, &
            -0.654653670707978D0, &
            0D0, &
            0.654653670707977D0, &
            1D0 /)


    ! ##########################################################
    ! w4 is the weights used in the GLL quadrature 
    ! with w4(i) being the weight of the ith term in the sum
    ! ##########################################################

    REAL(KIND = DOUBLE), DIMENSION(0:4) :: w4 = (/ &
            1D0/10D0, &
            0.5444444444444456D0, &
            32D0/45D0,&
            0.5444444444444456D0,&
            1D0/10D0 /)
                   
    CONTAINS

    ! ########################################################################
    ! N-choose-k Function
    ! ########################################################################
       REAL(KIND=DOUBLE) FUNCTION choose(alpha,k)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: k
            REAL(KIND=DOUBLE), INTENT(IN) :: alpha
            INTEGER :: i
            REAL(KIND=DOUBLE) :: HOLDER
            
            HOLDER = 1D0

            DO i = 1,k
                HOLDER = HOLDER*((alpha-DBLE(k-i))/(DBLE(i)))
            END DO
            choose = HOLDER
        END FUNCTION choose
    
    ! ########################################################################
    ! Legendre Polynomial function of degree N
    ! ########################################################################
        REAL(KIND=DOUBLE) FUNCTION legendre(x,N)
            IMPLICIT NONE
            REAL(KIND=DOUBLE), INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: N
            INTEGER :: k
    
            HOLDER = 0.D0
            DO k = 0,N
                HOLDER = HOLDER + choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**k
            END DO

            legendre = HOLDER*(2**N)

        END FUNCTION legendre

    ! ########################################################################
    ! Derivative of Legendre Polyomial of degree N
    ! ########################################################################
        REAL(KIND=DOUBLE) FUNCTION dlegendre(x,N)
            IMPLICIT NONE
            REAL(KIND=DOUBLE),INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: N
            INTEGER :: k

            HOLDER = 0.D0
            DO k = 1,N
                HOLDER = HOLDER + k*choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**(k-1)
            END DO

            dlegendre = HOLDER*(2**N)

        END FUNCTION dlegendre

    ! ########################################################################
    ! 2nd Derivative of Legendre Polyomial of degree N
    ! ########################################################################
        REAL(KIND=DOUBLE) FUNCTION ddlegendre(x,N)
            IMPLICIT NONE
            REAL(KIND=DOUBLE), INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: N
            INTEGER :: k

            HOLDER = 0.D0
            DO k = 2,N
                HOLDER = HOLDER + k*(k-1)*choose(DBLE(N),k)*choose((N+k-1)/2D0,N)*x**(k-2)
            END DO

            ddlegendre = HOLDER*(2**N)

        END FUNCTION ddlegendre

    ! ########################################################################
    ! Subroutine for computing GLL nodes based on the derivative of N'th Order Legendre Polynomial
    ! ########################################################################
        SUBROUTINE gllnewton(N,nodes)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(OUT) :: nodes
            REAL(KIND=DOUBLE) :: xnew,xold,error,tol, PI
            INTEGER :: k

            PI = DACOS(-1D0)

            tol  = 10.D0**(-8)

            nodes(0) = -1D0
            nodes(N) = 1D0

            DO k = 1,N-1
                error = 1D0
                xold = -1D0*DCOS( ((2*k-1)/(2D0*(N-1)))*PI)

                DO WHILE (error>tol)
                    xnew = xold - (dlegendre(xold,N))/(1D0*ddlegendre(xold,N))
                    error = DABS(xnew-xold)
                    xold = xnew
                END DO
                nodes(k) = xold
            END DO
        END SUBROUTINE gllnewton

    ! ########################################################################
    ! Computing weights associated with N+1 nodes for quadratures on [-1,1]
    ! ########################################################################
        SUBROUTINE weights(N,nodes,wghts)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(OUT) :: wghts
            INTEGER :: k

            DO k = 0,N
                wghts(k) = 2D0/(N*(N+1)*(legendre(nodes(k),N))**2)
            END DO

        END SUBROUTINE weights

    ! ########################################################################
    ! Subroutine for filling in the C-matrix. Used for interchanging solution 
    ! between DG and PPM, computed using GLL quadrature
    ! ########################################################################
        SUBROUTINE Cmat(N,nodes,wghts,dx,dxelem,output)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes, wghts
            REAL(KIND=DOUBLE), DIMENSION(0:N,0:N), INTENT(OUT) :: output
            REAL(KIND=DOUBLE), INTENT(IN) :: dx,dxelem
            INTEGER :: m,j,k
            REAL(KIND=DOUBLE) :: HOLDER,dz,x

            dz = (2D0*dx)/(dxelem)
 
            HOLDER = 0.D0
            DO m=0,N
                DO j=0,N
                    HOLDER = 0D0
                    DO k=0,N
                        x = dz*(m+(nodes(k)+1)/2D0)-1D0
                        HOLDER = HOLDER + wghts(k)*phi(x,j,N,nodes)
                    END DO
                    output(m,j) = HOLDER/2D0
                END DO
            END DO

        END SUBROUTINE Cmat

    ! ########################################################################
    ! p4(x) computes the 4th order Legendre polynomial value 
    ! at x. dp4(x) compute the derivatives of these polys at x
    ! ########################################################################
        REAL FUNCTION p4(x)
            IMPLICIT NONE
            REAL(KIND = DOUBLE), INTENT(IN) :: x
            
            p4 = (35*x**4 - 30*x**2 + 3)/8
        END FUNCTION p4

        REAL FUNCTION dp4(x)
            IMPLICIT NONE
            REAL(KIND = DOUBLE), INTENT(IN) :: x
        
            dp4 = (35*x**3 - 15*x)/2

        END FUNCTION dp4

    ! ###########################################################
    ! phi computes the k'th basis function, a Lagrange interpolating
    ! polynomial, for a given set of nodes
    ! ###########################################################
        REAL(KIND=DOUBLE) FUNCTION phi(x,k,N,nodes)
            IMPLICIT NONE
            REAL(KIND=DOUBLE), INTENT(IN) :: x
            REAL(KIND=DOUBLE) :: HOLDER
            INTEGER, INTENT(IN) :: k, N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            INTEGER :: i = 0
            
            HOLDER = 1D0
            DO i = 0, N
                IF(i /= k) HOLDER = HOLDER * ( (x-nodes(i))/(nodes(k)-nodes(i)) )
            END DO
            
            phi = HOLDER

       END FUNCTION phi


    ! #############################################################################
    ! Subroutine D(N,nodes,output) computes the matrix of values of diff(phi(k,x),x) 
    ! evaluated at x = nodes(n); Used in calc of the Galerkin step
    ! #############################################################################
        SUBROUTINE Dmat(N,nodes,output)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            REAL(KIND=DOUBLE), DIMENSION(0:N, 0:N), INTENT(OUT) :: output
            INTEGER :: i = 0, j = 0
        
            DO i = 0, N
                DO j = 0, N
                    IF (i /= j) THEN
                      output(i,j) = (legendre(nodes(j),N))/(legendre(nodes(i),N)*(nodes(j)-nodes(i)))
                    ELSE
                      output(i,j) = 0.D0
                    END IF
                END DO
            END DO

            output(0,0) = -(N)*(N+1)/4D0
            output(N,N) = (N)*(N+1)/4D0

        END SUBROUTINE Dmat

        REAL(KIND = DOUBLE) FUNCTION phitld(x,N,nodes,j,Ain)
            REAL(KIND = DOUBLE), INTENT(IN) :: x
            INTEGER, INTENT(IN) :: j, N
            REAL(KIND = DOUBLE), DIMENSION(0:,1:), INTENT(IN) :: Ain
            REAL(KIND=DOUBLE), DIMENSION(0:N), INTENT(IN) :: nodes
            INTEGER :: i

            phitld = 0D0

            DO i = 0, N
                 phitld = phitld + Ain(i,j)*phi(x,i,N,nodes)
            END DO

        END FUNCTION phitld

        SUBROUTINE DGuvel(u,x,y,nx,ny,ntest) ! Computes initial u velocity distributions at DG nodes for each level of y
            IMPLICIT NONE
    
            INTEGER, INTENT(IN) :: nx,ny,ntest
            REAL(KIND=DOUBLE), INTENT(OUT), DIMENSION(1:nx,1:ny) :: u
            REAL(KIND=DOUBLE), INTENT(IN), DIMENSION(1:nx) :: x
            REAL(KIND=DOUBLE), INTENT(IN), DIMENSION(1:ny) :: y

            ! Local Parameters
            INTEGER :: i,j
            REAL(KIND=DOUBLE), DIMENSION(1:nx,1:ny) :: psi ! Stream function
            REAL(KIND=DOUBLE) :: PI

            PI = DACOS(-1D0)

            SELECT CASE(ntest)
                CASE(1:2) 
                ! Uniform flow: u=v=1, no t dependence
                    DO j=1,ny
                        DO i=1,nx
                            psi(i,j) = -x(i) + y(j)
                        END DO
                    END DO
                CASE(5:6,9)
                ! LeVeque (1996) Deformation
                    DO j=1,ny
                        DO i=1,nx
                            psi(i,j) = (1D0/PI)*(DSIN(PI*x(i))**2)*(DSIN(PI*y(j))**2)
                        END DO
                    END DO
                CASE(99)
                ! Uniform flow: u=1,v=0, no t dependence
                    DO j=1,ny
                        DO i=1,nx
                            psi(i,j) = y(j)
                        END DO
                    END DO
            END SELECT

            ! Compute u velocities from stream function
            DO j = 2,ny
               u(:,j) = (psi(:,j) - psi(:,j-1))/(y(j) - y(j-1))
            END DO
            ! For the first row, use forward divided difference to stay within domain
            u(:,1) = (psi(:,2) - psi(:,1))/(y(2) - y(1))

        END SUBROUTINE DGuvel

        SUBROUTINE GEPP_INV (M,N,Minv)
        !
        ! Subroutine to perform the partial-pivoting Gaussian elimination.
        ! A(N,N) is the original matrix in the input and transformed matrix
        ! plus the pivoting element ratios below the diagonal in the output.
        ! INDX(N) records the pivoting order. 
        ! Ainv is found by performing the same operations on I(N,N) and then
        ! solving the related systems for each column of Ainv. [A|I] -> [A|Y]
        !
          IMPLICIT NONE
          INTEGER, INTENT (IN) :: N
          INTEGER :: I,J,K,ITMP
          INTEGER, DIMENSION (N) :: INDX
          REAL :: C1,PI,PI1,PJ
          REAL(KIND=DOUBLE), INTENT (IN), DIMENSION (N,N) :: M
          REAL(KIND=DOUBLE), INTENT(OUT), DIMENSION(N,N) :: Minv
          REAL(KIND=DOUBLE), DIMENSION(N,N) :: Y,Tmp1,Tmp2,A
          REAL, DIMENSION (N) :: C

        ! Initialize A and Minv
        A(:,:) = M(:,:)
        Minv(:,:) = 0.D0

        ! Initialize Y as I(N,N)
          Y(:,:) = 0D0
          DO I = 1,N
            Y(I,I) = 1D0
          END DO

        !
        ! Initialize the index
        !
          DO I = 1, N
            INDX(I) = I
          END DO
        !
        ! Select largest absval element, one from each row
        !
!          DO I = 1, N
!            C1= 0.0
!            DO J = 1, N
!              C1 = DMAX1(C1,ABS(A(I,J)))
!            END DO
!            C(I) = C1
!          END DO


          DO J = 1, N-1
            
            ! Select pivoting (largest) element from each column
            PI1 = 0.0
            DO I = J, N
              PI = DABS(A(INDX(I),J)) !/C(INDX(I))
              IF (PI.GT.PI1) THEN
                PI1 = PI
                K   = I
              END IF
            END DO
        !
        ! Interchange the rows via INDX(N) to record pivoting order
        !
            ITMP    = INDX(J)
            INDX(J) = INDX(K)
            INDX(K) = ITMP
            DO I = J+1, N
              PJ  = A(INDX(I),J)/A(INDX(J),J)
        !
        ! Record pivoting ratios below the diagonal
        !
              A(INDX(I),J) = 0D0!PJ
        !
        ! Modify other elements accordingly
        !
              DO K = J+1, N
                A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
              END DO
              DO K = 1,N
                Y(INDX(I),K) = Y(INDX(I),K)-PJ*Y(INDX(J),K)
              END DO
            END DO
          END DO

        ! Swap rows to get it back to the correct form
        Tmp1 = A
        Tmp2 = Y
        DO I=1,N
            A(I,:) = Tmp1(INDX(I),:)
            Y(I,:) = Tmp2(INDX(I),:)
        END DO

        ! To find Minv (for the PIVOTED matrix), solve n-systems using back substitution
        ! Anew*Ainv(:,k) = Y(:,k) k=1..n
          DO K = 1,N
            Minv(N,K) = Y(N,K)/A(N,N)
            DO I = N-1,1,-1
                Minv(I,K) = Y(I,K)
                DO J = I+1,N
                    Minv(I,K) = Minv(I,K) - A(I,J)*Minv(J,K)
                END DO
            Minv(I,K) = Minv(I,K)/A(I,I)
            END DO
          END DO

        END SUBROUTINE GEPP_INV

END MODULE DGmod
