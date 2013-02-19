PROGRAM test

    USE DGmod

    IMPLICIT NONE
    INTEGER :: N,nx,num_elem,M
    REAL(KIND=8) dx,dxelem,V,x
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: nodes,wghts,HOLDER
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: C,Cinv,D1,D2
    INTEGER :: i, npad,k,j


   	nx = 48
    dx = 1D0/48
    dxelem = 1D0/8
    N = 5

    allocate(nodes(0:N),wghts(0:N),C(0:N,0:N),D1(0:N,0:N),D2(0:N,0:N),HOLDER(0:N),Cinv(0:N,0:N))

    num_elem = 8

    M = nx/num_elem

    write(*,*) 'N=',N
    write(*,*) 'M=',M
    write(*,*) 'M should be N+1!!'

    CALL gllnewton(N,nodes)
    write(*,*) 'nodes =', nodes

    CALL weights(N,nodes,wghts)
    write(*,*) 'w = ', wghts

    CALL Cmat(N,nodes,wghts,dx,dxelem,C)
    CALL GEPP_INV(C,N+1,Cinv)

    write(*,*) 'C'
    DO i=0,N
	    write(*,*) C(i,:)
    END DO
    
    write(*,*) 'Cinv'
    DO i=0,N
	    write(*,*) Cinv(I,:)
    END DO

    write(*,*) 'Checking Inverse..'
    HOLDER = 0D0
    DO I=0,N
        DO J=0,N
            DO K=0,N
    		        HOLDER(K) = C(I,K)*Cinv(K,J)
            END DO
            D1(I,J) = SUM(HOLDER)
        END DO
    END DO

    D2(:,:) = 0D0
    DO I=0,N
        DO J=0,N
            IF(I==J) THEN
                D2(I,J) = 1D0-D1(I,J)
            ELSE
                D2(I,J) = 0D0-D1(I,J)
            END IF
        END DO
    END DO

    write(*,*) ''
    write(*,*) 'Errors in C*Cinv'
    DO I=0,N
    write(*,*) D2(I,:)
    END DO

    write(*,*) 'C*Cinv'
    DO I=0,N
    write(*,*) D1(I,:)
    END DO

    HOLDER = 0D0
    DO I=0,N
        DO J=0,N
            DO K=0,N
	            HOLDER(K) = Cinv(I,K)*C(K,J)
            END DO
            D1(I,J) = SUM(HOLDER)
        END DO
    END DO

    D2(:,:) = 0D0
    DO I=0,N
        DO J=0,N
            IF(I==J) THEN
                D2(I,J) = 1D0-D1(I,J)
            ELSE
                D2(I,J) = 0D0-D1(I,J)
            END IF
        END DO
    END DO

    write(*,*) ''
    write(*,*) 'Errors in Cinv*C'
    DO I=0,N
    write(*,*) D2(I,:)
    END DO

    write(*,*) 'Cinv*C'
    DO I=0,N
    write(*,*) D1(I,:)
    END DO
    
    write(*,*) DACOS(-1D0)

    deallocate(nodes,wghts,C,Cinv,D1,D2)
END PROGRAM test
