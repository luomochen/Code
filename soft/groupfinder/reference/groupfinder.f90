!************************************************************
! SPACE GROUP FINDER FOR PERIODIC STRUCTURES, BASED ON THE IDEA
! THAT A GIVEN BRAVAIS LATTICE CAN ONLY SUPPORT A FINITE NUMBER 
! OF ROTATION OPERATIONS, WHICH ARE THEN VERIFIED ONE BY ONE.
!************************************************************

program groupfinder
    use matrix_operation
    ! Declare variables.
    implicit none
    integer I, J, K, N1, N2, N3, MMA, NPA, MML
    integer NT(3)
    real(kind=8) VOLUME, EPS, DD, DDV   
    real(kind=8) X(3), H(3,3), HH(3,3), G(3,3)
    real(kind=8),allocatable :: TAU(:,:), S(:,:), XX(:,:,:)
    PARAMETER (MMA=10,MML=2,EPS=1D-4)

    ! Read Bravais lattice vectors.
    read(*,*) ((H(I,J),J=1,3),I=1,3)
    ! Form the metric matrix: HH = HH+.
    do I=1,3
        do J=1,3
            HH(I,J) = 0.D0
            do K=1,3
                HH(I,J) = HH(I,J) + H(I,K)*H(J,K)
            end do
        end do
    end do
    ! Calculate the inverse matrix: G = H^(-1)
    call MATINV (H,G,VOLUME)
    ! Read in atom type and real space coordinates.
    write(*,*) 'Please input the atom number of your unit cell:'
    read(*,*) NPA
    if ( NPA>MMA ) then
        write(*,*) 'Number of atoms in unit cell greater than', MMA
        stop
    end if
    allocate(TAU(NPA,3))
    allocate(S(NPA,3))
    read(*,*) ((TAU(I,J),J=1,3),I=1,NPA)
    do I = 1,NPA
        do J = 1,3
            S(I,J) = G(1,J)*TAU(I,1)+G(2,J)*TAU(I,2)+G(3,J)*TAU(I,3)
        end do
    end do
    write(*,*) S
    !---------------------------------------------------
    ! Find all possible rotational parts and then verify 
    ! whether they could form space group operations.
    !---------------------------------------------------
    
    ! Step 1: select lattice vector transformations which 
    ! preserve length for any of the Bravais lattice vectors:
    do I = 1, 3
        do J = 1, 3
            if ( J==I ) then
                X(J)=1.D0
            else
                X(J)=0.D0
            end if  
        end do
        call INPRODUCT(X,X,HH,DDV)
        K = 0
        do N1 = -MML, MML
            do N2 = -MML, MML
                do N3 = -MML, MML
                    X(1)=N1
                    X(2)=N2
                    X(3)=N3
                call INPRODUCT(X,X,HH,DD)
                if ( abs(DD-DDV) < EPS ) then
                    K=K+1
                end if
                end do
            end do
        end do
        allocate(XX(3,K,3))
        K = 0
        do N1 = -MML, MML
            do N2 = -MML, MML
                do N3 = -MML, MML
                    X(1)=N1
                    X(2)=N2
                    X(3)=N3
                call INPRODUCT(X,X,HH,DD)
                if ( abs(DD-DDV) < EPS ) then
                    K=K+1
                    XX(I,K,1) = N1
                    XX(I,K,2) = N2
                    XX(I,K,3) = N3
                end if
                end do
            end do
        end do
        NT(I)=K
    end do
    ! Step 2: select unitary transformation matrices
    ! among NT(1)*NT(2)*NT(3) candidates:
end program groupfinder