C ------------------------------------------------------------------
C     PROGRAM Group V1.1
C     SPACE GROUP FINDER FOR PERIODIC STRUCTURES, BASED ON THE IDEA
C     THAT A GIVEN BRAVAIS LATTICE CAN ONLY SUPPORT A FINITE NUMBER 
C     OF ROTATION OPERATIONS, WHICH ARE THEN VERIFIED ONE BY ONE.
C     FROM THESE OPERATIONS WHICH SHOULD COMMUTE WITH THE SINGLE 
C     ELECTRON HAMILTONIAN*, ONE CAN COUNT THE ENERGY DEGENERACY OF
C     A REGULARLY MESHED FIRST BZ, THUS IN TOTAL ENERGY SUMMATION**
C     ONE ONLY NEEDS TO EVALUATE FOR A SMALLER SET OF K-POINTS.
C     
C                      AUG. 4, 1998  
C                      ORIGINAL PROGRAM BY AMES LAB
C                      MINOR REVISIONS BY JU LI (MIT)
C ------------------------------------------------------------------      
      
      PROGRAM Group
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMA=10,MML=2,EPS=1D-4)
C     **********************************************************
C     MMA: MAXIMUM NUMBER OF ATOMS PER UNIT CELL
C     MML: MAXIMUM NUMBER OF SHELLS TO FULLY ENCLOSE THE CENTER
C     warning: IF PARAMETERS MMA AND EPS ARE CHANGED HERE, 
C     OTHER PLACES MUST CHANGE TOO.
C     **********************************************************
      DIMENSION H(3,3), G(3,3), HH(3,3), RO(3,3)
      DIMENSION TAU(MMA,3),S(MMA,3),SS(MMA,3),SSS(MMA,3)
      DIMENSION X(3),XX(3,(2*MML+1)**3,3)
      DIMENSION NTYPE(MMA),NT(3)
C     MAXIMUM ORDER OF SPACE GROUP = 48
      DIMENSION ROT(48,3,3),ROTATION(48,3,3),TRANSLATION(48,3)
      CHARACTER *50 BUF, NAME
      DATA LP_GROUP /27/
      LOGICAL UNITARY, EQUIVALENT

      print *, 'Reading instructions...'
      print *, ' '
      
C     Read file header: "Name:", NAME, "Lattice:"
      READ (*, '(A50,/,A50,/,A50)') BUF,NAME,BUF

C     Read in Bravais lattice vectors in certain "characteristic" 
C     length unit L. Notice that a vector is a row in the matrix, 
C     and so (x, y, z) = (s1, s2, s3) * A.
      READ *, ((H(I,J),J=1,3),I=1,3)
      
C     Form the metric matrix: HH = HH+
      DO I=1,3
      DO J=1,3
      HH(I,J) = 0.D0
      DO K=1,3
      HH(I,J) = HH(I,J) + H(I,K)*H(J,K)
      ENDDO
      ENDDO
      ENDDO 
      
C     Calculate the inverse matrix: G = H^(-1)
      CALL MATINV (H,G,VOLUME)
C     Now G(:,m) is the reciprocal to lattice vector A(m,:)
      
      READ (*, '(A50)') BUF
C     "Number of atoms in the unit cell:"
      READ *, NPA
      IF (NPA.GT.MMA) THEN 
      PRINT *, 'Number of atoms in unit cell greater than', MMA
      STOP
      ENDIF
      
      READ (*, '(A50)') BUF
C     "Type (     X            Y            Z   ):"
      READ *, (NTYPE(N),(TAU(N,I),I=1,3),N=1,NPA)
C     Read in atom type and real space coordinates (in L)
      
C     Transform real space coordinates into coefficients
C     of lattice vectors (reduced coordinates)
      DO N = 1,NPA
      DO I = 1,3
      S(N,I) = G(1,I)*TAU(N,1)+G(2,I)*TAU(N,2)+G(3,I)*TAU(N,3)
      ENDDO
      ENDDO
      
C     ---------------------------------------------------
C     Find all possible rotational parts and then verify 
C     whether they could form space group operations.
C     ---------------------------------------------------
      
C     Step 1: select lattice vector transformations which 
C     preserve length for any of the Bravais lattice vectors:
      
      DO NX = 1,3
C     Pick lattice vector NX
      DO J = 1,3
      IF (J.EQ.NX) THEN 
      X(J)= 1.D0
      ELSE
      X(J) = 0.D0
      ENDIF
      ENDDO
C     Calculate its norm
      DDV = PPRODUCT(X, X, HH)
C     Go over surrounding shells
      IC = 0
      DO 10 N1 = -MML,MML 
      DO 10 N2 = -MML,MML 
      DO 10 N3 = -MML,MML 
      X(1) = N1
      X(2) = N2
      X(3) = N3
      DD = PPRODUCT(X, X, HH)
C     If norm not equal, look for new ones;
      IF (ABS(DD-DDV).GT.EPS) GOTO 10
C     else, increase the index
      IC = IC+1
C     and save this record for vector NX
      XX(NX,IC,1) = N1
      XX(NX,IC,2) = N2
      XX(NX,IC,3) = N3
 10   CONTINUE
C     total number of records for vector NX
      NT(NX) = IC
      WRITE (*, '(i2," lattice vectors has the same",
     A   " length as Bravais vector ",i1)') NT(NX), NX
      ENDDO

C     Step 2: select unitary transformation matrices
C     among NT(1)*NT(2)*NT(3) candidates:
      NTRANS = 0
      DO 20 IX = 1,NT(1)
      DO 20 IY = 1,NT(2)
      DO 20 IZ = 1,NT(3)
      DO J = 1,3
      RO(1,J) = XX(1,IX,J)
      RO(2,J) = XX(2,IY,J)
      RO(3,J) = XX(3,IZ,J)
      ENDDO
      IF (.NOT.UNITARY(RO,HH)) GOTO 20
      NTRANS = NTRANS + 1
      DO I=1,3
      DO J=1,3
      ROT(NTRANS,I,J) = RO(I,J)
      ENDDO
      ENDDO
 20   CONTINUE
      
      WRITE (*, '(/," Among the", i5, " possible candidates,",
     A " we find ", i2, " orthogonal matrices.")')
     A NT(1)*NT(2)*NT(3), NTRANS

C     Number of group operations:
      NGROUP = 0
C     Number of non-symorphic operations:
      NONSYM = 0

C     Loop over rotation
      DO 30 NG = 1,NTRANS
      
      DO N = 1,NPA
C     Operate on the atoms in the cluster
      DO J = 1,3
      SS(N,J) = 0.D0
      DO K=1,3
      SS(N,J) = SS(N,J) + S(N,K)*ROT(NG,K,J)
      ENDDO
      ENDDO
      ENDDO
      
      DO M = 1,NPA
C     For all non-primitive translations S(M,:) - SS(1,:)
      DO N = 1,NPA
      DO K = 1,3
      SSS(N,K) = SS(N,K) + S(M,K) - SS(1,K)
      ENDDO
      ENDDO
C     If the two structures are equivalent:
      IF (EQUIVALENT(NPA,NTYPE,SSS,NTYPE,S)) GOTO 40
      ENDDO
C     No, none of the translations are good
      GOTO 30
C     We have found a valid group operation
 40   NGROUP = NGROUP+1
      DO K = 1,3
      TRANSLATION(NGROUP,K) = 
     A  S(M,K)-SS(1,K)-NINT(S(M,K)-SS(1,K))
      DO J = 1,3
      ROTATION(NGROUP,K,J) = ROT(NG,K,J)
      ENDDO
      ENDDO
C     If it is also an non-symmorphic operation:
      IF (ABS(TRANSLATION(NGROUP,1)).GT.EPS.OR.
     A    ABS(TRANSLATION(NGROUP,2)).GT.EPS.OR.
     A    ABS(TRANSLATION(NGROUP,3)).GT.EPS) NONSYM = NONSYM+1
 30   CONTINUE
      
      write (*, '(/, " Among them, we found ", i2, 
     A " symmorphic operations,", /, 18x, "and ", i2, 
     A " non-symmorphic operations.",/)') NGROUP-NONSYM,NONSYM

      print *,'The group operations will be saved in file \"group-op\",'
      print *, 'both in terms of lattice vector coefficients (left),'
      print *, 
     A 'and in Cartesian frame and characteristic length L (right).'
      print *, ' '
      
C     Write group operations in file "group-op", in both 
C     reduced and real coordinates:
      OPEN (UNIT=LP_GROUP, STATUS='UNKNOWN', FORM='FORMATTED',
     A      FILE="group-op")
      WRITE (LP_GROUP,
     A '("Number of group operations found = ",I2)') NGROUP
      WRITE (LP_GROUP,
     A '("Number of non-symmorphic group operations =",I2)') NONSYM
      DO L = 1,NGROUP 
      DO I = 1,3
      X(I) = 0.
      DO J = 1,3
      X(I) = X(I) + TRANSLATION(L,J)*H(J,I)
      RO(I,J) = 0.D0
      DO K1=1,3
      DO K2=1,3
      RO(I,J) = RO(I,J) + G(I,K1)*ROTATION(L,K1,K2)*H(K2,J)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      WRITE (LP_GROUP,'(/,"No.",i2)') L
      WRITE (LP_GROUP,'("     |",3f9.5," |     |",3f9.5," |")') 
     A  ((ROTATION(L,I,J),J=1,3),(RO(I,J),J=1,3),I=1,3)
      WRITE (LP_GROUP,'(3x,"+ [",3f9.5," ]",3x,"+ [",3f9.5," ]")')
     A  (TRANSLATION(L,J),J=1,3), (X(J),J=1,3)
      ENDDO
      CLOSE (LP_GROUP)

      print *, 
     A 'Generating IBZ k-points (in reciprocal vector coeff.)...'
      print *, ' '
C     Generate IBZ k-points from a regularly meshed BZ:
      READ (*, '(/,A50)') BUF
C     "Number of IBZ k-point files:"
      READ *, NFILE
      DO I = 1,NFILE
      READ (*, '(/,A50)') BUF
C     "Mesh density + shift:"
      READ *, NX,NY,NZ,SX,SY,SZ
      READ (*, '(A50)') BUF
C     "filename:"
      READ (*, '(A50)') BUF
      NRK = IBZ_K_POINTS(NX,NY,NZ,SX,SY,SZ,NGROUP,ROTATION,BUF)
      PRINT *, NRK,  ' k-points saved in ', BUF(1:INDEX(BUF,' ')-1)
      ENDDO
      
      STOP
      END
      
      
C     GENERATE IBZ K-POINT GRID FROM A REGULARLY MESHED
C     BZ: NX, NY, NZ ARE THE NUMBER OF K-POINTS IN THREE
C     DIRECTIONS. SX, SY, SZ SHIFT THE GRID FROM THE ORIGIN.
      INTEGER FUNCTION IBZ_K_POINTS
     A (NX,NY,NZ,SX,SY,SZ,NGROUP,ROTATION,FILENAME)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MAXK=40000,EPS=1D-4)
      DIMENSION ROTATION(48,3,3)
      DIMENSION XK(3,MAXK),RKTRAN(3),KDEGEN(MAXK)
      CHARACTER *50 FILENAME
      DATA LP_KPTS /28/
      
      IF (NX*NY*NZ.GT.MAXK) THEN 
      PRINT *, 'IBZ_K_POINTS(): increase NMAX'
      STOP
      ENDIF
      
C     A uniform grid of k points:
      NK = 0
      DO I = 1,NX
      DO J = 1,NY
      DO K = 1,NZ
      NK = NK+1
      XK(1,NK) = (I-1+SX)/NX
      XK(2,NK) = (J-1+SY)/NY
      XK(3,NK) = (K-1+SZ)/NZ
      KDEGEN(NK) = 1
      ENDDO
      ENDDO
      ENDDO
      
      NRK = 0
      DO 10 I = 1,NK
C     Only RK will be recorded:
      IF (KDEGEN(I).EQ.0) GOTO 10
      NRK = NRK+1
      DO 20 L = 1,NGROUP
C     Operate on RK with all symmetry operations
      DO J = 1,3
      RKTRAN(J) =  ROTATION(L,J,1)*XK(1,I) 
     A           + ROTATION(L,J,2)*XK(2,I) 
     A           + ROTATION(L,J,3)*XK(3,I) 
C     Translate to first BZ
      RKTRAN(J) = RKTRAN(J) - NINT(RKTRAN(J))
C     Make it [0,1) for comparison:
      IF (RKTRAN(J).LT.0) RKTRAN(J) = RKTRAN(J)+1.D0
      ENDDO
C     Check for upward k-points equivalence 
      DO 30 K = I+1,NK
      IF (KDEGEN(K).EQ.0) GOTO 30
C     Both the rotated k-vector and its inverse
C     (time-reversal symmetry if no spin-polarization)
      IF ( (ABS(RKTRAN(1)-XK(1,K)).LT.EPS.AND.
     A      ABS(RKTRAN(2)-XK(2,K)).LT.EPS.AND.            
     A      ABS(RKTRAN(3)-XK(3,K)).LT.EPS) .OR.
     A     (ABS(1.D0-RKTRAN(1)-XK(1,K)).LT.EPS.AND.
     A      ABS(1.D0-RKTRAN(2)-XK(2,K)).LT.EPS.AND.            
     A      ABS(1.D0-RKTRAN(3)-XK(3,K)).LT.EPS) ) THEN 
      KDEGEN(I) = KDEGEN(I)+1
      KDEGEN(K) = 0
      ENDIF
 30   CONTINUE
 20   CONTINUE
 10   CONTINUE

C     THE K-POINTS WILL BE STORED IN RECIPROCAL LATTICE 
C     VECTOR COEFFICIENTS [0,1)
      OPEN (UNIT=LP_KPTS, STATUS='UNKNOWN', FORM='FORMATTED',
     A      FILE=FILENAME(1:INDEX(FILENAME,' ')-1))
      WRITE (LP_KPTS, *) NRK
      DO N = 1,NK
      IF (KDEGEN(N).GT.0) THEN 
C     RKTRAN(1) = 2.D0*(XK(1,N)*G(1,1)+XK(2,N)*G(1,2)+XK(3,N)*G(1,3))
C     RKTRAN(2) = 2.D0*(XK(1,N)*G(2,1)+XK(2,N)*G(2,2)+XK(3,N)*G(2,3))
C     RKTRAN(3) = 2.D0*(XK(1,N)*G(3,1)+XK(2,N)*G(3,2)+XK(3,N)*G(3,3))
      WRITE (LP_KPTS,'(4f17.14)') 
     A   XK(1,N),XK(2,N),XK(3,N),dble(KDEGEN(N))/NK
C    A   RKTRAN(1),RKTRAN(2),RKTRAN(3),dble(KDEGEN(N))/NK
      ENDIF
      ENDDO
      CLOSE(LP_KPTS)
      
      IBZ_K_POINTS = NRK
      RETURN
      END

      
C     SUBROUTINE TO EVALUATE THE INNER PRODUCT OF
C     TWO VECTORS UNDER ARBITRARY METRIC.
      DOUBLE PRECISION FUNCTION PPRODUCT (X, Y, HH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(3), HH(3,3), Y(3)
      PPRODUCT = 0.D0
      DO I = 1,3
      DO J = 1,3
      PPRODUCT = PPRODUCT + X(I) * HH(I,J) * Y(J)
      ENDDO
      ENDDO
      RETURN
      END
      
      
C     CHECK IF THE TRANSFORMATION R, EXPRESSED FROM AND INTO
C     LATTICE VECTOR COEFFICIENTS, IS UNITARY:  RHH+R+ = HH+.
      LOGICAL FUNCTION UNITARY (R, HH)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (EPS=1D-4)
      DIMENSION R(3,3),HH(3,3),AR(3,3)
      DO I=1,3
      DO J=1,3
      AR(I,J) = 0.D0
      DO K1=1,3
      DO K2=1,3
      AR(I,J) = AR(I,J) + R(I,K1) * HH(K1,K2) * R(J,K2)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      DO I=1,3
      DO J=1,3
      IF (ABS(HH(I,J)-AR(I,J)).GT.EPS) THEN 
      UNITARY = .FALSE.
      RETURN
      ENDIF
      ENDDO
      ENDDO
      UNITARY = .TRUE.
      RETURN
      END

      
C     Compare two periodic group of atoms to see if they 
C     are structurally equivalent, indices notwithstanding.
      LOGICAL FUNCTION EQUIVALENT (NPA,NTYPE1,S1,NTYPE2,S2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (MMA=10,EPS=1D-4)
      DIMENSION S1(MMA,3),S2(MMA,3)
      DIMENSION NTYPE1(MMA),NTYPE2(MMA)
C     Assume that no atoms in each group are close within EPS
      DO 10 M = 1,NPA
      DO N = 1,NPA
      IF (NTYPE1(M).EQ.NTYPE2(N).AND.
     A ABS(S1(M,1)-S2(N,1)-NINT(S1(M,1)-S2(N,1))).LT.EPS.AND.
     A ABS(S1(M,2)-S2(N,2)-NINT(S1(M,2)-S2(N,2))).LT.EPS.AND.
     A ABS(S1(M,3)-S2(N,3)-NINT(S1(M,3)-S2(N,3))).LT.EPS) GOTO 10
      ENDDO
      EQUIVALENT = .FALSE.
      RETURN
 10   CONTINUE
      EQUIVALENT = .TRUE.
      RETURN
      END

      
      SUBROUTINE MATINV(A,B,C)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(3,3),B(3,3)
      D11=A(2,2)*A(3,3)-A(2,3)*A(3,2)
      D22=A(3,3)*A(1,1)-A(3,1)*A(1,3)
      D33=A(1,1)*A(2,2)-A(1,2)*A(2,1)
      D12=A(2,3)*A(3,1)-A(2,1)*A(3,3)
      D23=A(3,1)*A(1,2)-A(3,2)*A(1,1)
      D31=A(1,2)*A(2,3)-A(1,3)*A(2,2)
      D13=A(2,1)*A(3,2)-A(3,1)*A(2,2)
      D21=A(3,2)*A(1,3)-A(1,2)*A(3,3)
      D32=A(1,3)*A(2,1)-A(2,3)*A(1,1)
      C=A(1,1)*D11+A(1,2)*D12+A(1,3)*D13
      B(1,1)=D11/C
      B(2,2)=D22/C
      B(3,3)=D33/C
      B(1,2)=D21/C
      B(2,3)=D32/C
      B(3,1)=D13/C
      B(2,1)=D12/C
      B(3,2)=D23/C
      B(1,3)=D31/C
      RETURN
      END
      
C ------------------------------------------------------------------
C     Convention:  row vectors, and {R|t}r = rR + t
C     * All rotational and translational operations commute with 
C     the total Hamiltonian H(r1,r2,..R1,R2..) where {R1,R2..} are 
C     nuclear coordinates. However if O{R1,R2..} = {R1,R2..}, then 
C     OH=HO even if O operates on {r1,r2,..} only. Furthermore, one 
C     assume the ground state charge density rho(r) to be the identity
C     representation of the symmetry group, which is self-consistent
C     but in principle not necessary.
C     ** Although IBZ is enough for total energy summation, it can 
C     lead to fictitious results in LDOS calculation if the symmetry 
C     operations contains non-primitive translations, such as for the
C     two Si atoms in 2H SiC, though physically equivalent, may have
C     different IBZ summed charges. The reason is simply although 
C     O\psi(k)=\psi(Ok) has the same energy as \psi(k), the charge
C     density on atom i in \psi(k) can shift to an equivalent but 
C     different site after O\psi(k): so although \sum_{BZ}E(k) = 
C     N\sum_{IBZ}E(k) and \sum_{BZ}\rho(k,r) has the full symmetry of 
C     the crystal, \sum_{IBZ}\rho(k,r) may not. This problem has a 
C     simple solution: one shows that charges belonging to one 
C     class of atoms (two Si above) never goes other classes for 
C     any k, so one just needs to symmetrize within each class.
C ------------------------------------------------------------------      
