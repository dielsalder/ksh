PROGRAM KSUPERHELIXT
   IMPLICIT NONE
   
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: oldX, oldY, oldZ
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: newX, newY, newZ
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: newXYZ_AMAT, newXYZ, MINIMAT
   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: newXYZ_BMAT, XCROSS, Xsoln,SING

   CHARACTER(LEN=100)                      :: INFILE,ARG_MAXITR,ARG_NATOMS
   CHARACTER(LEN=100)                      :: ARG_NFRAMES
   INTEGER                                 :: I, COUNT
   INTEGER                                 :: ierr
   INTEGER                                 :: NFRAMES,MAXITR,FRSH
   INTEGER                                 :: IFRAME
   INTEGER, PARAMETER                      :: MAXRECS = 1024
   DOUBLE PRECISION                        :: PI, RHO, RTD
   
   DOUBLE PRECISION ::  COEF1, COEF2, COEF3, COEF4, COEF5, COEF6, COEF7, COEF8, COEF9
   DOUBLE PRECISION ::  O1, O2, O3, R
   DOUBLE PRECISION :: dPHI,dTHETA, PHI_L,THETA_M, PHID, THED, VPHI, VTHETA
   INTEGER :: IATOM, IPHI, ITHETA, J

   INTEGER(KIND=4)                 :: NR
   DOUBLE PRECISION                :: START, FINISH
   INTEGER,PARAMETER               :: MX = 3
   CHARACTER(LEN=8), DIMENSION(1)  :: NAMES = (/' '/)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                                            ON  PRBLOCK  U  W
   INTEGER,SAVE,DIMENSION(4)       :: KPVEC = (/ 1, 000000, -1,69 /)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   DOUBLE PRECISION,DIMENSION(MX)                    :: D(MX)
   DOUBLE PRECISION,DIMENSION(MX)                    :: WORK(2*MX)

   DOUBLE PRECISION :: X0, Y0, X02, Y02, R_SVD
   DOUBLE PRECISION :: RES
   DOUBLE PRECISION :: newZ_MAX, newZ_MIN
   DOUBLE PRECISION :: BEST_RES, BEST_RAD, BEST_PIT, BEST_PHI,  &
                                              BEST_THETAPIT, &
                                              BEST_THE, &
                                    BEST_S1, BEST_S2, BEST_S3, BEST_S0
   DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: RESMAT,      &
                                              RSVDMAT,     &
                                              PITMAT,      &
                                              THETAPITMAT, &
                                              PHIDMAT,     &
                                              THEDMAT,     &
                                              SING1,SING2,SING3,SING0 
                                              
   INTEGER,     DIMENSION(:)  ,ALLOCATABLE :: LOWESTRES
   INTEGER      :: MIN_IPHI, MIN_ITHETA

   DOUBLE PRECISION :: PHISTART,PHIEND, THESTART,THEEND
   DOUBLE PRECISION :: PHIS,THES,PHIN,THEN,PHIRANGE,THERANGE

   DOUBLE PRECISION :: VECIX, VECIY, VECIp1X, VECIp1Y, MAGI, MAGIp1
   DOUBLE PRECISION :: DOTIIp1, DOTOVERMAGS, THETASTEP, &
           THETAHELIX, ZETAPITCH, KAPPA, PITCH, THEDSTEP
   
   CHARACTER(LEN=100) :: ARG_PHISTART, ARG_PHIEND, ARG_THEEND, ARG_THESTART

   PI                = 3.141592653589793 ! Comment
   RHO               = 0.001D+00         ! FOR SVD
   RTD               = 180.0/PI          ! FOR PITCH


!                                                                           !
! 7777777777777777777777777777777777777777777777777777777777777777777777777 !
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL !
!
 
   CALL CPU_TIME(START)

   ! open out input trajectory that only contains the coordinates
   !  being analyzed. Cannot contain the first cpptraj comments line.
   CALL getarg(1,INFILE) 
   OPEN(51,FILE=INFILE,ACTION='READ',STATUS='OLD')
   CLOSE(1)        

   ! this is the practical resolution, ie the number of degree steps
   !  that the ranges, phi and theta, are split into. 
   ! IF the ranges are different for phi and theta, the resolution of 
   !  phi and theta will be different, ie different increments rotated.
   CALL getarg(2,ARG_MAXITR)
   READ(ARG_MAXITR,*) MAXITR

   ! this is just the number of atoms in each frame of the trajectory.
   ! this only gets used to index the frames, so if only one frame,
   ! used.
   CALL getarg(3,ARG_NATOMS)
   READ(ARG_NATOMS,*) NR  !Here, we no longer need NR to be set by

   ! this is used to set how many frames you are going to analyze,
   ! in sequence. cannot exceed tha maximum number of frames of the 
   ! trajetory. gets used to index frames.
   CALL getarg(4,ARG_NFRAMES)
   READ(ARG_NFRAMES,*) NFRAMES  

   ! set variables that control the bounds of the search over the unit sphere
   !  noting that instead of reading these in from prompt, eventually we can
   !  set these values via BFGS, for example...
   !  Vagelis HK...
   CALL getarg(5,ARG_PHISTART)
   READ(ARG_PHISTART,*) PHISTART
   
   CALL getarg(6,ARG_PHIEND)
   READ(ARG_PHIEND,*) PHIEND

   CALL getarg(7,ARG_THESTART)
   READ(ARG_THESTART,*) THESTART

   CALL getarg(8,ARG_THEEND)
   READ(ARG_THEEND,*) THEEND

   ! allocate the array that stores the coordinates
   ! 3 = x,y,z; NR = number of atoms per frame NFRAMES = number of frames
   ALLOCATE(x          (3*NR*NFRAMES),stat=ierr)
   IF(IERR /= 0) WRITE(*,*) "Allocate error in x, dude..."
   READ(51,*) x
   
   ! allocate array to store components, e_X, e_Y, e_Z
   ALLOCATE(oldX       (NR))
   ALLOCATE(oldY       (NR))
   ALLOCATE(oldZ       (NR))
   ! allocate array to store rotated components, e_x, e_y, e_z
   ALLOCATE(newX       (NR))
   ALLOCATE(newY       (NR))
   ALLOCATE(newZ       (NR))
   
   ! allocate a single array to store one rotated set of coordinates x,y,z
   !  that have been linearized per our linear problem AX=B
   ALLOCATE(newXYZ_AMAT(NR,3))
   ! allocate array for our linear problem, B matrix (area of circle)
   ALLOCATE(newXYZ_BMAT(NR))
   ! allocate array for rotation trajectory printing
   ALLOCATE(newXYZ     (NR,3))
   !allocate array for SVA and SVD routines
   ALLOCATE(SING       (3*MX))
   ALLOCATE(Xsoln      (3*MX))
   ALLOCATE(XCROSS     (3*NR))
   ! allocate arrays for finding the best fit parameters, the dimensions of
   !  these arrays precisely map to the grid over the unit sphere that is 
   !  searched, controlled by PhiStart, PhiEnd, ThetaStart, ThetaEnd and 
   !  just as importantly, MaxItr from above...
   ALLOCATE(RESMAT (MAXITR,MAXITR))
   ALLOCATE(RSVDMAT(MAXITR,MAXITR))
   ALLOCATE(PITMAT (MAXITR,MAXITR))
   ALLOCATE(THETAPITMAT(MAXITR,MAXITR))
   ALLOCATE(PHIDMAT(MAXITR,MAXITR))
   ALLOCATE(THEDMAT(MAXITR,MAXITR))
   ALLOCATE(SING1 (MAXITR,MAXITR))
   ALLOCATE(SING2 (MAXITR,MAXITR))
   ALLOCATE(SING3 (MAXITR,MAXITR))
   ALLOCATE(SING0 (MAXITR,MAXITR))
   ! allocated array that stores the values after MINLOC function finds
   !  the lowest residual in ResMat above...
   ALLOCATE(LOWESTRES(2))
   ! allocated array that finds local minima in the solution surface...
   ALLOCATE(MINIMAT(MAXITR,3))


!                                                                           !
!                                                                           !
! 6666666666666666666666666666666666666666666666666666666666666666666666666 !
! 9999999999999999999999999999999999999999999999999999999999999999999999999 !
!                                                                           !
!                                                                           !


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  MASTER CONTROL LOOP: ONE FRAME AT A TIME
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set loop counter up
    COUNT = 0

DO IFRAME = 1,NFRAMES   ! For each frame, for each atom, load XYZs

    ! For each structure/frame in the input coordinates, build arrays
    FRSH = ( IFRAME -1 ) * ( 3 ) * ( NR )
    DO IATOM = 1,NR
     oldX(IATOM) = x( 3*(IATOM-1) + 1 + FRSH )
     oldY(IATOM) = x( 3*(IATOM-1) + 2 + FRSH )
     oldZ(IATOM) = x( 3*(IATOM-1) + 3 + FRSH )
    END DO

!  Initialize the rotated components arrays
   DO IATOM = 1,NR
      newX(IATOM)=oldX(IATOM)
      newY(IATOM)=oldY(IATOM)
      newZ(IATOM)=oldZ(IATOM)
   END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           
!  START LOOP FOR ROTATIONS OVER THE UNIT SPHERE ALONG PHI AND THETA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !Convert input values from degrees to radians
   PHIS     = PHISTART * PI / 180.0
   THES     = THESTART * PI / 180.0

   PHIN     = PHIEND   * PI / 180.0
   THEN     = THEEND   * PI / 180.0
   
   !Find the range in Phi and Theta over which to search, and for loop control
   PHIRANGE = PHIN - PHIS
   THERANGE = THEN - THES

   !Find the step size through Phi and Theta, and for loop control
   !  Remind: Spherical coordinates: 0 < phi < PI ; 0 < theta < 2*PI
   dTHETA      = THERANGE / MAXITR
   dPHI        = PHIRANGE / MAXITR

 !Now launch the loops over the appropriate ranges in Phi and Theta, 
 ! incremented based on the resolution (step size) based on MAXITR
 DO IPHI    = 2,MAXITR-1   ! Loop-PHI
  DO ITHETA = 2,MAXITR-1   ! Loop-THETA

     !Get the Phi and Theta angles, in radians, back from the loop counters
     ! so we can apply the rotation matrices to the input coordinates
      PHI_L   = ( IPHI   * dPHI )   + PHIS
      THETA_M = ( ITHETA * dTHETA ) + THES

      !Convert PHI_L and THETA_M to degrees for writing
      PHID = PHI_L   * 180.0 /PI
      THED = THETA_M * 180.0 /PI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ROTATE COORDINATES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Setup factors of elements of the transformation matix
        O1  = SIN(PHI_L) * COS(THETA_M)
        O2  = SIN(PHI_L) * SIN(THETA_M)
        O3  = COS(PHI_L)
        R   = SQRT(1-(O1*O1))

        !Setup elements of the transformation matrix. 
        ! See EQN(4) in accompanying PDF documentation.
        !
        !     ( x )   ( COEF1  COEF2  COEF3 ) ( X )
        !     ( y ) = ( COEF4  COEF5  COEF6 ) ( Y )
        !     ( z )   ( COEF7  COEF8  COEF9 ) ( Z )
        !
        COEF1 =   SQRT(1.0-(O1*O1))
        COEF2 =   0
        COEF3 =   O1
        COEF4 =  -O1*O2 / R
        COEF5 =   O3    / R
        COEF6 =   O2
        COEF7 =  -O1*O3 / R
        COEF8 =  -O2    / R
        COEF9 =   O3

! Linear Algebra: loop over IATOMs, we transform the coordinates
   DO IATOM = 1,NR         ! Loop-Atoms
      newX(IATOM) = COEF1*oldX(IATOM) + COEF2*oldY(IATOM) + COEF3*oldZ(IATOM)
      newY(IATOM) = COEF4*oldX(IATOM) + COEF5*oldY(IATOM) + COEF6*oldZ(IATOM)
      newZ(IATOM) = COEF7*oldX(IATOM) + COEF8*oldY(IATOM) + COEF9*oldZ(IATOM)
   END DO
   DO IATOM = 1,NR         ! Loop-Atoms
      ! Build arrays for printing rotated coordinates trajectory...
      newXYZ(IATOM,1) = newX(IATOM)
      newXYZ(IATOM,2) = newY(IATOM)
      newXYZ(IATOM,3) = newZ(IATOM)
   END DO ! IATOM_complete transformation of coordinate set(phi,theta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SOLVE THE LINEAR PROBLEM, A * X = B
!
!   See EQN(9) in accompanying documentation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   !B matrix
   DO IATOM = 1,NR
      newXYZ_BMAT(IATOM) = newX(IATOM)*newX(IATOM) + newY(IATOM)*newY(IATOM)
   END DO
   
   !A matrix
   newXYZ_AMAT(:,1) = 1
   DO IATOM = 1,NR
      newXYZ_AMAT(IATOM,2) = newX(IATOM)
      newXYZ_AMAT(IATOM,3) = newY(IATOM)
   END DO

! FIND THE BEST ESTIMATE OF X USING SVD:
! Call the SVA routine from lawson.f90 to get our parameters
!  See accompanying documentation for the variables below.
   CALL SVA (newXYZ_AMAT,NR,NR,3,NR,newXYZ_BMAT,SING,KPVEC,NAMES,1,D,WORK)

! Calculate the radius and circle center from the output  parameters >>>
!  Back-out the actual parameters, R, X0, and Y0 from the fit results
    X0  = newXYZ_AMAT(2,3) /2  
    Y0  = newXYZ_AMAT(3,3) /2 
    X02 = X0*X0
    Y02 = Y0*Y0
    R_SVD = SQRT( X02 + Y02 + (newXYZ_AMAT(1,3)) )

! Calculate the pitch >>>
! We-ll take a trivial approach, where the angle of the helix
!  we need to find the pitch comes from the sum of the individual
!  angles between successive steps, ie. steps along the helix.
   THETAHELIX = 0.0D+00   
   DO IATOM = 1, NR-1
      VECIX    = newX(IATOM)   - X0                !ith x crd
      VECIY    = newY(IATOM)   - Y0                !ith y crd
      VECIp1X  = newX(IATOM+1) - X0                !i+1th x crd
      VECIp1Y  = newY(IATOM+1) - Y0                !i+1th y crd
       DOTIIp1 = VECIX*VECIp1X + VECIY*VECIp1Y     !dot i,i+1
       MAGI    = SQRT(VECIX*VECIX + &              !mag of ith vector
                      VECIY*VECIY)
       MAGIp1  = SQRT(VECIp1X*VECIp1X + &     !mag of i+1th vector
                      VECIp1Y*VECIp1Y)

       ! ALPHA = ARCCOS(A*B/|A||B|), A is ith vec, B is i+1th vec...
       DOTOVERMAGS = DOTIIp1/(MAGI*MAGIp1)        !value holder
       THETASTEP   = ACOS(DOTOVERMAGS)            !angle, in rad., bet steps
       THEDSTEP    = THETASTEP*RTD
      
       !testing algorithm...
       !WRITE (*,*) 'PHID: ', PHID, 'THETAD: ', THED
       !WRITE (*,*) 'IATOMS: ',IATOM, IATOM+1,'THETA STEP (DEG.): ',THEDSTEP
             
       THETAHELIX  = THETAHELIX + THEDSTEP       !sum over steps, total angle

       !! SPIT = MAGI*THETASTEP     !s is distance along arc bet endpoints
                                    ! of vectors IVECI and IVECIplus1, where
                                    ! s=r*theta.
   END DO
   ! now that we have the angle of the helix, calculate Kappa
   !
   ! HK: Try MAXVAL intrinsic function to find ZETAPITCH 
         newZ_MAX = MAXVAL( newZ(:), NR )![,MASK] )
         newZ_MIN = MINVAL( newZ(:), NR )![,MASK] )
   !   see: http://gcc.gnu.org/onlinedocs/gfortran/MAXVAL.html
   !
       ZETAPITCH = ABS (newZ_MAX - newZ_MIN)   !height, e_z, of helix
       KAPPA     = ZETAPITCH * RTD / THETAHELIX   !this needs to be in rad.
   ! calculate the pitch, finally, using Kappa
       PITCH     = ABS ( KAPPA * 2*PI )          

       !testing algorithm...
       !WRITE(*,*) 'PHID:      ', PHID,     'THED:     ', THED      
       !WRITE(*,*) 'newZ(NR):  ', newZ(NR), 'newZ(1):  ', newZ(1)
       !WRITE(*,*) 'ZETAPITCH: ', ZETAPITCH
       !WRITE(*,*) 'KAPPA:     ', KAPPA
       !WRITE(*,*) 'PITCH:     ', PITCH
       !WRITE(*,*) 'RESIDUAL:  ', RES
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   CALCULATE THE RESIDUAL: HELIX
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculate the residual as per Vageli >>>
   RES = 0.0D+00
   DO I = 1,NR
    RES = RES + ABS( &
                     & (newX(I)-X0)*(newX(I)-X0) + &
                     & (newY(I)-Y0)*(newY(I)-Y0) + &
                     & (R_SVD)*(R_SVD)           - &
                     & (2*R_SVD)*( SQRT (          &
                     & (newX(I)-X0)*(newX(I)-X0) + &
                     & (newY(I)-Y0)*(newY(I)-Y0) ) ) )
   END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   MONITOR THE RESIDUAL DURING THE SCAN...
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!   WRITE (*,FMT=57) 'RES:   ', RES,  &
!                    'PHID:  ', PHID, &
!                    'THED:  ', THED, &
!                    'RADUS: ', R_SVD,&
!                    'PITCH: ', PITCH

!57 FORMAT(A7,F8.4,2X,A7,F8.4,2X,A7,F8.4,2X,A7,F8.4,2X,A7,F8.4)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PRINT TRAJECTORY OF ROTATIONS, PDB FMT, WITH RES IN BETA COLUMN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
!
!   Make rotation trajectory: each rotation coordinate set is a separate frame
!    while RESidual information is being printed into the BETA column of PDB,
!    VMD does not appear to be displaying the colors....
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    COUNT = COUNT + 1  !Set each rotation as a frame in trajectory, PDB FMT
!!    WRITE(*,FMT=44) 'MODEL ', COUNT  !Print the MODEL Number
!!    !Write information for each atom in rotation frame
!!    DO IATOM = 1,NR
!!      WRITE(*,FMT=43)  'ATOM  ', IATOM, '  CA', 'ALA', '   1', &
!!                        newXYZ(IATOM,1), newXYZ(IATOM,2), newXYZ(IATOM,3),   &
!!                        '       ', RES
!!    END DO
    !Format for the ATOMic properties, including coordinates and RES as BETA 
!!43  FORMAT(A6,I5,A4,2X,A3,2X,A4,3X,F8.3,1X,F8.3,F8.3,A6,F6.2)
    !Append at the end of rotation frame ENDMDL to signify end of frame to VMD
!!    WRITE (*,FMT=45) 'ENDMDL'
    !Format for the MODEL number information
!!44  FORMAT(A6,1X,I4)
    !Format for the ENDMDL append information
!!45  FORMAT(A6)

!
!   Make rotation trajectory a single "molecule" in the PDB output:
!    proper color-coding display of RESidual information for each rotation
!    COMMENT OUT ABOVE TRAJECTORY WRITING FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    COUNT = COUNT + 1  !Set each rotation as a frame in trajectory, PDB FMT
!    DO IATOM = 1,NR
!      WRITE(*,FMT=43)  'ATOM  ', COUNT, '  CA', 'ALA', '   1', &
!                        newXYZ(IATOM,1), newXYZ(IATOM,2), newXYZ(IATOM,3),   &
!                        '       ', RES
!    END DO
!    !Format for the ATOMic properties, including coordinates and RES as BETA
!43  FORMAT(A6,I5,A4,2X,A3,2X,A4,3X,F8.3,1X,F8.3,F8.3,A6,F6.2)
!    !Append at the end of rotation frame ENDMDL to signify end of frame to VMD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Build the singular values matrices for plucking. There will be four, for each
!  singular value s1 >= s2 >= s3 > s0 = 0
   SING1(IPHI:IPHI,ITHETA:ITHETA) = SING(1)
   SING2(IPHI:IPHI,ITHETA:ITHETA) = SING(2)
   SING3(IPHI:IPHI,ITHETA:ITHETA) = SING(3)
   SING0(IPHI:IPHI,ITHETA:ITHETA) = SING(4)

!  Build the RESMAT array, indexed IPHI,ITHETA to keep track of which plane
!   gives the best-fit circle and thus represents the superhelical normal plane
   RESMAT(IPHI:IPHI,ITHETA:ITHETA)  = RES
   !Clear the scalar for the next rotation
   RES = 0.0D+00

!  Build the Pitch matrix to store all Pitch-s indx IPHI, ITHETA to pluck from
   PITMAT(IPHI:IPHI,ITHETA:ITHETA)  =  PITCH
   !Clear the scalar for the next rotation
   PITCH = 0.0D+00

!  Build the total sweep angle matrix.
   THETAPITMAT(IPHI:IPHI,ITHETA:ITHETA) =  THETAHELIX
   !Clear the scalar for the next rotation
   THETAHELIX = 0.0D+00
   THETASTEP  = 0.0D+00
   THEDSTEP   = 0.0D+00

!  Build the Radius matrix to store all the Radii indx IPHI, ITHETA to pluck..
   RSVDMAT(IPHI:IPHI,ITHETA:ITHETA) = R_SVD
   R_SVD = 0.0D+00

!  Build the PHID and THED matrices to store the angles(deg) of transformation
   PHIDMAT(IPHI:IPHI,ITHETA:ITHETA) = PHID
   THEDMAT(IPHI:IPHI,ITHETA:ITHETA) = THED
   !Clear scalars for next rotation
   PHID    = 0.0D+00
   THED    = 0.0D+00
   PHI_L   = 0.0D+00
   THETA_M = 0.0D+00

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! FINISH LOOP OVER UNIT SPHERE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                                   !
  END DO ! ITHETA_complete the loop of the phi-s   !
 END DO ! IPHI_complete the loop of the theta-s    !
                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   VAGELIS TEST RESIDUAL PLANE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    COUNT = 0

!    DO I    = 2,MAXITR-1   ! Loop-PHI
!     DO J   = 2,MAXITR-1   ! Loop-THETA

!     COUNT = COUNT + 1

     !See line 211 for index to angle converter...
!     VPHI   =  ( ( I * dPHI )   + PHIS ) * 180.0 / PI
!     VTHETA =  ( ( J * dTHETA ) + THES ) * 180.0 / PI

!     IF ( RESMAT(I,J) .LT. RESMAT(I+1,J) .AND. &
!          RESMAT(I,J) .LT. RESMAT(I-1,J) .AND. &
!          RESMAT(I,J) .LT. RESMAT(I,J+1) .AND. &
!          RESMAT(I,J) .LT. RESMAT(I,J-1) ) THEN

     !First column, values of PHI; Second column, THETA; Third column,
     !Residual at the values of PHI and THETA...
!       MINIMAT(COUNT,1) = VPHI
!       MINIMAT(COUNT,2) = VTHETA
!       MINIMAT(COUNT,3) = RESMAT(I,J)

!        WRITE(*,FMT=77) 'MINPHI:   ', VPHI, &
!                        'MINTHETA: ', VTHETA, &
!                        'MINRESID: ', RESMAT(I,J)

!!!!         WRITE(*,*) 'TEST OF IF STATEMENT'
!!!!     ELSE 
!!!!         WRITE(*,*) 'TEST OF ELSE STATEMENT'
!     END IF

!     END DO
!   END DO

!77 FORMAT(A10,2X,F8.4,2X,A10,F8.4,2X,A10,F8.4)


!   DO I = 1, COUNT
!     WRITE(*,*) 'MIN_PHI:   ', MINIMAT(I,1), &
!                'MIN_THETA: ', MINIMAT(I,2), &
!                'MIN_RES:   ', MINIMAT(I,3)
!   END DO



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! AFTER ALL ROTATIONS OVER SPHERE, FIND BEST ROTATION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !Find the combination of Phi and Theta (indices of *MAT below) that
   ! yielded the lowest residual, ie the best fit
   LOWESTRES   = MINLOC ( RESMAT , MASK = RESMAT .GT. 0.1D-11)

   !Set the index values as variables for which to use to select parameters
   ! corresponding to the "Best" Phi and Theta rotations.
   MIN_IPHI    = LOWESTRES(1)
   MIN_ITHETA  = LOWESTRES(2)
     ! Check what the index is, IPHI, ITHETA, of our residuals matrix
     !!WRITE(*,*)    'LOWEST RESIDUAL LOCATION IN RESMAT: ', LOWESTRES(1:2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! OBTAIN THE PARAMETERS AND PROPERTIES FROM THE PLANE OF BEST FIT
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   BEST_RES = RESMAT (MIN_IPHI,MIN_ITHETA)
   BEST_RAD = RSVDMAT(MIN_IPHI,MIN_ITHETA)
   BEST_PIT = PITMAT (MIN_IPHI,MIN_ITHETA)
   BEST_THETAPIT = THETAPITMAT(MIN_IPHI,MIN_ITHETA)
   BEST_PHI = PHIDMAT(MIN_IPHI,MIN_ITHETA)
   BEST_THE = THEDMAT(MIN_IPHI,MIN_ITHETA)
   BEST_S1  = SING1  (MIN_IPHI,MIN_ITHETA)
   BEST_S2  = SING2  (MIN_IPHI,MIN_ITHETA)
   BEST_S3  = SING3  (MIN_IPHI,MIN_ITHETA)
   BEST_S0  = SING0  (MIN_IPHI,MIN_ITHETA)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! WRITING CONTROL
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Print the best residual, BEST_RES, and the associated
!  radius, BEST_RAD, and pitch, BEST_PIT:
   WRITE(*,FMT=99) IFRAME, BEST_PHI, BEST_THE,BEST_RES,BEST_RAD,BEST_PIT,&
                   BEST_THETAPIT
99 FORMAT(1X,I6,2X,F8.4,2X,F8.4,2X,F11.6,2X,F11.6,2X,F11.6,2X,F11.6)


!   WRITE(*,*) 'BEST_S3:   BEST_S2:   BEST_S1:  BEST_S0:'
!   WRITE(*,FMT=98) BEST_S1, BEST_S2, BEST_S3, BEST_S0
!98 FORMAT(F8.4,2X,F8.4,2X,F8.4,2X,F8.4)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  ZERO-OUT ARRAYS AND VARIABLES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Clear variables that are summed within an IPHI ITHETA loop
   MIN_IPHI   = 0  !V -1
   MIN_ITHETA = 0  !V -1

!  Clear the arrays
   DO I = 1, MAXITR-2
      RESMAT(I,:)  = 0.0D+00
      RESMAT(:,I)  = 0.0D+00
      PITMAT(I,:)  = 0.0D+00
      PITMAT(:,I)  = 0.0D+00
      RSVDMAT(I,:) = 0.0D+00
      RSVDMAT(:,I) = 0.0D+00
      SING1  (:,I) = 0.0D+00
      SING1  (I,:) = 0.0D+00
      SING2  (:,I) = 0.0D+00
      SING2  (I,:) = 0.0D+00
      SING3  (:,I) = 0.0D+00
      SING3  (I,:) = 0.0D+00
      SING0  (:,I) = 0.0D+00
      SING0  (I,:) = 0.0D+00
   END DO

!  Clean up the coordinate arrays
   DO IATOM = 1,NR
      newX(IATOM)          = 0.0D+00
      newY(IATOM)          = 0.0D+00
      newZ(IATOM)          = 0.0D+00
      newXYZ_AMAT(IATOM,3) = 0.0D+00
      newXYZ_AMAT(IATOM,2) = 0.0D+00
      newXYZ_AMAT(IATOM,1) = 0.0D+00
      newXYZ_BMAT(IATOM)   = 0.0D+00
   END DO

!  End of the IFRAME loop, so we-re done with the trajectory!
END DO


! CPU TIMING
   CALL CPU_TIME(FINISH)
   PRINT '("Time = ",f10.3, " seconds")',finish-start

! End program..
END PROGRAM KSUPERHELIXT
