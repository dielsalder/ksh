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

   integer, parameter	:: n = 2, m = 5, iprint = 1
   integer, parameter	:: dp = kind(1.0d0)
   real(dp), parameter	:: factr = 1.d+1

   character(len=60)	:: task, csave
   logical		:: lsave(4)
   integer		:: isave(44)
   real(dp)             :: f
   real(dp)             :: dsave(29)
   integer,  allocatable:: nbd(:), iwa(:)
   real(dp), allocatable:: x(:), l(:), u(:), g(:), wa(:)

   allocate ( nbd(n), x(n), l(n), u(n), g(n)

   ! find out what this means
   allocate ( iwa(3*n) )
   allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

   l(1) = 0
   u(2) = 180

   do i = 1, n
	   x(i) = 0.0
   task = 'START'

   do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. & task.eq.'START'
	   call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
		       wa, iwa, task, iprint,&
		       csave, lsave, isave, dsave )

	   if (task(1:2) .eq. 'FG') then

	   f = CALL SVA (newXYZ_AMAT,NR,NR,3,NR,newXYZ_BMAT,SING, &
	                 KPVEC,NAMES,1,D,WORK)

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
