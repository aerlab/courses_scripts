

! **********************************************************************
! coefficient calculations for four first-order differential equations.
! **********************************************************************
	subroutine coeff1(ib, w, w1, w2, w3, t0, t1, u0, f0, &
	                   b, c)   
        implicit none
	include 'Gaussian-Legendre.h'
	include 'Rad4S.h'

! input 

       integer ib
       real w, w1, w2, w3, t0, t1, u0, f0

! output
       
       real, dimension(4,3):: b
       real, dimension(4,5):: c
 
! internal
       integer i, j
       real x, w0w, w1w, w2w, w3w
       real q1, q2, q3
       real fq, fw
       
	x = 0.5 * w
	w0w = x
	w1w = x * w1
	w2w = x * w2
	w3w = x * w3
	if ( ib .le. mbs ) then
	  fw = u0 * u0
	  q1 = - w1w * u0
	  q2 = w2w * ( 1.5 * fw - 0.5 )
	  q3 = - w3w * ( 2.5 * fw - 1.5 ) * u0
	endif
	fq = 0.5 * w0w
	do 10 i = 3, 4
	   do 20 j = 1, 4
	      c(i,j) = fq + w1w * p11d(i,j) + &
     	      w2w * p22d(i,j) + w3w * p33d(i,j) 
	      if ( i .eq. j ) then 
   	        c(i,j) = ( c(i,j) - 1.0 ) / u(i)
	      else
	        c(i,j) = c(i,j) / u(i)
	      endif
20	   continue
10	continue
	do 30 i = 1, 4
	   if ( ib .le. mbs ) then
	     c(i,5) = w0w + q1 * p1d(i) + &
     	     q2 * p2d(i) + q3 * p3d(i) 
	   else
	     c(i,5) = 1.0
           endif
	   c(i,5) = c(i,5) / u(i)
30	continue
	b(1,1) = c(4,4) - c(4,1)
	b(1,2) = c(4,4) + c(4,1)
	b(2,1) = c(4,3) - c(4,2)
	b(2,2) = c(4,3) + c(4,2)
	b(3,1) = c(3,4) - c(3,1)
	b(3,2) = c(3,4) + c(3,1)
	b(4,1) = c(3,3) - c(3,2)
	b(4,2) = c(3,3) + c(3,2)
	b(1,3) = c(4,5) - c(1,5)
	b(2,3) = c(3,5) - c(2,5)
	b(3,3) = c(3,5) + c(2,5)
	b(4,3) = c(4,5) + c(1,5)
	return
	end

! **********************************************************************
! coefficient calculations for second order differential equations.
! **********************************************************************
	subroutine coeff2(ib, w, w1, w2, w3, t0, t1, u0, f0, b, &
	                   a, d)
	implicit none
	include 'Rad4S.h'

! input
       real w, w1, w2, w3, t0, t1, u0, f0
       integer ib
       real, dimension(4,3)::b

! output
       real, dimension(2,2,2)::a
       real, dimension(4):: d
 
! internal
       real fw1, fw2, fw3, fw4
              
              	 
!	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
!	common /coedf1/ b(4,3)
!	common /coedf2/ a(2,2,2), d(4)

	fw1 = b(1,1) * b(1,2)
	fw2 = b(2,1) * b(3,2)
	fw3 = b(3,1) * b(2,2)
	fw4 = b(4,1) * b(4,2)
	a(2,2,1) = fw1 + fw2
	a(2,1,1) = b(1,1) * b(2,2) + b(2,1) * b(4,2)
	a(1,2,1) = b(3,1) * b(1,2) + b(4,1) * b(3,2)
	a(1,1,1) = fw3 + fw4
	a(2,2,2) = fw1 + fw3
	a(2,1,2) = b(1,2) * b(2,1) + b(2,2) * b(4,1)
	a(1,2,2) = b(3,2) * b(1,1) + b(4,2) * b(3,1)
	a(1,1,2) = fw2 + fw4
	d(1) = b(3,2) * b(4,3) + b(4,2) * b(3,3) + b(2,3) / u0
	d(2) = b(1,2) * b(4,3) + b(2,2) * b(3,3) + b(1,3) / u0
	d(3) = b(3,1) * b(1,3) + b(4,1) * b(2,3) + b(3,3) / u0
	d(4) = b(1,1) * b(1,3) + b(2,1) * b(2,3) + b(4,3) / u0
	return
	end

! **********************************************************************
! coefficient calculations for fourth-order differential equations.
! **********************************************************************
	subroutine coeff4 (ib, w, w1, w2, w3, t0, t1, u0, f0, a, d, &
	                   b1, c1, z)
	implicit none

! input
        integer ib
	real 	w, w1, w2, w3, t0, t1, u0, f0
        real, dimension(2,2,2)::a
        real, dimension(4):: d

! output
        real, dimension(4):: z
	real b1, c1


! internal
        real x
	
!	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
!	common /coedf2/ a(2,2,2), d(4)
!	common /coedf4/ b1, c1, z(4)
	

	x = u0 * u0
	b1 = a(2,2,1) + a(1,1,1)
	c1 = a(2,1,1) * a(1,2,1) - a(1,1,1) * a(2,2,1)
	z(1) = a(2,1,1) * d(3) + d(4) / x - a(1,1,1) * d(4)
	z(2) = a(1,2,1) * d(4) - a(2,2,1) *d(3) + d(3) / x
	z(3) = a(2,1,2) * d(1) + d(2) / x - a(1,1,2) * d(2)
	z(4) = a(1,2,2) * d(2) - a(2,2,2) * d(1) + d(1) / x
	return
	end
	
	
! **********************************************************************
! fk1 and fk2 are the eigenvalues.
! **********************************************************************
	subroutine coeffl ( ib, w, w1, w2, w3, t0, t1, u0, f0, &
	                     b, a, d, b1, c1, z, &
			     aa, zz, a1, z1, fk1, fk2)
	
	implicit none
	include 'Rad4S.h'

! input
        real    w, w1, w2, w3, t0, t1, u0, f0
	integer ib
	real b(4,3), a(2,2,2), d(4), b1, c1, z(4)
	
!output
        real 	aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
	 
!internal
        real dt, x, fw, y, zx, fq0, fq1 , a2, b2 , fw1, fw2 
        integer i

!        use RadParams 
!	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
!	common /coedf1/ b(4,3)
!	common /coedf2/ a(2,2,2), d(4)
!	common /coedf4/ b1, c1, z(4)
!	common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
	
	
	dt = t1 - t0
	x = sqrt ( b1 * b1 + 4.0 * c1 )
	fk1 = sqrt ( ( b1 + x ) * 0.5 )
	fk2 = sqrt ( ( b1 - x ) * 0.5 )
	fw = u0 * u0
	x = 1.0 / ( fw * fw ) - b1 / fw - c1
!by JW	
	if ( abs(x) .lt. 1.0E-20) then
	  print*, 'JW: x is too small, change from', x 
	  if ( x .lt. 0 ) then 
	   x = -1.0e-20
	  else
	   x = 1.0e-20
	  endif 
          print*, 'JW: TTTTTTO x', x 
	endif  
!by JW

	fw = 0.5 * f0 / x
	z(1) = fw * z(1) 
	z(2) = fw * z(2) 
	z(3) = fw * z(3) 
	z(4) = fw * z(4)
	z1(1) = 0.5 * ( z(1) + z(3) )
	z1(2) = 0.5 * ( z(2) + z(4) )
	z1(3) = 0.5 * ( z(2) - z(4) )
	z1(4) = 0.5 * ( z(1) - z(3) )
	a2 = ( fk1 * fk1 - a(2,2,1) ) / a(2,1,1)
	b2 = ( fk2 * fk2 - a(2,2,1) ) / a(2,1,1)
	x = b(1,1) * b(4,1) - b(3,1) * b(2,1)
	fw1 = fk1 / x
	fw2 = fk2 / x
	y = fw2 * ( b2 * b(2,1) - b(4,1) ) 
	zx = fw1 * ( a2 * b(2,1) - b(4,1) )
	a1(1,1) = 0.5 * ( 1 - y )
   	a1(1,2) = 0.5 * ( 1 - zx )
	a1(1,3) = 0.5 * ( 1 + zx )
	a1(1,4) = 0.5 * ( 1 + y )
	y = fw2 * ( b(3,1) - b2 * b(1,1) ) 
	zx = fw1 * ( b(3,1) - a2 * b(1,1) ) 
	a1(2,1) = 0.5 * ( b2 - y )
	a1(2,2) = 0.5 * ( a2 - zx )
	a1(2,3) = 0.5 * ( a2 + zx )
	a1(2,4) = 0.5 * ( b2 + y )
        a1(3,1) = a1(2,4)
        a1(3,2) = a1(2,3)
        a1(3,3) = a1(2,2)
        a1(3,4) = a1(2,1)
        a1(4,1) = a1(1,4)
        a1(4,2) = a1(1,3)
        a1(4,3) = a1(1,2)
        a1(4,4) = a1(1,1)
	if ( ib .le. mbs ) then
	  fq0 = exp ( - t0 / u0 )
          fq1 = exp ( - t1 / u0 )
	else
	  fq0 = 1.0
	  fq1 = exp ( - dt / u0 )
	endif
	x = exp ( - fk1 * dt )
	y = exp ( - fk2 * dt )
	do 40 i = 1, 4
	   zz(i,1) = z1(i) * fq0
	   zz(i,2) = z1(i) * fq1
	   aa(i,1,1) = a1(i,1)
	   aa(i,2,1) = a1(i,2)
	   aa(i,3,1) = a1(i,3) * x
	   aa(i,4,1) = a1(i,4) * y
	   aa(i,3,2) = a1(i,3)
	   aa(i,4,2) = a1(i,4)
	   aa(i,1,2) = a1(i,1) * y
	   aa(i,2,2) = a1(i,2) * x
40	continue
	return
	end

! **********************************************************************
! See the paper by Liou, Fu and Ackerman (1988) for the formulation of
! the delta-four-stream approximation in a homogeneous layer.
! **********************************************************************
	subroutine coefft(ib, w, w1, w2, w3, t0, t1, u0, f0, &
	                  aa, zz, a1, z1, fk1, fk2) 
        implicit none

! input 
        integer ib
	real w, w1, w2, w3, t0, t1, u0, f0

!output
        real aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2
	 

!internal
        real b(4,3), a(2,2,2), d(4), b1, c1, z(4)
	real, dimension(4,5):: c

	call coeff1 (ib, w, w1, w2, w3, t0, t1, u0, f0, &
	                   b, c) 
	call coeff2 (ib, w, w1, w2, w3, t0, t1, u0, f0, b, &
	                   a, d) 
	call coeff4 (ib, w, w1, w2, w3, t0, t1, u0, f0, a, d, &
	                   b1, c1, z)
	call coeffl ( ib, w, w1, w2, w3, t0, t1, u0, f0, &
	                     b, a, d, b1, c1, z, &
			     aa, zz, a1, z1, fk1, fk2)
	return
	end
	
	
	
! **********************************************************************
! In the limits of no scattering ( Fu, 1991 ), fk1 = 1.0 / u(3) and
! fk2 = 1.0 / u(4).
! **********************************************************************
	subroutine coefft0 (ib, w, w1, w2, w3, t0, t1, u0, f0, &
	                    aa, zz, a1, z1, fk1, fk2)
	
	implicit none
	include 'Rad4S.h'
	include 'Gaussian-Legendre.h'

! input
 	integer ib
	real w, w1, w2, w3, t0, t1, u0, f0

!output
        real aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2


! internal
        real  wjuu0, y, dt, x
	real fw
	integer i, j, jj, k

	
!	common /point/ u(4)
!	common /coedfi/ ib, w, w1, w2, w3, t0, t1, u0, f0
!	common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2

	wjuu0 = 0.0

	fk1 = 4.7320545
	fk2 = 1.2679491
	y = exp ( - ( t1 - t0 ) / u0 )
	fw = 0.5 * f0
	do 10 i = 1, 4
	   if ( ib .le. mbs ) then
             z1(i) = 0.0
             zz(i,1) = 0.0
             zz(i,2) = 0.0
           else
	     jj = 5 - i
!by JW             
	     wjuu0 =  u(jj) / u0 

	     if ( wjuu0 .eq. -1. ) then
		if (u(jj) .gt. 0 ) then 
		 wjuu0 = -1.0+1.0E-6 
		else
		 wjuu0 = -1.- 1.0E-6 
		endif
                write(*,111), 'JW: wjju0 changes from',  u(jj) , u0, 'to ', wjuu0  
111                format(1x, a, 1x, e17.10, e17.10, a, e17.10, a, I4, a, I4)    
!	      stop
	     endif	
	     
!	     z1(i) = fw / ( 1.0 + u(jj) / u0 )
	     z1(i) = fw / ( 1.0 + wjuu0 )

!endofchange	     
	     zz(i,1) = z1(i) 
	     zz(i,2) = z1(i) * y
	   endif
	do 11 j = 1, 4
	   a1(i,j) = 0.0
	do 12 k = 1, 2
	   aa(i,j,k) = 0.0
12	continue
11      continue
10      continue

	do 20 i = 1, 4
	   j = 5 - i
	   a1(i,j) = 1.0
20	continue
	dt = t1 - t0
	x = exp ( - fk1 * dt )
	y = exp ( - fk2 * dt )
	aa(1,4,1) = y
	aa(2,3,1) = x
	aa(3,2,1) = 1.0
	aa(4,1,1) = 1.0
	aa(1,4,2) = 1.0
	aa(2,3,2) = 1.0
	aa(3,2,2) = x
	aa(4,1,2) = y
	return
	end
	
	

! **********************************************************************
! In the solar band  asbs is the surface albedo, while in the infrared
! band asbs is  blackbody intensity emitted at the surface temperature
! times surface emissivity.  In this subroutine, the delta-four-stream
! is applied to nonhomogeneous atmospheres. See comments in subroutine
! 'qcfel' for array AB(13,4*n).
! **********************************************************************
!        use RadParams 
!	common /dis/ a(4)
!	common /point/ u(4)
!	common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx),& 
!    	                t(ndfsx), u0(ndfsx), f0(ndfsx) 
!	common /coedfi/ ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
!	common /coedfl/ aa(4,4,2), zz(4,2), a1(4,4), z1(4),&
!                       fk1t, fk2t
!	common /qccfeo/ fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx), &
!     	                z4(4,ndfsx), g4(4,ndfsx)
!	common /qcfelc/ ab(13,ndfs4x), bx(ndfs4x), xx(ndfs4x)
!	dimension fu(4,4), wu(4)



	subroutine qccfe ( ndfs, ndfs4, ib, asbs, ee, &
	                   w1, w2, w3, w, t, u0, f0 , &
			   fk1, fk2, a4, z4, g4)
        implicit none
	include 'Rad4S.h'
	include 'Gaussian-Legendre.h'
	
!input
        integer ib, ndfs, ndfs4
	real asbs
	real ee
	real, dimension(ndfsx):: w1, w2, w3, w, t, u0, f0
	

!output
	real, dimension(ndfsx):: fk1, fk2
	real a4(4, 4, ndfsx), z4(4, ndfsx), g4(4, ndfsx)
	

!internal
	real fu(4,4), wu(4)
	integer ibn
	real wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n
	real aa(4, 4, 2), zz(4, 2), a1(4, 4), z1(4), fk1t, fk2t
	integer n, n4, i, j, k, kf, i1, i2, j1, j2, i3, j3, i8
	integer m1, m2, m18, m28
	real fw1, fw2
	real v1, v2, v3
	real ab(13, ndfs4x), bx(ndfs4x), xx(ndfs4x)

	
! starts
	n = ndfs
	n4 = ndfs4
	do 333 i = 1, n4
	do 333 j = 1, 13
	   ab(j,i) = 0.0
333	continue

!  JW calculate coefficients at the top layer 
	ibn = ib
	wn = w(1)
	w1n = w1(1)
	w2n = w2(1)
	w3n = w3(1)
	t0n = 0.0
	t1n = t(1)
	u0n = u0(1)
	f0n = f0(1)

	if ( wn .ge. 0.999999 ) then
          wn = 0.999999
        endif
	if ( wn .le. 1.0e-4 ) then

 	  call coefft0 (ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n, &
	                    aa, zz, a1, z1, fk1t, fk2t)
	  fk1(1) = fk1t
	  fk2(1) = fk2t
        else
       	  call coefft( ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n, &
	                  aa, zz, a1, z1, fk1t, fk2t) 
	  fk1(1) = fk1t
	  fk2(1) = fk2t
	endif

! continue to calculate coefficients  
	do 10 i = 1, 4
	   z4(i,1) = z1(i)
	do 10 j = 1, 4
	   a4(i,j,1) = a1(i,j)
10	continue
	do 20 i = 1, 2
	   bx(i) = - zz(i+2,1)
           i8 = i + 8
	do 20 j = 1, 4
	   ab(i8-j,j) = aa(i+2,j,1)
20	continue
	do 30 i = 1, 4
	   wu(i) = zz(i,2)
	do 30 j = 1, 4
	   fu(i,j) = aa(i,j,2)
30	continue

! calculate for other layers. 
	do 40 k = 2, n
	   wn = w(k)
	   w1n = w1(k)
	   w2n = w2(k)
	   w3n = w3(k)
	   t0n = t(k-1)
	   t1n = t(k)
	   u0n = u0(k)
	   f0n = f0(k)
	   if ( wn .ge. 0.999999 ) then
             wn = 0.999999
           endif
	   if ( wn .le. 1.0e-4 ) then
 
	     call coefft0 (ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n, &
	                    aa, zz, a1, z1, fk1t, fk2t)
	     fk1(k) = fk1t
	     fk2(k) = fk2t
	   else
	     call coefft (ibn, wn, w1n, w2n, w3n, t0n, t1n, u0n, f0n, &
	                  aa, zz, a1, z1, fk1t, fk2t) 
	     fk1(k) = fk1t
	     fk2(k) = fk2t
           endif
	   do 50 i = 1, 4
	      z4(i,k) = z1(i)
	   do 50 j = 1, 4
	      a4(i,j,k) = a1(i,j)
50	   continue
	   kf = k + k + k + k
	   i1 = kf - 5
	   i2 = i1 + 3
	   j1 = kf - 7
	   j2 = j1 + 3
	   i3 = 0
	   do 55 i = i1, i2
	      i3 = i3 + 1
	      bx(i) = - wu(i3) + zz(i3,1)
	      j3 = 0
              i8 = i + 8
	      do 60 j = j1, j2
	         j3 = j3 + 1
	         ab(i8-j,j) = fu(i3,j3)
60	      continue
	      j3 = 0
	      do 65 j = j2 + 1, j2 + 4
	         j3 = j3 + 1
	         ab(i8-j,j) = - aa(i3,j3,1)
65	      continue
55	   continue
	   do 70 i = 1, 4
	      wu(i) = zz(i,2)
	   do 70 j = 1, 4
	      fu(i,j) = aa(i,j,2)
70	   continue
40	continue
	if ( ib .le. mbs ) then
	  v1 = 0.2113247 * asbs
	  v2 = 0.7886753 * asbs
	  v3 = asbs * u0(1) * f0(1) * exp ( - t(n) / u0(1) )
	  m1 = n4 - 1
	  m2 = n4
          m18 = m1 + 8
          m28 = m2 + 8
	  fw1 = v1 * wu(3)
	  fw2 = v2 * wu(4)
	  bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
	  bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
	  do 80 j = 1, 4
	     j1 = n4 - 4 + j
	     fw1 = v1 * fu(3,j)
	     fw2 = v2 * fu(4,j)
	     ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
	     ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
80	  continue
        else
	  v1 = 0.2113247 * ( 1.0 - ee )
	  v2 = 0.7886753 * ( 1.0 - ee )
	  v3 = asbs
	  m1 = n4 - 1
	  m2 = n4
          m18 = m1 + 8
          m28 = m2 + 8
	  fw1 = v1 * wu(3)
	  fw2 = v2 * wu(4)
	  bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
	  bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
	  do 85 j = 1, 4
	     j1 = n4 - 4 + j
	     fw1 = v1 * fu(3,j)
	     fw2 = v2 * fu(4,j)
	     ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
	     ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
85	  continue
	endif
	call qcfel (ndfs, ndfs4, ab, bx, xx)
	do 90 k = 1, n
	   j = k + k + k + k - 4
	do 90 i = 1, 4
	   j = j + 1
	   g4(i,k) = xx(j)
90	continue
	return
	end

! **********************************************************************
	subroutine qcfel (ndfs, ndfs4, ab, b, x)
! **********************************************************************
! 1. `qcfel' is the abbreviation of ` qiu constants for each layer'.
! 2. The inhomogeneous atmosphere is divided into n adjacent homogeneous
!    layers where the  single scattering properties are constant in each
!    layer and allowed to vary from one to another. Delta-four-stream is
!    employed for each homogeneous layer. The boundary conditions at the
!    top and bottom of the atmosphere,  together with  continuity condi-
!    tions  at  layer interfaces lead to a system of algebraic equations
!    from which 4*n unknown constants in the problom can be solved.
! 3. This subroutine is used for solving the 4*n unknowns of A *X = B by
!    considering the fact that the coefficient matrix is a sparse matrix
!    with the precise pattern in this special problom.
! 4. The method is not different in principle from the general scheme of
!    Gaussian elimination with backsubstitution, but carefully optimized
!    so as to minimize arithmetic operations.  Partial  pivoting is used
!    to quarantee  method's numerical stability,  which will  not change
!    the basic pattern of sparsity of the matrix.
! 5. Scaling special problems so as to make  its nonzero matrix elements
!    have comparable magnitudes, which will ameliorate the stability.
! 6. a, b and x present A, B and X in A*X=B, respectively. and n4=4*n.
! 7. AB(13,4*n) is the matrix A in band storage, in rows 3 to 13; rows 1
!    and 2 and other unset elements should be set to zero on entry.
! 8. The jth column of A is stored in the jth column of the array AB  as
!    follows:
!            AB(8+i-j,j) = A(i,j) for max(1,j-5) <= i <= min(4*n,j+5).
!    Reversedly, we have
!            A(ii+jj-8,jj) = AB(ii,jj).
! **********************************************************************
!common /qcfelc/ ab(13,ndfs4x), b(ndfs4x), x(ndfs4x)

        implicit none
	include 'Rad4S.h'
	
!input    
        integer ndfs, ndfs4

! output and input 
         real 	ab(13,ndfs4x), b(ndfs4x), x(ndfs4x)
        
	 

!internal
         integer n, n4, n44, i, j, m
	 integer k, l, m1, k44, m18, i0m1, i0, m1f, i0f, im1, ifq 
	 integer n1, n2, n3, m2, m3, m4, m28, m38, m48
	 real p, yy, xx
	 real t
	 	 	
! starts
	
	n = ndfs
	n4 = ndfs4
	do 5 k = 1, n - 1
           k44 = 4 * k - 4
	   do 3 l= 1, 4
	      m1 = k44 + l
	      p = 0.0
	      do 10 i = 8, 14 - l
	         if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
	           p = ab(i,m1)
	           i0 = i
                 endif
10	      continue
              i0m1 = i0 + m1
              m18 = m1 + 8
	      if ( i0 .eq. 8 ) goto 20
	      do 15 j = m1, m1 + 8 - l
                 i0f = i0m1 - j
                 m1f = m18 - j
	         t = ab(i0f,j)
	         ab(i0f,j) = ab(m1f,j)
  	         ab(m1f,j) = t
15            continue
              i0f = i0m1 - 8
	      t = b(i0f)
	      b(i0f) = b(m1)
	      b(m1) = t
20	      continue
	      yy = ab(8,m1)
	      ab(8,m1) = 1.0
	      do 25 j = m1 + 1, m1 + 8 - l
                 m1f = m18 - j
	         ab(m1f,j) = ab(m1f,j) / yy
25	      continue
	      b(m1) = b(m1) / yy
	      do 30 i = 9, 14 - l
	         xx = ab(i,m1)
                 ab(i,m1) = 0.0
                 im1 = i + m1
	         do 35 j = m1 + 1, m1 + 8 - l
                    ifq = im1 - j
                    m1f = m18 - j
	            ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
35	         continue
	         ifq = im1 - 8
	         b(ifq) = b(ifq) - b(m1) * xx
30	      continue
3	   continue
5	continue
	n44 = n4 - 4
	do 40 l = 1, 3
	   m1 = n44 + l
	   p = 0.0
	   do 45 i = 8, 12 - l
	      if ( abs ( ab(i,m1) ) .gt. abs ( p ) ) then
	        p = ab(i,m1)
	        i0 = i
              endif
45	   continue
           i0m1 = i0 + m1
           m18 = m1 + 8

!JW	   if( i0 .eq. 8 ) goto 55
	 
	   if ( i0 .ne. 8 ) then  
	   do 50 j = m1, m1 + 4 - l
              i0f = i0m1 - j
              m1f = m18 - j
	      t = ab(i0f,j)
	      ab(i0f,j) = ab(m1f,j)
  	      ab(m1f,j) = t
50	   continue
           i0f = i0m1 - 8
	   t = b(i0f)
	   b(i0f) = b(m1)
	   b(m1) = t
	   endif
	   
!55	   continue



	   yy = ab(8,m1)
           ab(8,m1) = 1.0
	   do 60 j = m1 + 1, m1 + 4 - l
              m1f = m18 - j
	      ab(m1f,j) = ab(m1f,j) / yy
60	   continue
	   b(m1) = b(m1) / yy
	   do 65 i = 9, 12 - l
	      xx = ab(i,m1)
              ab(i,m1) = 0.0
              im1 = i + m1
	      do 70 j = m1 + 1, m1 + 4 - l
                 ifq = im1 - j
                 m1f = m18 - j
	         ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
70	      continue
              ifq = im1 - 8
	      b(ifq) = b(ifq) - b(m1) * xx
65	   continue
40	continue
	yy = ab(8,n4)
	ab(8,n4) = 1.0
	b(n4) = b(n4) / yy
	n3 = n4 - 1
	n2 = n3 - 1
	n1 = n2 - 1
	x(n4) = b(n4)
	x(n3) = b(n3) - ab(7,n4) * x(n4)
	x(n2) = b(n2) - ab(7,n3) * x(n3) - ab(6,n4) * x(n4)
	x(n1) = b(n1) - ab(7,n2) * x(n2) - ab(6,n3) * x(n3) -&
     	        ab(5,n4) * x(n4)
	do 80 k = 1, n - 1
	   m4 = 4 * ( n - k )
	   m3 = m4 - 1
	   m2 = m3 - 1
	   m1 = m2 - 1
           m48 = m4 + 8
           m38 = m3 + 8
           m28 = m2 + 8
           m18 = m1 + 8
	   x(m4) = b(m4)
	   do 85 m = m4 + 1, m4 + 4
  	      x(m4) = x(m4) - ab(m48-m,m) * x(m)
85	   continue
	   x(m3) = b(m3)
	   do 90 m = m3 + 1, m3 + 5
  	      x(m3) = x(m3) - ab(m38-m,m) * x(m)
90	   continue
	   x(m2) = b(m2)
	   do 95 m = m2 + 1, m2 + 6
  	      x(m2) = x(m2) - ab(m28-m,m) * x(m)
95	   continue
	   x(m1) = b(m1)
	   do 100 m = m1 + 1, m1 + 7
   	      x(m1) = x(m1) - ab(m18-m,m) * x(m)
100	   continue
80	continue
	return
	end




! **********************************************************************
! In this subroutine, we incorporate a delta-function adjustment to
! account for the  forward  diffraction  peak in the context of the 
! four-stream approximation ( Liou, Fu and Ackerman, 1988 ). w1(n),
! w2(n), w3(n), w(n), and t(n) are the adjusted parameters.
! **********************************************************************
!	common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),&
!     	               ww(ndfsx), tt(ndfsx)
!	common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), &
!     	                t(ndfsx), u0a(ndfsx), f0a(ndfsx)


	subroutine adjust (ndfs, ww1, ww2, ww3, ww4, ww, tt, &
	                        w1, w2, w3, w, t)
	 
        implicit none
	include 'Rad4S.h'

! input 
        integer ndfs 
        real, dimension(ndfsx):: ww1, ww2, ww3, ww4, ww, tt

! output
        real, dimension(ndfsx):: w1, w2, w3, w, t 	
		       
! internal
        real, dimension(ndfsx) :: dtt, dt
	integer i, n
	real tt0, f, fw
	
	
	n = ndfs
	tt0 = 0.0
	do 10 i = 1, n
	   f = ww4(i) / 9.0
	   fw = 1.0 - f * ww(i) 
	   w1(i) = ( ww1(i) - 3.0 * f ) / ( 1.0 - f )
	   w2(i) = ( ww2(i) - 5.0 * f ) / ( 1.0 - f )
	   w3(i) = ( ww3(i) - 7.0 * f ) / ( 1.0 - f )
	   w(i) = ( 1.0 - f ) * ww(i) / fw
	   dtt(i) = tt(i) - tt0
	   tt0 = tt(i)
	   dt(i) = dtt(i) * fw
10	continue
	t(1) = dt(1)
	do 20 i = 2, n
	   t(i) = dt(i) + t(i-1)
20	continue
	return
	end

! **********************************************************************
! The delta-four-stream approximation for nonhomogeneous atmospheres
! in the solar wavelengths (Fu, 1991). The input parameters are ndfs,
! mdfs, and ndfs4 through 'param.ipt',  ib, as, u0, f0 for solar and
! ib, bf, bs, ee for IR through arguments of  'qfts' and 'qfti', and
! ww1(ndfs), ww2(ndfs), ww3(ndfs), ww4(ndfs), ww(ndfs), and tt(ndfs)
! through common statement 'dfsin'.
! **********************************************************************
!	common /dis/ a(4)
!	common /point/ u(4)
!	common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),&
!     	               ww(ndfsx), tt(ndfsx)
!	common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), &
!     	                t(ndfsx), u0a(ndfsx), f0a(ndfsx)
!	common /qccfeo/ fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx),& 
!     	                z4(4,ndfsx), g4(4,ndfsx)
!	common /dfsout/ ffu(mdfsx), ffd(mdfsx)
!-------------------------------------------------------------------------	
	
	subroutine qfts ( ndfs, mdfs, ib, as, u0, f0,  &
	                  ww1, ww2, ww3, ww4, ww, tt, &
			  ffu, ffd )
	
	implicit none
	include 'Rad4S.h'
	include 'Gaussian-Legendre.h'

! input
        integer ib, ndfs, mdfs
	real as
	real u0, f0
	real, dimension(ndfsx):: ww1, ww2, ww3, ww4, ww, tt

! output
        real, dimension(mdfsx):: ffu, ffd

! internal	
	real x(4), fi(4)
	integer n, m, k, jj, i , ndfs4, ii
	real y1 
	real asbs, ee, fw1, fw2, fw3, y, fw4
	real, dimension(ndfsx):: w1, w2, w3, w, t, u0a, f0a
	real, dimension(ndfsx):: fk1, fk2
	real a4(4, 4, ndfsx), z4(4, ndfsx), g4(4, ndfsx)

! initial	
	n = ndfs
	m = mdfs
	ndfs4 = ndfs*4
        ee = 0.0
	asbs = as
	
	
! first adjust	
	
	call adjust (n, ww1, ww2, ww3, ww4, ww, tt, &
	                        w1, w2, w3, w, t)
	
	do 5 i = 1, n
	   u0a(i) = u0
	   f0a(i) = f0
5	continue
       

! second calcualte coefficients

	call qccfe ( ndfs, ndfs4, ib, asbs, ee, &
	                   w1, w2, w3, w, t, u0a, f0a , &
			   fk1, fk2, a4, z4, g4)
        
	fw1 = 0.6638961
	fw2 = 2.4776962
        fw3 = u0 * 3.14159 * f0 
	do 10 i = 1, m
	   if ( i .eq. 1 ) then
             x(1) = 1.0
	     x(2) = 1.0
	     x(3) = exp ( - fk1(1) * t(1) )
	     x(4) = exp ( - fk2(1) * t(1) )
             k = 1
	     y = 1.0
	   elseif ( i .eq. 2 ) then
             x(1) = exp ( - fk2(1) * t(1) )
	     x(2) = exp ( - fk1(1) * t(1) )
	     x(3) = 1.0
	     x(4) = 1.0
	     k = 1
	     y = exp ( - t(1) / u0 )
	   else
	     k = i - 1
	     y1 = t(k) - t(k-1)
             x(1) = exp ( - fk2(k) * y1 )
             x(2) = exp ( - fk1(k) * y1 )
	     x(3) = 1.0
	     x(4) = 1.0
	     y = exp ( - t(k) / u0 )
	   endif
	   do 37 jj = 1, 4
	      fi(jj) = z4(jj,k) * y
37	   continue
	   do 40 ii = 1, 4
	      fw4 = g4(ii,k) * x(ii)
	      do 45 jj = 1, 4
	         fi(jj) = fi(jj) + a4(jj,ii,k) * fw4
45            continue
40	   continue
	   ffu(i)= fw1 * fi(2) + fw2 * fi(1) 
	   ffd(i)= fw1 * fi(3) + fw2 * fi(4) + fw3 * y
10	continue
	return
	end




! **********************************************************************
! The exponential approximation for the Planck function in optical depth
! is used for the infrared ( Fu, 1991). Since the direct solar radiation
! source has an exponential function form in terms of optical depth, the
! formulation of the delta-four-stream approximation for infrared  wave-
! lengths is the same as that for solar wavelengths. 
! **********************************************************************
!	common /dis/ a(4)
!	common /point/ u(4) 
!	common /dfsin/ ww1(ndfsx), ww2(ndfsx), ww3(ndfsx), ww4(ndfsx),&
!                    ww(ndfsx), tt(ndfsx)
!	common /qccfei/ w1(ndfsx), w2(ndfsx), w3(ndfsx), w(ndfsx), &
!     	                t(ndfsx), u0(ndfsx), f0(ndfsx)
!	common /qccfeo/ fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx), &
!     	                z4(4,ndfsx), g4(4,ndfsx)
!	common /dfsout/ ffu(mdfsx), ffd(mdfsx)
!	common /planci/ bf(nv1x), bs
!	common /grids/ ix, iy!


	subroutine qfti ( ndfs, mdfs, ib, ee, &
	                  bs, bf, ww1, ww2, ww3, ww4, ww, tt, &
			  ffu, ffd)

	implicit none
	include 'Rad4S.h'
	include 'Gaussian-Legendre.h'


!input 
        integer ib, ndfs, mdfs
	real ee 
	real bs
	real, dimension(ndfsx):: ww1, ww2, ww3, ww4, ww, tt
	real bf(nv1x)	

!output
        real, dimension(mdfsx):: ffu, ffd

!internal
        
	integer n, m, i, k, ndfs4, jj, ii
        real x(4), fi(4) 
	real asbs
	real, dimension(ndfsx):: w1, w2, w3, w, t, u0, f0
	real t0
	real timt0, q2, y1
	real fk1(ndfsx), fk2(ndfsx), a4(4,4,ndfsx), &
     	                z4(4,ndfsx), g4(4,ndfsx)
	real fw1, fw2, xy, fw3		
        real q1	

! initial
	
	n = ndfs
	m = mdfs
	asbs = bs * ee
	ndfs4 = 4 * ndfs


! adjust
	call adjust(ndfs, ww1, ww2, ww3, ww4, ww, tt,  &
	                        w1, w2, w3, w, t)
		 
	t0 = 0.0
	do 3 i = 1, n
	   q1 = alog ( bf(i+1) / bf(i) )
	   
!followd by Paul in later FUcodes and changed by Jun 	   
	   timt0 = t(i) - t0
	   if ( timt0 .lt. 1.0e-12) then
	     !print*, 'JWtimt0 is too small, change from ', t(i) - t0, 'to',  timt0 
	     !print*, 'adjust t ', t(i), tt(i)
	     timt0 = 1.0e-12
           endif
	     
	   q2 = 1.0/timt0
	   
	   
! change ends	   
!	   q2 = 1.0 / ( t(i) - t0 )
	   
	   f0(i) = 2.0 * ( 1.0 - w(i) ) * bf(i)
           if (abs(q1) .le. 1.0e-10) then
             u0(i) = -1.0e+10 / q2
           else
	     u0(i) = - 1.0 / ( q1 * q2 )
           end if
	   t0 = t(i)
!JW       
!       if ( ix .eq. 41 .and. iy .eq. 42 ) then
!        print*, 'q1  =', q1, 'q2 = ', q2 
!        print*, 'JWbefore call qccfe: u0= ', u0(i) 
!        print*, 'JWbefore call qccfe: bf= ',  bf(i+1), bf(i) 
!        print*, 'JWbefore call qccfe: t= ',  t(i) 
!       endif


3	continue

        call qccfe (ndfs, ndfs4, ib, asbs, ee, &
	            w1, w2, w3, w, t, u0, f0 , &
			fk1, fk2, a4, z4, g4     )

	fw1 = 0.6638958
	fw2 = 2.4776962
	do 10 i = 1, m
	   if ( i .eq. 1 ) then
             x(1) = 1.0
	     x(2) = 1.0
	     x(3) = exp ( - fk1(1) * t(1) )
	     x(4) = exp ( - fk2(1) * t(1) )
             k = 1
             xy = 1.0
	   elseif ( i .eq. 2 ) then
             x(1) = exp ( - fk2(1) * t(1) )
	     x(2) = exp ( - fk1(1) * t(1) )
	     x(3) = 1.0
	     x(4) = 1.0
             k = 1
	     xy =  exp ( - t(1) / u0(1) )
	   else
             k = i - 1
	     y1 = t(k) - t(k-1)
             x(1) = exp ( - fk2(k) * y1 )
             x(2) = exp ( - fk1(k) * y1 )
	     x(3) = 1.0
	     x(4) = 1.0
	     xy =  exp ( - y1 / u0(k) )
	   endif
	   do 37 jj = 1, 4
	      fi(jj) = z4(jj,k) * xy
37	   continue
           do 40 ii = 1, 4
	      fw3 = g4(ii,k) * x(ii)
	      do 45 jj = 1, 4
	         fi(jj) = fi(jj) + a4(jj,ii,k) * fw3
45	      continue
40	   continue
	   ffu(i)= fw1 * fi(2) + fw2 * fi(1)
	   ffd(i)= fw1 * fi(3) + fw2 * fi(4)
10	continue
	return
	end
 	 
