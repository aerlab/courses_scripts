!
! 	subroutine rayle ( ib, u0 )
! *********************************************************************
! tr, wr, and wwr are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and
! 4 ) due to the Rayleigh scattering for a given layer.
! *********************************************************************
	
	subroutine rayle ( trp, ib, u0, nv, tr, wr, wwr )
        
	implicit none
	include 'Rad4S.h' 

! input
        integer ib, nv
	real u0
	real trp(nvx)
	

! internal	
	real ri(mbs)
	real x
	integer i

! output, tau, ssa and g for rayleigh scattering
        real, dimension(nvx):: tr, wr
	real wwr(nvx, 4)	

! initial
        data ri / 0.9022e-5, 0.5282e-6, 0.5722e-7, &
     	          0.1433e-7, 0.4526e-8, 0.1529e-8 /
	
! starts, ray vlues are zero for in infrared bands
	if ( ib .le. mbs ) then
	  if ( ib .eq. 1 ) then
	    x = -3.902860e-6*u0*u0+6.120070e-6*u0+4.177440e-6
	  else
	    x = ri(ib)
	  endif
	  do 100 i = 1, nv
	     tr(i) = trp(i) * x
	     wr(i) = 1.0
	     wwr(i,1) = 0.0
	     wwr(i,2) = 0.5
	     wwr(i,3) = 0.0
	     wwr(i,4) = 0.0
100       continue
	else
	  do 200 i = 1, nv
	     tr(i) = 0.0
	     wr(i) = 0.0
	     wwr(i,1) = 0.0
	     wwr(i,2) = 0.0
	     wwr(i,3) = 0.0
	     wwr(i,4) = 0.0
200       continue
	endif
	return
	end
