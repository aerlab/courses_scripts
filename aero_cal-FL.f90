
!
! purpose: to include aerosol properties into the FuLiou code.
!          to let the algorithms simple and efficient, so
!          it can be used in the RAMS model

! Input:  FULiou band index, AOT at 0.55um., aerosol type
! Output:  tau profiles, single scattering albedo, and asymeteric 
!          factor in each layers for the specific band


!=======================================================================
!   main routines
!=======================================================================

	subroutine aero_cal (ib, nv, ataup, atype, ta, wa, wwa)
	implicit none
	include 'Rad4S.h'
	include 'Aerosol.h'
	
	
! input
	integer ib,  nv
        real ataup(nvx), atype(nvx)

	
! output	
	real ta(nvx), wa(nvx), wwa(nvx, 4)

	
! other tmp vars
         integer i, j, itype 
	 real tmptau           ! total tau, tau @ 0.5um
	 	 
! first get opt properties in the specifid band
! 4 values for phase function expansion assuming H-G phase function.
! total tau in this band
	 do i = 1, nv
	   itype = atype(i)
	   tmptau = ataup(i) *  a_ext(ib, itype)/a_ext(1, itype) 
	   if ( tmptau .lt. 1.0e-6) then
	      ta(i) = 0.0
	      wa(i) = 0
	      wwa(i, 1) = 0.0
	      wwa(i, 2) = 0.0
	      wwa(i, 3) = 0.0
	      wwa(i, 4) = 0.0
	   else
	      ta(i) = tmptau  
	      wa(i) = a_ssa(ib, itype)
	      wwa(i, 1) = 3.* a_asy(ib, itype) 
	      wwa(i, 2) = 5 * a_asy(ib, itype)**2 
	      wwa(i, 3) = 7 * a_asy(ib, itype)**3
	      wwa(i, 4) = 9 * a_asy(ib, itype)**4
	   !   print*, 'ta', ta(i), wa(i)
	   endif
	 enddo       
        
	END

