
	subroutine water_cal ( ib, dz, pre, plwc, tw, ww, www )
! *********************************************************************
! tw, ww, and www are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and
! 4) due to the Mie scattering of water clouds for a given layer. 
! By using the mean single scattering properties of the eight drop
! size distributions in each spectral band, the single scattering
! properties of a water cloud with the given liquid water content
! and effective radius are obtained by interpolating (Eqs. 4.25 -
! 4.27 of Fu, 1991). 
! *********************************************************************
       implicit none 
       include 'Rad4S.h'
       include 'cloud_rain.h'

!********** INPUTVAR *************************************
        integer:: ib             !  spectral band
                  
        real, dimension(nvx) :: pre, plwc,  & ! water effective, ice water content
                                dz           ! height of vertical layer

!*******  output var ***********************************
        real ::    tw(nvx), &    ! optical depth 
                   ww(nvx), &    ! single scattering albedo 
                   www(nvx, 4)   ! expansion of phase function 

!*** *** tmp var
        real::  x1, x2, ggtmp, x3, x4
	integer :: i, j 
	
	do 10 i = 1, nvx
	   if ( plwc(i) .lt. 1.0e-5 ) then
             tw(i) = 0.0
             ww(i) = 0.0
             www(i,1) = 0.0
             www(i,2) = 0.0
             www(i,3) = 0.0
             www(i,4) = 0.0
           else
	     if ( pre(i) .lt. re(1) ) then
! A cloud with the effective radius smaller than 4.18 um is assumed
! to have an effective radius of 4.18 um with respect to the single
! scattering properties.  
	       tw(i) = dz(i) * plwc(i) * bz(1,ib) / fl(1)
 	       ww(i) = wz(1,ib)
               x1 = gz(1,ib)
               x2 = x1 * gz(1,ib)
               x3 = x2 * gz(1,ib)
	       x4 = x3 * gz(1,ib)
	       www(i,1) = 3.0 * x1
	       www(i,2) = 5.0 * x2
	       www(i,3) = 7.0 * x3
	       www(i,4) = 9.0 * x4
	     elseif ( pre(i) .gt. re(nc) ) then
! A cloud with the effective radius larger than 31.23 um is assumed
! to have an effective radius of 31.18 um with respect to the single
! scattering properties.  
	       tw(i) = dz(i) * plwc(i) * bz(nc,ib) / fl(nc)
	       ww(i) = wz(nc,ib)
	       x1 = gz(nc,ib)
               x2 = x1 * gz(nc,ib)
	       x3 = x2 * gz(nc,ib)
               x4 = x3 * gz(nc,ib)
	       www(i,1) = 3.0 * x1
	       www(i,2) = 5.0 * x2
	       www(i,3) = 7.0 * x3
	       www(i,4) = 9.0 * x4
	     else
	       j = 1
	    
	    do while ((pre(i) .lt. re(j)) .or. & 
	                (pre(i).gt. re(j+1))) 
	     j = j + 1
            enddo
	       tw(i) = dz(i) * plwc(i) * ( bz(j,ib) / fl(j) + &
     	       ( bz(j+1,ib) / fl(j+1) - bz(j,ib) / fl(j) ) /  &
     	       ( 1.0 / re(j+1) - 1.0 / re(j) ) * ( 1.0 / pre(i) &
     	       - 1.0 / re(j) ) ) 
	       ww(i) = wz(j,ib) + ( wz(j+1,ib) - wz(j,ib) ) / &
     	       ( re(j+1) - re(j) ) * ( pre(i) - re(j) ) 
	       ggtmp = gz(j,ib) + ( gz(j+1,ib) - gz(j,ib) ) / &
                    ( re(j+1) - re(j) ) * ( pre(i) - re(j) ) 
               x1 = ggtmp
                x2 = x1 * ggtmp
               x3 = x2 * ggtmp
               x4 = x3 * ggtmp
	       www(i,1) = 3.0 * x1
	       www(i,2) = 5.0 * x2
	       www(i,3) = 7.0 * x3
	       www(i,4) = 9.0 * x4
	     endif
           endif
10	continue
	return
	end subroutine water_cal


        subroutine ice_cal ( ib, dz, pde, piwc, ti, wi, wwi )
! *********************************************************************
! ti, wi, and wwi are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and
! 4) due to the scattering of ice clouds for a given layer.
! *********************************************************************
       implicit none
       include 'Rad4S.h'
       include 'cloud_rain.h'

!********** INPUTVAR *************************************
        integer:: ib             !  spectral band
                  
        real, dimension(nvx) :: pde, piwc,  & ! water effective, ice water content
                                dz           ! height of vertical layer

!*******  output var ***********************************
        real ::    ti(nvx), &    ! optical depth 
                   wi(nvx), &    ! single scattering albedo 
                   wwi(nvx, 4)   ! expansion of phase function 

!***** tmp var ************************
        integer :: i, ibr,j 
        real :: fw1, fw2, fw3, fd, wf1, wf2, wf3, wf4
        real :: x1, x2, ggtmp, x3, x4
        do 10 i = 1, nvx
           if ( piwc(i) .lt. 1.0e-5 ) then
             ti(i) = 0.0
             wi(i) = 0.0
             wwi(i,1) = 0.0
             wwi(i,2) = 0.0
             wwi(i,3) = 0.0
             wwi(i,4) = 0.0
           else
! The constant 1000.0 below is to consider the units of dz(i) is km.
             fw1 = pde(i)
             fw2 = fw1 * pde(i)
             fw3 = fw2 * pde(i)
             ti(i) = dz(i) * 1000.0 * piwc(i) * ( ap(1,ib) + &
             ap(2,ib) / fw1 + ap(3,ib) / fw2 )
             wi(i) = 1.0 - ( bp(1,ib) + bp(2,ib) * fw1 + &
             bp(3,ib) * fw2 + bp(4,ib) * fw3 )
             if ( ib .le. mbs ) then
               fd = dps(1,ib) + dps(2,ib) * fw1 + &
              dps(3,ib) * fw2 + dps(4,ib) * fw3
               wf1 = cps(1,1,ib) + cps(2,1,ib) * fw1 + &
              cps(3,1,ib) * fw2 + cps(4,1,ib) * fw3
               wwi(i,1) = ( 1.0 - fd ) * wf1 + 3.0 * fd
               wf2 = cps(1,2,ib) + cps(2,2,ib) * fw1 + &
              cps(3,2,ib) * fw2 + cps(4,2,ib) * fw3
               wwi(i,2) = ( 1.0 - fd ) * wf2 + 5.0 * fd
               wf3 = cps(1,3,ib) + cps(2,3,ib) * fw1 + &
              cps(3,3,ib) * fw2 + cps(4,3,ib) * fw3
               wwi(i,3) = ( 1.0 - fd ) * wf3 + 7.0 * fd
               wf4 = cps(1,4,ib) + cps(2,4,ib) * fw1 + &
              cps(3,4,ib) * fw2 + cps(4,4,ib) * fw3
               wwi(i,4) = ( 1.0 - fd ) * wf4 + 9.0 * fd
             else
               ibr = ib - mbs
               ggtmp = cpir(1,ibr) + cpir(2,ibr) * fw1 + &
               cpir(3,ibr) * fw2 + cpir(4,ibr) * fw3
               x1 = ggtmp
               x2 = x1 * ggtmp
               x3 = x2 * ggtmp
               x4 = x3 * ggtmp
               wwi(i,1) = 3.0 * x1
               wwi(i,2) = 5.0 * x2
               wwi(i,3) = 7.0 * x3
               wwi(i,4) = 9.0 * x4
             endif
           endif
10      continue
        return
        end subroutine ice_cal



