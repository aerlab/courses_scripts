
	subroutine rain_cal ( ib, dz, prwc, trn, wrn, wwrn )
! *********************************************************************
! trn, wrn, and wwrn are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
! due to the Mie scattering of rain for a given layer. 
!                        Jan. 19, 1993
! *********************************************************************
	implicit none
        include 'Rad4S.h'
        include 'cloud_rain.h'
!********** INPUTVAR *************************************
        integer:: ib             !  spectral band
        real::  prwc(nvx),  & ! rain content 
                 dz(nvx)     ! height of vertical layer
!*******  output var ***********************************
        real ::    trn(nvx), &   ! optical depth 
                   wrn(nvx), &   ! single scattering albedo 
                wwrn(nvx, 4)   ! expansion of phase function 

!*** *** tmp var
        real::  x1, x2, x3,  x4, y1, y2, y3, y4
        integer :: i, j

        x1 = grn(ib)
        x2 = x1 * grn(ib)
        x3 = x2 * grn(ib)
	x4 = x3 * grn(ib)
	y1 = 3.0 * x1
	y2 = 5.0 * x2
	y3 = 7.0 * x3
	y4 = 9.0 * x4
	do 10 i = 1, nvx
	   if ( prwc(i) .lt. 1.0e-5 ) then
             trn(i) = 0.0
             wrn(i) = 0.0
             wwrn(i,1) = 0.0
             wwrn(i,2) = 0.0
             wwrn(i,3) = 0.0
             wwrn(i,4) = 0.0
           else
	     trn(i) = dz(i) * prwc(i) * brn(ib) / rwc
	     wrn(i) = wrnf(ib)
	     wwrn(i,1) = y1
	     wwrn(i,2) = y2
	     wwrn(i,3) = y3
	     wwrn(i,4) = y4
           endif
10	continue
	return
	end subroutine rain_cal 


	subroutine graup_cal ( ib, dz, pgwc, tgr, wgr, wwgr )
! *********************************************************************
! tgr, wgr, and wwgr are the optical depth, single scattering albedo,
! and expansion coefficients of the phase function ( 1, 2, 3, and 4 )
! due to the Mie scattering of graupel for a given layer. 
!                        Jan. 19, 1993
! *********************************************************************
	implicit none
        include 'Rad4S.h'
        include 'cloud_rain.h'
!********** INPUTVAR *************************************
        integer:: ib             !  spectral band
        real, dimension(nvx) :: pgwc,  & ! rain content 
                                dz
!*******  output var ***********************************
        real ::    tgr(nvx), &   ! optical depth 
                   wgr(nvx), &   ! single scattering albedo 
                wwgr(nvx, 4)   ! expansion of phase function 

!*** *** tmp var
        real::  x1, x2, x3,  x4, y1, y2, y3, y4
        integer :: i, j

        x1 = gg(ib)
        x2 = x1 * gg(ib)
        x3 = x2 * gg(ib)
	x4 = x3 * gg(ib)
	y1 = 3.0 * x1
	y2 = 5.0 * x2
	y3 = 7.0 * x3
	y4 = 9.0 * x4
	do 10 i = 1, nvx
	   if ( pgwc(i) .lt. 1.0e-5 ) then
             tgr(i) = 0.0
             wgr(i) = 0.0
             wwgr(i,1) = 0.0
             wwgr(i,2) = 0.0
             wwgr(i,3) = 0.0
             wwgr(i,4) = 0.0
           else
	     tgr(i) = dz(i) * pgwc(i) * bg(ib) / gwc
             wgr(i) = wgf(ib)
	     wwgr(i,1) = y1
	     wwgr(i,2) = y2
	     wwgr(i,3) = y3
	     wwgr(i,4) = y4
           endif
10	continue
	return
	end subroutine graup_cal
