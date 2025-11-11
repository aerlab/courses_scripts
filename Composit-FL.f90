	subroutine comscp(nv, tr, wr, tgm, tg, ta, wa, & 
                              ti, wi, tw, ww, trn, wrn, & 
                              tgr, wgr, &
                              wwr, wwa, wwi, www, wwrn, wwgr, & 
                              tt,  wc, wc1, wc2, wc3, wc4)
	                  
			  
! *********************************************************************
! This subroutine is used to  COMbine Single-Scattering Properties  due
! to  ice crystals,  water droplets, and  Rayleigh molecules along with
! H2O continuum absorption and nongray gaseous absorption.  See Section
! 3.4 of Fu (1991). wc, wc1, wc2, wc3, and wc4, are total (or combined)
! single - scattering  albedo,  and   expansion   coefficients  of  the
! phase function ( 1, 2, 3, and 4 ) in nv layers. tt(nv) are the normal
! optical depth ( from the top of the atmosphere to a given level ) for
! level 2 - level nv1( surface ). The single-scattering  properties  of
! rain and graupel are also incorporated in ( Jan. 19, 1993 ).
! *********************************************************************
!	common /dfsin/ wc1(nvx), wc2(nvx), wc3(nvx), wc4(nvx),& 
!     	               wc(nvx), tt(nvx)
	implicit none
	include 'Rad4S.h'

!input
        real, dimension(nvx):: tr, wr,  &   ! Rayleigh opt, ssa 
	                       tgm, tg, &   ! gas
                               ta, wa , &   ! aerosol
                               ti, wi,  &   !ice 
                               trn, wrn, &   ! rain
                               tw, ww  , &   ! water
                               tgr, wgr

        real, dimension(nvx, 4):: wwr, wwa, www, wwi, wwrn, wwgr ! phase	 
	integer nv
	
!output
        real, dimension(nvx):: tt, wc, wc1, wc2, wc3, wc4

!internal
        real tc(nvx)	
	integer i
	real tas, fw, tis, tws, trns, tgrs 

	do 10 i = 1, nv
	   tc(i) =  tr(i) + tgm(i) + tg(i) + ta(i) + & 
                    ti(i) + tw(i) + trn(i) + tgr(i)

           tis = ti(i) * wi(i)
           tws = tw(i) * ww(i)
           trns = trn(i) * wrn(i)
           tgrs = tgr(i) * wgr(i) 
	   tas  = ta(i) * wa(i)                      !scattering
           fw = tr(i) +  tas + tis + tws  + trns + tgrs !total scattering
	   wc(i) =  fw / tc(i)                       ! SSA
           if ( fw .lt. 1.0e-20 ) then
             wc1(i) = 0.0
             wc2(i) = 0.0
             wc3(i) = 0.0
             wc4(i) = 0.0
	   else
             wc1(i) = ( tr(i) * wwr(i,1) + tas * wwa(i,1) + & 
                        tis * wwi(i,1) + tws * www(i,1)   + & 
                        trns * wwrn(i,1) + tgrs*wwgr(i,1) )/fw
             wc2(i) = ( tr(i) * wwr(i,2) + tas * wwa(i,2) + & 
                        tis * wwi(i,2) + tws * www(i,2)   + & 
                        trns * wwrn(i,2) +  tgrs*wwgr(i,2) )/fw
             wc3(i) = ( tr(i) * wwr(i,3) + tas * wwa(i,3) + &  
                        tis * wwi(i,3) + tws * www(i,3)   + & 
                        trns * wwrn(i,3) +  tgrs*wwgr(i,3)  )/fw
             wc4(i) = ( tr(i) * wwr(i,4) + tas * wwa(i,4) + & 
                        tis * wwi(i,4) + tws * www(i,4)  +  & 
                        trns * wwrn(i,4) + tgrs*wwgr(i,4)  )/fw
	   endif
10	continue
	tt(1) = tc(1)
	do 20 i = 2, nv
	   tt(i) = tt(i-1) + tc(i)
20	continue
	return
	end
