	subroutine gascon ( ib, nv, nv1, pp, pt, ph, tgm )
! *********************************************************************
! tgm(nv) are the optical depthes due to water vapor continuum absorp-
! tion in nv layers for a given band ib. We include continuum absorp-
! tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
! wavenumbers of each band in this region. 
! *********************************************************************
!	
        implicit none
	include 'Rad4S.h'

! input
        integer ib, nv, nv1
	real, dimension(nv1x):: pp, pt, ph

! output
        real tgm(nvx)
	
! internal
        integer i
        real vv(mb)	
 	data vv / 10*0.0, 1175.0, 1040.0, 890.0, 735.0, &
     	          605.0, 470.0, 340.0, 0.0 /
	
	
	if ( ib .gt. 10 .and. ib .lt. mb ) then
	   call qopcon (nv, nv1, vv(ib), pp, pt, ph, tgm )
	else
	   do 10 i = 1, nv
              tgm(i) = 0.0
10	   continue
	endif
	return
	end


      subroutine qopcon ( nv, nv1, vv, pp, pt, ph, tg )
      implicit none
      include 'Rad4S.h' 
 
! input
      real vv
      integer nv1, nv
      real, dimension(nv1x):: pp, pt, ph
      

! output 
     real tg(nvx)       
     
    
! internal      
      integer i
      real ff(nv1x),pe(nv1x)
      real x,y,z,r,s,w

! initial 
      x = 4.18
      y = 5577.8
      z = 0.00787
      r = 0.002
      s = ( x + y * exp ( - z * vv ) ) / 1013.25

      do 3 i = 1, nv1
       pe(i) = pp(i) * ph(i) / ( 0.622 + 0.378 * ph(i) )
       w = exp ( 1800.0 / pt(i) - 6.08108 )
       ff(i) = s * ( pe(i) + r * pp(i) ) * w
3      continue

      do 5 i = 1, nv
       tg(i) = ( ff(i) * ph(i) + ff(i+1) * ph(i+1) )*   &
               ( pp(i+1) - pp(i) ) * 0.5098835
5      continue

      return
      end




      subroutine gases (nv, nv1, ib, ig, pp, pt, ph, po, umco2, umch4, umn2o, hk, tg )
!  *****************************************************************
!  tg(nv) are the optical depthes due to nongray gaseous absorption, 
!  in nv layers for a given band ib and cumulative probability ig. 
!  *****************************************************************
      implicit none
      include 'Rad4S.h'
      include 'CKD.h'


! input 
      real umco2, umch4, umn2o
      real, dimension(nv1x):: pp, pt, ph, po        
      integer ib, ig 
      integer nv, nv1

! output
      real hk
      

! internal
      real fkg(nv1x), fkga(nv1x), fkgb(nv1x), pq(nv1x)
      real tg1(nvx), tg2(nvx), tg3(nvx)
      real, dimension(nvx):: tg
      real fk


       ZERO:   select case  (ib)

      
       case (1) 	 
!  ---------------------------------------------------------------------
!  In this band ( 50000 - 14500 cm**-1 ), we have considered the nongray
!  gaseous absorption of O3.    619.618 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
         fk = fk1o3(ig)
         call qopo3s (nv, fk, pp, po, tg )
         hk = 619.618 * hk1(ig)
      
       case (2)
!  ---------------------------------------------------------------------
!  In this band ( 14500 - 7700 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.  484.295 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
!        call qks ( c2h2o(1,1,ig), fkg )
	call  qks (nv1, pp, pt, c2h2o(1,1,ig) , fkg )
        call qoph2o (nv, fkg,pp, ph,  tg )
        hk = 484.295 * hk2(ig)
      
      case (3)  
!  ---------------------------------------------------------------------
!  In this band ( 7700 - 5250 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 149.845 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
!      call qks ( c3h2o(1,1,ig), fkg )
      call  qks (nv1, pp, pt, c3h2o(1,1,ig) , fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = 149.845 * hk3(ig)
	
       case (4)
!  ---------------------------------------------------------------------
!  In this band ( 5250 - 4000 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 48.7302 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
!      call qks ( c4h2o(1,1,ig), fkg )
      call  qks (nv1, pp, pt, c4h2o(1,1,ig) , fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = 48.7302 * hk4(ig)


      case (5)
!  ---------------------------------------------------------------------
!  In this band ( 4000 - 2850 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 31.6576 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
!      call qks ( c5h2o(1,1,ig), fkg )
      call  qks (nv1, pp, pt, c5h2o(1,1,ig) , fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = 31.6576 * hk5(ig)


     case (6)  
!  ---------------------------------------------------------------------
!  In this band ( 2850 - 2500 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O. 5.79927 is the solar energy contained in
!  the band in units of Wm**-2.
!  ---------------------------------------------------------------------
!      call qks ( c6h2o(1,1,ig), fkg )
      call  qks (nv1, pp, pt, c6h2o(1,1,ig) , fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = 5.79927 * hk6(ig)



      case (7)
!  ---------------------------------------------------------------------
!  In this band ( 2200 - 1900 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
      call  qki (nv1, pp, pt,  c7h2o(1,1,ig), fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = hk7(ig)

     case (8)
!  ---------------------------------------------------------------------
!  In this band ( 1900 - 1700 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
!      call qki ( c8h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c8h2o(1,1,ig), fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = hk8(ig)

     case (9) 
!  ---------------------------------------------------------------------
!  In this band ( 1700 - 1400 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
!WJ 9    call qki ( c9h2o(1,1,ig), fkg )
!      call qki ( c9h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c9h2o(1,1,ig), fkg )
      call qoph2o (nv, fkg,pp, ph,  tg )
      hk = hk9(ig)
!WJ      goto 20

     case(10)
!  ---------------------------------------------------------------------
!  In this band ( 1400 - 1250 cm**-1 ), we have considered the 
!  overlapping absorption of H2O, CH4, and N2O by approach one of 
!  Fu(1991).
!  ---------------------------------------------------------------------
!WJ 10   call qki ( c10h2o(1,1,ig), fkg )
!      call qki ( c10h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c10h2o(1,1,ig), fkg )
      call qoph2o (nv, fkg,pp, ph,  tg1 )
!      call qki ( c10ch4, fkg )
      call  qki (nv1, pp, pt,  c10ch4, fkg )
      call qopch4 (nv,  fkg, pp, tg2 )
!      call qki ( c10n2o, fkg )
      call  qki (nv1, pp, pt,  c10n2o, fkg )
!      call qopn2o ( fkg, tg3 )
      call  qopn2o (nv, pp, fkg, tg3 )
      do i = 1, nv
       tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
       end do
      hk = hk10(ig)
!WJ      goto 20


      case (11)
!  ---------------------------------------------------------------------
!  In this band ( 1250 - 1100 cm**-1 ), we have considered the 
!  overlapping absorption of H2O, CH4, and N2O by approach one of 
!  Fu(1991).
!  ---------------------------------------------------------------------

!      call qki ( c11h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c11h2o(1,1,ig), fkg )
      call qoph2o (nv, fkg,pp, ph,  tg1 )

!      call qki ( c11ch4, fkg )
      call  qki (nv1, pp, pt,  c11ch4, fkg )
!      call qopch4 ( fkg, tg2 )
      call qopch4 (nv,  fkg, pp,  tg2 )
!      call qki ( c11n2o, fkg )
      call  qki (nv1, pp, pt,  c11n2o, fkg )
!      call eopn2o ( fkg, tg3 )
      call  qopn2o (nv, pp, fkg, tg3 )
      do i = 1, nv
       tg(i) = tg1(i) + tg2(i)/1.6*umch4 + tg3(i)/0.28*umn2o
       end do

       hk = hk11(ig)
	

        case (12)
!  ---------------------------------------------------------------------
!  In this band ( 1100 - 980 cm**-1 ), we have considered the overlapping
!  absorption of H2O and O3 by approach one of Fu(1991).
!  ---------------------------------------------------------------------
!        call qkio3 ( c12o3(1,1,ig), fkg )
	call  qkio3 (nv1, pp, pt,  c12o3(1,1,ig), fkg )
!          call qopo3i ( fkg, tg1 )
       call  qopo3i (nv, pp, po,  fkg, tg1 )
!          call qki ( c12h2o, fkg )
      call  qki (nv1, pp, pt,  c12h2o, fkg )
          call qoph2o (nv, fkg,pp, ph,  tg2 )

          do i = 1, nv
           tg(i) = tg1(i) + tg2(i)
          end do

           hk = hk12(ig)

	

      case (13)
!  ---------------------------------------------------------------------
!  In this band ( 980 - 800 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
!         call qki ( c13h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c13h2o(1,1,ig), fkg )
!         call qoph2o ( fkg, tg )
        call qoph2o (nv, fkg,pp, ph,  tg )
         hk = hk13(ig)
	

       case (14)
!  ---------------------------------------------------------------------
!  In this band ( 800 - 670 cm**-1), we have considered the overlapping
!  absorption of H2O and CO2 by approach two of Fu(1991).
!  ---------------------------------------------------------------------
!WJ 14   do i = 1, nv1
      do i = 1, nv1
       if ( pp(i) .ge. 63.1 ) then
         pq(i) = ph(i)
       else
         pq(i) = 0.0
       endif
       end do
!      call qki ( c14hca(1,1,ig), fkga )
!      call qki ( c14hcb(1,1,ig), fkgb )
      call  qki (nv1, pp, pt,  c14hca(1,1,ig), fkga )
      call  qki (nv1, pp, pt,  c14hcb(1,1,ig), fkgb )
      do i = 1, nv1
       fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
       end do
!      call qophc ( fkg, tg)
      call  qophc (nv, pp,  fkg, tg )
      hk = hk14(ig)



       case  (15)  
!  ---------------------------------------------------------------------
!  In this band ( 670 - 540 cm**-1), we have considered the overlapping
!  absorption of H2O and CO2 by approach two of Fu(1991).
!  ---------------------------------------------------------------------
!WJ 15   do i = 1, nv1
      do i = 1, nv1
       if ( pp(i) .ge. 63.1 ) then
         pq(i) = ph(i)
       else
         pq(i) = 0.0
       endif
       end do
!      call qki ( c15hca(1,1,ig), fkga )
!      call qki ( c15hcb(1,1,ig), fkgb )
      call  qki (nv1, pp, pt,  c15hca(1,1,ig), fkga )
      call  qki (nv1, pp, pt,  c15hcb(1,1,ig), fkgb )
      do i = 1, nv1
       fkg(i) = fkga(i)/330.0*umco2 + pq(i) * fkgb(i)
       end do
!      call qophc ( fkg, tg)
      call  qophc (nv, pp,  fkg, tg )
      hk = hk15(ig)
!WJ      goto 20

      case (16)
!  ---------------------------------------------------------------------
!  In this band ( 540 - 400 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
!WJ 16   call qki ( c16h2o(1,1,ig), fkg )
!      call qki ( c16h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c16h2o(1,1,ig), fkg )
!      call qoph2o ( fkg, tg )
        call qoph2o (nv, fkg,pp, ph,  tg )
      hk = hk16(ig)
!WJ      goto 20

      case (17)
!  ---------------------------------------------------------------------
!  In this band ( 400 - 280 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
!WJ 17   call qki ( c17h2o(1,1,ig), fkg )
!      call qki ( c17h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c17h2o(1,1,ig), fkg )
!      call qoph2o ( fkg, tg )
        call qoph2o (nv, fkg,pp, ph,  tg )
      hk = hk17(ig)
!WJ      goto 20


      case (18)
!  ---------------------------------------------------------------------
!  In this band ( 280 - 000 cm**-1 ), we have considered the nongray
!  gaseous absorption of H2O.
!  ---------------------------------------------------------------------
!WJ 18   call qki ( c18h2o(1,1,ig), fkg )
!      call qki ( c18h2o(1,1,ig), fkg )
      call  qki (nv1, pp, pt,  c18h2o(1,1,ig), fkg )
!      call qoph2o ( fkg, tg )
        call qoph2o (nv, fkg,pp, ph,  tg )
      hk = hk18(ig)

      case default
         print*, 'ib = ', ib, 'IB choice is wrong'
         stop 'in gases'

!WJ 20   continue
      end select ZERO
        
        
       return
      end


! 338kfix
	subroutine qks (nv1, pp, pt, coefks, fkg )
! *********************************************************************
! fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
! **-1 for a given cumulative probability in nv1 layers. coefks(3,11)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 11 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at conditions other than those eleven
! pressures is interpolated linearly with pressure (Fu, 1991).
! *********************************************************************
      implicit none
      include 'Rad4S.h'

! input
      real, dimension(nv1x):: pp, pt 
      real coefks(3,11)
      integer nv1

! output
     real fkg(nv1x)

! internal     
     real  stanp(11)
     integer i1, i
     real x1, y1, x2, tk(nv1x) 
     data stanp / 10.0, 15.8, 25.1, 39.8, 63.1, 100.0,   &
      	             158.0, 251.0, 398.0, 631.0, 1000.0 /

       do i = 1, nv1
         if (  pt(i) .le. 180.0 ) then
	   tk(i) = 180.
	 else if( pt(i) .ge. 320 ) then
	   tk(i) = 320.
	 else
	   tk(i) = pt(i)
	 endif  
	end do
        
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefks(1,1) + coefks(2,1) * ( pt(i) - 245.0 )   &
             + coefks(3,1) * ( pt(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	     
	   elseif ( pp(i) .ge. stanp(11) ) then
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefks(1,10) + coefks(2,10) * ( pt(i) - 245.0 )   &
      	     + coefks(3,10) * y1 )
    	     x2 = exp ( coefks(1,11) + coefks(2,11) * ( pt(i) - 245.0 )   &
      	     + coefks(3,11) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(11) - stanp(10) )   &
      	     * ( pp(i) - stanp(10) )
	   else
!30	     continue
!	     if ( pp(i) .ge. stanp(i1) ) goto 20

             do 20  while ( pp(i) .ge. stanp(i1) )
	      i1 = i1+1
20 	     continue 
             
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefks(1,i1-1) + coefks(2,i1-1) * (pt(i)-245.0)   &
      	     + coefks(3,i1-1) * y1 )
	     x2 = exp ( coefks(1,i1) + coefks(2,i1) * ( pt(i) - 245.0 )   &
      	     + coefks(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )   &
      	     * ( pp(i) - stanp(i1-1) )
!	     goto 5
!20           i1 = i1 + 1
!	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qki (nv1, pp, pt,  coefki, fkg )
! *********************************************************************
! fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
! **-1 for a given cumulative probability in nv1 layers. coefki(3,19)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!         ln k = a + b * ( t - 245 ) + c * ( t - 245 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).
! *********************************************************************
	implicit none
	include 'Rad4S.h'

! input
     real, dimension(nv1x):: pp, pt 
     integer nv1 
     real  coefki(3,19)

! output
    real fkg(nv1x)
    
! internal
	real  stanp(19)
	integer i1, i
	real x1, y1, x2
	data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51,    &
      	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,   &
      	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( pt(i) - 245.0 )   &
             + coefki(3,1) * ( pt(i) - 245.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( pt(i) - 245.0 )   &
      	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( pt(i) - 245.0 )   &
      	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )   &
      	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 245.0 ) * ( pt(i) - 245.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (pt(i)-245.0)   &
      	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( pt(i) - 245.0 )   &
      	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )   &
      	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

	subroutine qkio3 (nv1, pp, pt,  coefki, fkg )
! *********************************************************************
! fkg(nv1) are the gaseous absorption coefficients in units of (cm-atm)
! **-1 for a given cumulative probability in nv1 layers. coefki(3,19)
! are the coefficients to calculate the absorption coefficient at the
! temperature t for the 19 pressures by
!         ln k = a + b * ( t - 250 ) + c * ( t - 250 ) ** 2
! and the absorption coefficient at  conditions  other  than  those 19
! pressures is interpolated linearly with pressure (Fu, 1991).
! *********************************************************************
	implicit none
	include 'Rad4S.h'
! input
      real, dimension(nv1x):: pp, pt 
      real coefki(3,19)
      integer nv1

! output
     real fkg(nv1x)

! internal
     integer i1, i
     real x1, y1, x2
     real  stanp(19)
     data stanp / 0.251, 0.398, 0.631, 1.000, 1.58, 2.51,    &
      	             3.98, 6.31, 10.0, 15.8, 25.1, 39.8, 63.1,   &
      	             100.0, 158.0, 251.0, 398.0, 631.0, 1000.0 /
	
	i1 = 1
	do 5 i = 1, nv1
	   if ( pp(i) .lt. stanp(1) ) then
  	     x1 = exp ( coefki(1,1) + coefki(2,1) * ( pt(i) - 250.0 )   &
             + coefki(3,1) * ( pt(i) - 250.0 ) ** 2 )
	     fkg(i) = x1 * pp(i) / stanp(1)
	   elseif ( pp(i) .ge. stanp(19) ) then
	     y1 = ( pt(i) - 250.0 ) * ( pt(i) - 250.0 )
	     x1 = exp ( coefki(1,18) + coefki(2,18) * ( pt(i) - 250.0 )   &
      	     + coefki(3,18) * y1 )
    	     x2 = exp ( coefki(1,19) + coefki(2,19) * ( pt(i) - 250.0 )   &
      	     + coefki(3,19) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(19) - stanp(18) )   &
      	     * ( pp(i) - stanp(18) )
	   else
30	     continue
	     if ( pp(i) .ge. stanp(i1) ) goto 20
	     y1 = ( pt(i) - 250.0 ) * ( pt(i) - 250.0 )
	     x1 = exp ( coefki(1,i1-1) + coefki(2,i1-1) * (pt(i)-250.0)   &
      	     + coefki(3,i1-1) * y1 )
	     x2 = exp ( coefki(1,i1) + coefki(2,i1) * ( pt(i) - 250.0 )   &
      	     + coefki(3,i1) * y1 )
	     fkg(i) = x1 + ( x2 - x1 ) / ( stanp(i1) - stanp(i1-1) )   &
      	     * ( pp(i) - stanp(i1-1) )
	     goto 5
20           i1 = i1 + 1
	     goto 30
	   endif
5  	continue
	return
	end

! 338kfix



       subroutine qopo3s (nv, fk, pp, po, tg )

      implicit none
      include 'Rad4S.h' 

! input
      real fk
      real, dimension(nv1x):: pp, po 
      integer nv

! output
      real tg(nvx)

! internal
      real fq 
      integer i


      fq = 238.08 * fk
      do 10 i = 1, nv
       tg(i) = ( po(i) + po(i+1) ) * ( pp(i+1) - pp(i) ) * fq
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 476.16 * fk
!20      continue
! 476.16 = 2.24e4 / M * 10.0 / 9.8, where M = 48 for O3.

      return
      end

      subroutine qoph2o (nv, fkg,pp, ph,  tg )
      implicit none
      include 'Rad4S.h'

!input
      real fkg(nv1x)
      real, dimension(nv1x):: ph, pp 
     integer nv 
     
!output
      real tg(nvx)
    
!internal    
     integer i

      do 10 i = 1, nv
       tg(i) = ( fkg(i) * ph(i) + fkg(i+1) * ph(i+1) )   &
               * ( pp(i+1) - pp(i) ) * 634.9205
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 1269.841
!20      continue
! 1269.841 = 2.24e4 / M * 10.0 / 9.8, where M = 18 for H2O.

      return
      end


      subroutine qopch4 (nv,  fkg, pp,  tg )
      implicit none
      include 'Rad4S.h'
 
! input
      real, dimension(nv1x):: pp, fkg
      integer nv

! outout
      real tg(nv1x)

! internal
      integer i

      do 10 i = 1, nv
       tg(i) = ( fkg(i)+fkg(i+1) ) * ( pp(i+1)-pp(i) ) * 6.3119e-4
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 1.26238e-3
!20      continue
! 1.26238e-3 = 2.24e4 / M * 10.0 / 9.8 * 1.6e-6 * M / 28.97, where 
! M = 16 for CH4.

      return
      end

      subroutine qopn2o (nv, pp, fkg, tg )
      implicit none
      include 'Rad4S.h'

!input
      real, dimension(nv1x):: fkg, pp
      integer nv

!output
      real tg(nvx)

!internal
      integer i

      do 10 i = 1, nv
       tg(i) = ( fkg(i)+fkg(i+1) ) * ( pp(i+1)-pp(i) ) * 1.10459e-4
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 2.20918e-4
!20      continue
! 2.20918e-4 = 2.24e4 / M * 10.0 / 9.8 * 0.28e-6 * M / 28.97, where
! M = 44 for N2O.

      return
      end

      subroutine qopo3i (nv, pp, po, fkg, tg )
      implicit none
      include 'Rad4S.h'

! input
     integer nv
     real, dimension(nv1x):: pp, fkg, po

! output
     real tg(nvx)
      integer i

      do 10 i = 1, nv
       tg(i) = ( fkg(i) * po(i) + fkg(i+1) * po(i+1) )   &
               * ( pp(i+1) - pp(i) ) * 238.08
10     continue

!      do 20 i = 1, nv
!         tg(i) = tg(i) * 476.16
!20      continue

      return
      end

      subroutine qophc (nv, pp,  fkg, tg )
      implicit none
      include 'Rad4S.h'

!input
      real, dimension(nv1x):: fkg, pp
      integer nv

!output
      real tg(nvx)

!internal 
      integer i
      do 10 i = 1, nv
       tg(i) = ( fkg(i) + fkg(i+1) ) * ( pp(i+1) - pp(i) ) * 0.5
10     continue
!  ------------------------
!  See page 86 of Fu (1991).
!  ------------------------

      return
      end
