	subroutine rad( nv, nv1, ndfs, mdfs, &
	                pp, pt, ph, po, ataup, atype, &
                        cldwt, cldwtreff, cldice, cldicereff, rainwt, graupwt, & 
 			as, ee, pts, u0, ss,  &
	                umco2, umch4, umn2o, & 
		        fds, fus, fdir, fuir)

	
! *********************************************************************
! In this radiation scheme,  six  and  12 bands are selected for solar 
! and thermal IR regions, respectively. The spectral division is below: 
! 0.2 - 0.7 um, 0.7 - 1.3 um, 1.3 - 1.9 um, 1.9 - 2.5 um, 2.5 -3.5 um,
! 3.5 - 4.0 um, and 2200 - 1900 cm**-1, 1900 - 1700 cm**-1, 1700 -1400
! cm**-1,  1400 - 1250 cm**-1,  1250 - 1100 cm**-1, 1100 - 980 cm**-1,
! 980 - 800 cm**-1,  800 - 670 cm**-1,  670 - 540 cm**-1, 540 - 400 cm
! **-1,  400 - 280 cm**-1,  280 - 0 cm**-1,  where  the index  for the
! spectral band ( ib = 1, 2, ..., 18 ) is defined.
!
!                       **********************
!                       *  INPUT PARAMETERS  *
!                       **********************
!              as(mbs)   solar surface albedo, mbs = 6
!              u0        cosine of solar zenith angle
!              ss        solar constant ( W / m ** 2 )
!              pts       surface temperature ( K )
!              ee(mbir)  IR surface emissivity, mbir = 12
!              pp(nv1)   atmospheric pressure ( mb )
!              pt(nv1)   atmospheric temperature ( K )
!              ph(nv1)   water vapor mixing ratio ( kg / kg )
!              po(nv1)   ozone mixing ratio ( kg / kg )
!              pre(nv)   effective radius of water cloud ( um )
!              plwc(nv)  liquid water content ( g / m ** 3 )
!              pde(nv)   effective size of ice cloud ( um )
!              piwc(nv)  ice water content ( g / m ** 3 )
!              prwc(nv)  rain water content ( g / m ** 3 )
!              pgwc(nv)  graupel water content ( g / m ** 3 )
!
! Note:  (1)  as(mbs) and ee(mbir) consider the substantial wavelength
!             dependence of surface albedos and emissivities.
!        (2)  For CO2, CH4 and N2O, uniform mixing is assumed  through
!             the atmosphere with concentrations of 330, 1.6 and  0.28
!             ppmv, respectively.
!        (3)  nv, nv1, ndfs, mdfs, ndfs4, mb, mbs, mbir,  and  nc  are  
!             given through 'param.ipt'. 
!        (4)  nv1 and 1 are the surface and top levels, respectively.
!
!                       **********************
!                       *  OUTPUT PARAMETERS  *
!                       **********************
!              fds(nv1)   downward solar flux ( W / m ** 2 )
!              fus(nv1)   upward solar flux ( W / m **2 )
!              dts(nv)    solar heating rate ( K / day )
!              fdir(nv1)  downward IR flux ( W / m ** 2 )
!              fuir(nv1)  upward IR flux ( W / m **2 )
!              dtir(nv)   IR heating rate ( K / day )
!              fd(nv1)    downward net flux ( W / m ** 2 )
!              fu(nv1)    upward net flux ( W / m **2 )
!              dt(nv)     net heating rate ( K / day )
!
! Note:  Solar, IR, and net represent 0.2 - 0.4 um, 2200 - 0 cm**-1,
!        and  entire spectral regions, respectively.
!
! *********************************************************************
!	common /dfsout/ fu1(nv1x), fd1(nv1x)
!	common /planci/ bf(nv1x), bs
!	common /ic/ ti(nvx), wi(nvx), wwi(nvx,4)
!	common /wat/ tw(nvx), ww(nvx), www(nvx,4)
!	common /ray/ tr(nvx), wr(nvx), wwr(nvx,4)
!	common /aero/ ta(nvx), wa(nvx), wwa(nvx,4)
!	common  /aero_profile/ ataup(nvx), atau 
!	dimension as(mbs), ee(mbir)
!	integer gridx, gridy
 
        implicit none
	include 'Rad4S.h' 
	
! input
       real, dimension(nv1x):: pp, pt, ph, po
       real as(mbs), ee(mbir)
       real u0, ss                    ! sza and solar constant
       real umco2, umch4, umn2o
       integer nv, nv1, ndfs, mdfs
       real ataup(nvx), atype(nvx)
       real, dimension(nvx):: cldwt, cldwtreff, cldice, cldicereff, rainwt, graupwt  
       real pts
!output
       real, dimension(nv1x):: fds, fus, fdir, fuir

! internal variables

! kg(mb) is the number of intervals to perform the g-quadrature in
! each band to consider the nongray gaseous absorption.  In total,
! we need to perform 121 spectral calculations in  the  scattering
! problem for each atmospheric profile.
       real kg(mb) 
       real f0
       real dz(nvx)      ! thickness at each layer
       real trp(nvx)     !  p(mb)/t(k)*dz 
       integer mbn       ! starting band #. in night, mbn=7

! tau and SSA, and G for rayleight 
       real, dimension(nvx):: tr, wr    
       real wwr(nvx, 4)

! tau and SSA, and G for aerosol 
       real, dimension(nvx):: ta, wa    
       real wwa(nvx, 4)

! tau and SSA, and G for cld
       real, dimension(nvx):: tw, ww, ti, wi    
       real www(nvx, 4), wwi(nvx, 4)

! tau and SSA, and G for rain
       real, dimension(nvx):: trn, wrn, tgr, wgr    
       real wwrn(nvx, 4), wwgr(nvx, 4)

       
! tgcon, gas concentration
       real tgm(nvx)         

! plank function 
       real bf(nv1x), bs

! gases absoprtion
       real tg(nvx) 
       real hk               ! corresponding weight in each band

! after composite
       real, dimension(nvx):: tt, wc, wc1, wc2, wc3, wc4

! showrtwave and longwave flux in each band 
       real, dimension(nv1x):: fu1, fd1

! other
       real fuq1, fuq2
       integer i, ib, ig

! initizlize values

	data kg / 10, 8, 12, 7, 12, 5,& 
                 2, 3, 4, 4, 3, 5, 2, 10, 12, 7, 7, 8 /
	f0 = 1.0 / 3.14159

	do 10 i = 1, nv1
	   fds(i) = 0.0
	   fus(i) = 0.0
	   fdir(i) = 0.0
	   fuir(i) = 0.0
10	continue

!  calcualte thickness at each layer
	call thicks(pp, pt, nv, dz)
	print*, 'dz is ', dz


!  calcaulte p/(t*dz) for calculating raylae, 
!  14.6337 = R(287)/g(9.806)/2.
        do 90 i = 1, nv 
	  trp(i) = 14.6337 * ( pp(i) + pp(i+1) ) &
	           * alog( pp(i+1) / pp(i) ) 
90     continue

! 
! note 0.03 is the limit in surface albedo paramterization
! less than 0.03, the albedo values could be incorrect.
! see leaft 2 for details
	if ( u0 .le. 0.03 ) then
          mbn = mbs + 1
        else
          mbn = 1
        endif

! loop through each band
       
	do 20 ib = mbn, mb

          ! ray leight
	   call rayle ( trp, ib, u0, nv, tr, wr, wwr )

 400       format(1x, 3(2x, e10.4))	   
           
	  ! aerosol 
	   call aero_cal(ib,  nv, ataup, atype, ta, wa, wwa )
          
	  ! water  ! debug by Xi Chen 04/2022
	   call water_cal(ib,  dz, cldwtreff, cldwt, tw, ww, www )
	  
	  ! ice   ! debug by Xi Chen 04/2022
	   call ice_cal(ib,  dz, cldicereff, cldice, ti, wi, wwi )
         
	  ! rain  ! debug by Xi Chen 04/2022
	   call rain_cal(ib,  dz, rainwt, trn, wrn, wwrn )

	  ! graup  ! debug by Xi Chen 04/2022
	   call graup_cal(ib,  dz, graupwt, tgr, wgr, wwgr )

          ! gaes
	   call gascon ( ib, nv, nv1, pp, pt, ph, tgm )
          
	  ! for each band 
	  if ( ib .gt. mbs ) then
             call planck ( nv1, ib, pt, pts, bf, bs )
           endif
	   do 30 ig = 1, kg(ib)
	      call gases (nv, nv1,  ib, ig, pp, pt, ph, po, umco2, &
	                   umch4, umn2o, hk, tg )

              call comscp(nv, tr, wr, tgm, tg, ta, wa, &
                              ti, wi, tw, ww, trn, wrn, &
                              tgr, wgr, &
                              wwr, wwa, wwi, www, wwrn, wwgr, &
                              tt,  wc, wc1, wc2, wc3, wc4)
	       
	      if ( ib .le. mbs ) then
                call qfts ( ndfs, mdfs, ib, as(ib), u0, f0, &
		            wc1, wc2, wc3, wc4, wc, tt, &
			    fu1, fd1)
	        do 40 i = 1, nv1
                   fds(i) = fds(i) + fd1(i) * hk
                   fus(i) = fus(i) + fu1(i) * hk
40              continue
              else
                call qfti (ndfs, mdfs, ib, ee(ib-mbs), &
		            bs, bf, wc1, wc2, wc3, wc4, wc, tt, &
			    fu1, fd1)
	        do 50 i = 1, nv1
                   fdir(i) = fdir(i) + fd1(i) * hk
                   fuir(i) = fuir(i) + fu1(i) * hk
50              continue
              endif
30	   continue  
20	continue
	fuq1 = ss / 1340.0
! In this model, we used the solar spectral irradiance determined by
! Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained 
! in the spectral region 0.2 - 4.0 um.
	fuq2 = bs * 0.03 * 3.14159 * ee(12)
! fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
! hk of 0.03.
	do 60 i = 1, nv1
           fds(i) = fds(i) * fuq1
           fus(i) = fus(i) * fuq1
           fuir(i) = fuir(i) + fuq2
!	   fd(i) = fds(i) + fdir(i)
!	   fu(i) = fus(i) + fuir(i)
60	continue

!	do 70 i = 1, nv
!	   xx = fds(i) -fus(i) - fds(i+1) + fus(i+1)
!	   dts(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
!	   xx = fdir(i) -fuir(i) - fdir(i+1) + fuir(i+1)
!	   dtir(i) = 8.4392 * xx / ( pp(i+1) - pp(i) )
!	   dt(i) = dts(i) + dtir(i)
!70	continue
	return
	end
