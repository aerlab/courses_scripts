	implicit none
	integer, parameter:: ntau = 2, nratio = 11 
	include 'Rad4S.h'
	
	character*20 line
	character*41 a1,a2,a3,a4,a5,a6
	integer nlev(50), nv, nv1, ndfs, mdfs, ndfs4
! pressure, T, Water mixing   ratio, ozone, 
! downward, and upward sw and long flux                                       

  real, dimension(nv1x):: pp, pt, ph, po, fds, fus, fdir, fuir
  
! heating rate, sw, infrared, and total, thickness at each layer
  real, dimension(nvx) :: dts, dtir, dt, dz
      
! aerosol profile, aot vertical profiles, aerosol type
  real, dimension(nvx):: ataup, tmptaup, atype
  real, dimension(nvx):: cldwt, cldwtreff, cldice, cldicereff, rainwt, graupwt 

! other traces concentration
  real umco2, umch4, umn2o
  
! surface albedo
  real as(mbs), ee(mbir) 

! dust / smoke ratio
  real aotratio(nratio)


! others
  integer i, ii, jj
  real pts, u0, ss
	
	
! program testing
  real tau(ntau)
!  data tau / 0.0, 0.1, 0.3,  0.5, 0.6, 0.8, 1.2, 1.5, 1.8, 2.0/
  data tau /0.0,  1.0/
  data aotratio / 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0/

        nv = 30 
        nv1=nv+1
        ndfs = nv
        mdfs = nv1
        ndfs4= 4*ndfs
        
! set cloud = 0

    cldwt(1:nv1) = 0.0
    cldwtreff(1:nv1) = 0.0
    cldice(1:nv1) = 0.0 
    cldicereff(1:nv1) = 0.0
    rainwt(1:nv1) = 0.0
    graupwt(1:nv1)= 0.0 
    
! note, here maximum dimension of ataup = 40, I have let the AOT
! fraction in the
!bottom is equal to zero.

	print*, 'nvx = ', nvx, nv1
        open(8, file = "zprofile") 
	read ( 8, *) ( pp(i), pt(i), ph(i), po(i), ataup(i),& 
                      atype(i), i = 1, nv1 )
        close(8)

    
!     cldwt(1:nv1) = cldwt(1:nv1) * 2
!     do i = 1, nv1
!        if (cldwt(i) .gt. 1.0e-5)  cldwtreff(i) =  6.0
!     enddo   
!     
!     print*, 'cld wt:', cldwt
!     print*, 'cld wtreff:', cldwtreff
   

	
!** ss = solar flux
!** uo = cosine of solar zenith angle

        ss = 1365.0
        u0 = 0.5
	umCo2 = 330
	umch4 = 1.6
	umn2o = 0.28 
       
     
	pts = pt(nv1)
     
	as(1) = 0.1
	do 31 i=2,mbs
	  as(i) = as(1)
 31    continue

	ee(1) = 1.0
	do 32 i=mbs+1,mb
	  ee(i) = ee(1)
 32    continue
       
       
!       write(*,*) 'Atmospheric Profile' 
!       write(*, 33) (pp(i), pt(i), ph(i), po(i), ataup(i), i = 1, nv1)
! 33    format(1x, 5(2x, E10.4))
    
    print*, 'dust fraction ', sum(ataup(nv-3:nv))
    print*, 'smk fraction ', sum(ataup(nv-11: nv-4))
    
     open(9, file = 'forcing.dat'  )       
     do ii = 1, ntau 
      do jj = 1, nratio 
       do i  = 1, nv
         tmptaup(i) = 0.0
         if ( i >= nv-3 ) then    ! dust 
          tmptaup(i) = tau(ii)* aotratio(jj) * ataup(i)/sum(ataup(nv-3:nv))
	 else if ( i <= nv-4 .and. i >= nv - 11 ) then ! smoke
	  tmptaup(i) = tau(ii)* (1-aotratio(jj)) *  ataup(i)/sum(ataup(nv-11:nv-4))
	 endif

       enddo 
        
       print*, 'tmptaup = ', tmptaup(1:nv)
       print*, 'totaltau = ', sum(tmptaup(1:nv))

! ADD SHORTWAVE PART
! FULIOU CODE CONVERT CLOUDTAU TO CLDWT 
!       CLDOPT = 20.0
!       cddwt(10) = 20.0/


       call rad (nv, nv1, ndfs, mdfs,  &
                 pp, pt, ph, po, tmptaup, atype, &
                 cldwt, cldwtreff, cldice, cldicereff, rainwt, graupwt, &
		 as, ee, pts, u0, ss,  &
                 umco2, umch4, umn2o, & 
                 fds, fus, fdir, fuir)      
  
        write(9,'(5(2x, f10.4))')  tau(ii), fds(nv1), fus(nv1), fds(1),	fus(1)  
!	stop
    enddo
    enddo
     close(9)
     
	end
	
	
	



