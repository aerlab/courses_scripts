
 	subroutine planck (nv1, ib, pt, pts, bf, bs )
! **********************************************************************
! bf and bs are the blackbody intensity function integrated over the
! band ib at the nv1 levels and at the surface, respectively.    The
! units of bf and bs are W/m**2/Sr. nd*10 is the band width from ve.
! **********************************************************************
	implicit none
	include 'Rad4S.h'

! input
       real, dimension(nv1x):: pt
       integer ib, nv1
       real pts

! output
       real bf(nv1x)
       real bs       	
	
	
! internal	
	real ve(mbir), nd(mbir), bt(nv1x)
	integer nv11, irb, i, j , ibr
	real bts, v1, v2, w, fq1, fq2, x
	

! initial	
	data ve / 2200.0, 1900.0, 1700.0, 1400.0, 1250.0, 1100.0,&
                 980.0, 800.0, 670.0, 540.0, 400.0, 280.001 /
	data nd / 30, 20, 30, 15, 15, 12, &
                 18, 13, 13, 14, 12, 28 /
		 
	 
        nv11 = nv1 + 1
	ibr = ib - mbs
        bts = 0.0
	do 10 i = 1, nv1
           bt(i) = 0.0
10	continue


! starts
	v1 = ve(ibr)
	do 20 j = 1, nd(ibr)
	   v2 = v1 - 10.0
	   w = ( v1 + v2 ) * 0.5
	   fq1 = 1.19107e-8 * w * w * w
	   fq2 = 1.43884 * w
	   do 30 i = 1, nv11
              if ( i .eq. nv11 ) then
	        x = fq1 / ( exp ( fq2 / pts ) - 1.0 )
	        bts = bts + x
              else
	        x = fq1 / ( exp ( fq2 / pt(i) ) - 1.0 )
                bt(i) = bt(i) + x
              endif
30	    continue
              v1 = v2
20	continue
	do 40 i = 1, nv1
          bf(i) = bt(i) * 10.0
40	continue
        bs = bts * 10.0
	return
	end
