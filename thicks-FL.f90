	subroutine thicks (pp, pt, nv, dz)
! *********************************************************************
! dz is the thickness of a layer in units of km.
! *********************************************************************
        implicit none
	include 'Rad4S.h'
	
! input
        real, dimension(nv1x):: pp, pt
	integer nv

! output
        real dz(nvx)

! internal
        integer i
			
    dz(:) = 0.0
	do 100 i = 1, nv
	   dz(i) = 0.0146337 * ( pt(i) + pt(i+1) ) & 
     	   * alog( pp(i+1) / pp(i) )
100	continue


        return
	end
