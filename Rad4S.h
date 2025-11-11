! 
!4S header files
!
  integer, parameter :: nvx = 40, nv1x = nvx + 1
  integer, parameter:: ndfsx = nvx, mdfsx = nvx + 1, ndfs4x = 4*ndfsx
  integer, parameter:: mb = 18, mbir = 12 , nc=8, mbs = 6 
  integer, parameter :: naer=1, nreff = 13
  
  real, dimension(nreff) :: Reff
  
  data  Reff  /0.10, 0.12, 0.15, 0.22, 0.25, 0.28, 0.32, 0.35, 0.38, 0.42, 0.45, 0.48, 0.50/
