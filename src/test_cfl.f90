program cfl0
implicit none

integer*4::i,j,m
real*8:: cfl,cfl_max,fn,nu


cfl_max=1000
open(3,file='cfl.dat')
do i=1,1000
   !cfl = 0.5d0*cfl_max*(1+dtanh(i*7.d0/500.0 - 3.5))
   !nu=i*7.d0/500.0 - 5.0
   nu=i*3.d0/500.0 - 5.0
   fn=1+dtanh(nu)
   cfl = 0.5d0*cfl_max*fn
   write(3,*)i,cfl,nu,fn
enddo
close(3)

end program cfl0

