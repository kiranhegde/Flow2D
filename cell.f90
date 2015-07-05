! cell.f.old
!===========================================
subroutine cell_area
!===========================================
use flux
use st_end
!calculates the cell area for a given cell

!===============================================================
!ndim  = space dimension (2)
!x1,x2,x3,x4 = four corner coordinates of the cell
!area = area of the cell
!===============================================================
implicit none

integer(kind=i2):: i,j,k,bk
real(kind=dp):: x1(ndim),x2(ndim),x3(ndim),x4(ndim)
real(kind=dp):: a1,a2,a3,a4

k=1
do bk=1,nblk
   call start_end(bk)
do j = 1,cy
  do i = 1,cx
    x1(1) = blk(bk)%mesh(i  ,j,k)%x
    x1(2) = blk(bk)%mesh(i  ,j,k)%y
    x2(1) = blk(bk)%mesh(i+1,j,k)%x
    x2(2) = blk(bk)%mesh(i+1,j,k)%y
    x3(1) = blk(bk)%mesh(i+1,j+1,k)%x
    x3(2) = blk(bk)%mesh(i+1,j+1,k)%y
    x4(1) = blk(bk)%mesh(i  ,j+1,k)%x
    x4(2) = blk(bk)%mesh(i  ,j+1,k)%y
     
    a1 = x1(1)*(x2(2)-x4(2))
    a2 = x2(1)*(x3(2)-x1(2))
    a3 = x3(1)*(x4(2)-x2(2))
    a4 = x4(1)*(x1(2)-x3(2))

    blk(bk)%fv(i,j,k)%vol = 0.5d0*dabs(a1 + a2 + a3 + a4)
    if(blk(bk)%fv(i,j,k)%vol==0.d0) then
       print*,'-Vol'
       print*,i,j,blk(bk)%fv(i,j,k)%vol
       stop
    endif      
  enddo
enddo
enddo


!========================
end subroutine cell_area
!========================



!=================================================
subroutine normal_zhi
!=================================================

!calculates normals on zhi faces for all the computational cells

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!x1    = x coordinates
!y1    = y coordinates
!n_zhi = normals on the xi faces
!===============================================================
use flux
use st_end
implicit none

integer(kind=i2):: i,j,k,bk
k=1
do bk=1,nblk
   call start_end(bk)
do i = 1,is-one
do j = 1,cy
  blk(bk)%zhi(i,j,k)%nx = -(blk(bk)%mesh(i,j+1,k)%y - blk(bk)%mesh(i,j,k)%y)
  blk(bk)%zhi(i,j,k)%ny =  (blk(bk)%mesh(i,j+1,k)%x - blk(bk)%mesh(i,j,k)%x)
enddo
enddo

do j = 1,cy
  do i = is,cx+one
    blk(bk)%zhi(i,j,k)%nx =  (blk(bk)%mesh(i,j+1,k)%y - blk(bk)%mesh(i,j,k)%y)
    blk(bk)%zhi(i,j,k)%ny = -(blk(bk)%mesh(i,j+1,k)%x - blk(bk)%mesh(i,j,k)%x)
  enddo
enddo
enddo
!=========================
end subroutine normal_zhi
!=========================




!=================================================
subroutine normal_eta
!=================================================

!calculates normals on eta faces for all the computational cells

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!x1    = x coordinates
!y1    = y coordinates
!n_eta = normals on the eta faces
!===============================================================
use st_end
use flux
implicit none

integer(kind=i2):: i,j,k,bk
real(kind=dp):: dx,dy,ds,x1,y1,nx,ny 

k=1
do bk=1,nblk
   call start_end(bk)
do j = 1,js-one
do i = 1,cx
    blk(bk)%eta(i,j,k)%nx =  (blk(bk)%mesh(i+1,j,k)%y - blk(bk)%mesh(i,j,k)%y)
    blk(bk)%eta(i,j,k)%ny = -(blk(bk)%mesh(i+1,j,k)%x - blk(bk)%mesh(i,j,k)%x)
enddo
enddo

do j = js,cy+one
  do i = 1,cx
    blk(bk)%eta(i,j,k)%nx = -(blk(bk)%mesh(i+1,j,k)%y - blk(bk)%mesh(i,j,k)%y)
    blk(bk)%eta(i,j,k)%ny =  (blk(bk)%mesh(i+1,j,k)%x - blk(bk)%mesh(i,j,k)%x)
  enddo
enddo
enddo

k=1
open(3,file='prof.dat')
do bk=1,nblk
   call start_end(bk)
  j=js
  do i = is,cx
    x1 = blk(bk)%mesh(i,j,k)%x
    y1 = blk(bk)%mesh(i,j,k)%y 
    write(3,101)x1,y1
  enddo
    write(3,*)

  j=je+one
  do i = is,cx
    x1 = blk(bk)%mesh(i,j,k)%x
    y1 = blk(bk)%mesh(i,j,k)%y 
    write(3,101)x1,y1
  enddo
    write(3,*)

  i=is   
  do j = js,cy
    x1 = blk(bk)%mesh(i,j,k)%x
    y1 = blk(bk)%mesh(i,j,k)%y 
    write(3,101)x1,y1
  enddo
    write(3,*)

  i=ie+one
  do j = js,cy
    x1 = blk(bk)%mesh(i,j,k)%x
    y1 = blk(bk)%mesh(i,j,k)%y 
    write(3,101)x1,y1
  enddo

enddo
close(3)

!100 format(1x,4(f15.6,1x))
101 format(1x,2(f15.9,1x))

!=========================
end subroutine normal_eta
!=========================



!==================================================
subroutine cent_fac(bk,i,j,idir,xfc,yfc)
!==================================================

!calculates the face centres of a given cell face

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!i     = xi index
!j     = eta index
!idir  = 1 for xi direction, 2 for eta direction
!x     = x coordinates
!y     = y coordinates
!xfc   = x coordinate for the face centre
!yfc   = y coordinate for the face centre
!===============================================================
use flux

implicit none

integer(kind=i2):: bk,i,j,k,idir
real(kind=dp):: xfc,yfc

k=1

select case(idir)

case(1)
  xfc = 0.5d0*(blk(bk)%mesh(i+1,j,k)%x+blk(bk)%mesh(i+1,j+1,k)%x)
  yfc = 0.5d0*(blk(bk)%mesh(i+1,j,k)%y+blk(bk)%mesh(i+1,j+1,k)%y)
case(2)
  xfc = 0.5d0*(blk(bk)%mesh(i,j+1,k)%x+blk(bk)%mesh(i+1,j+1,k)%x)
  yfc = 0.5d0*(blk(bk)%mesh(i,j+1,k)%y+blk(bk)%mesh(i+1,j+1,k)%y)
end select


!=======================
end subroutine cent_fac
!=======================


!===================================================
subroutine centroid1(bk,i,j,idir,xfc,yfc)
!===================================================

!calculates the centroid of a given cell

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!i     = xi index
!j     = eta index
!idir  = 1 for xi direction, 2 for eta direction
!xfc   = x coordinate for the centroid
!yfc   = y coordinate for the centroid
!x     = x coordinates
!y     = y coordinates
!===============================================================
use dims
use flux
implicit none
integer(kind=i2):: i,j,k,bk,idir
real(kind=dp):: xfc,yfc

integer(kind=i2):: ip1,jp1
real(kind=dp):: xa,xb,xc,xd,ya,yb,yc,yd
real(kind=dp):: third
real(kind=dp):: xc1,yc1,termx,termy,dab,dbc,dca,speri,area1
real(kind=dp):: xc2,yc2,dda,dac,dcd,area2,area
k=1
third = 1.d0/3.d0
!
!initializing the vertices for zhi direction
!
ip1 = i+one
jp1 = j+one

xa  = blk(bk)%mesh(i  ,j,k)%x
xb  = blk(bk)%mesh(i  ,jp1,k)%x
xc  = blk(bk)%mesh(ip1,jp1,k)%x
xd  = blk(bk)%mesh(ip1,j,k)%x

ya  = blk(bk)%mesh(i  ,j,k)%y
yb  = blk(bk)%mesh(i  ,jp1,k)%y
yc  = blk(bk)%mesh(ip1,jp1,k)%y
yd  = blk(bk)%mesh(ip1,j,k)%y

!
!dividing into two triangles abc  dac
!
!first triangle
!
xc1   = (xa + xb + xc)*third
yc1   = (ya + yb + yc)*third

termx = xb - xa
termy = yb - ya
dab   = dsqrt(termx*termx +  termy*termy)

termx = xc - xb
termy = yc - yb
dbc   = dsqrt(termx*termx + termy*termy)

termx = xa - xc
termy = ya - yc
dca   = dsqrt(termx*termx + termy*termy)

speri = 0.5d0*(dab + dbc + dca)
area1 = dsqrt(speri*(speri - dab)*(speri - dbc)*(speri - dca))

!
!second triangle
!
xc2   = (xd + xa + xc)*third
yc2   = (yd + ya + yc)*third

termx = xa - xd
termy = ya - yd
dda   = dsqrt(termx*termx + termy*termy)

dac   = dca

termx = xd - xc
termy = yd - yc
dcd   = dsqrt(termx*termx + termy*termy)

speri = 0.5d0*(dda + dac + dcd)
area2 = dsqrt(speri*(speri - dda)*(speri - dac)*(speri - dcd))
!
!calculation of cell face centroids
!
area = area1 + area2

if(area .ne. 0.d0) then
  xfc = (xc1*area1 + xc2*area2)/area
  yfc = (yc1*area1 + yc2*area2)/area
else
  xfc = 0.5d0*(xc1 + xc2)
  yfc = 0.5d0*(yc1 + yc2)
endif

!========================
end subroutine centroid1
!========================




!===========================================================
subroutine cell_dist(bk,cx,cy,icell,jcell,idir,iside,r1,r2,r3,r4)!,area)
!===========================================================

!for a given interface, this subroutine computes the cell center
!distances of two cells on either sides from the given interface

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!cx    = same as n for this code
!cy    = same as n for this code
!icell = xi index
!jcell = eta index
!idir  = 1 for xi direction, 2 for eta direction
!iside = 1 for left interface, 2 for right interface
!r1    = distance of the interface to the cell centre of 2nd cell
!        to the left
!r2    = distance of the interface to the cell centre of 1st cell
!        to the left
!r3    = distance of the interface to the cell centre of 1st cell
!        to the right
!r4    = distance of the interface to the cell centre of 2nd cell
!        to the right
!n_zhi = normals on the xi faces
!n_eta = normals on the eta faces
!x     = x coordinates
!y     = y coordinates
!area  = area of the cell
!===============================================================
use flux
implicit none
integer(kind=i2):: bk,k,icell,jcell,idir,iside,cx,cy
real(kind=dp):: r1,r2,r3,r4

real(kind=dp):: amx,amy,grad_met,xfc,yfc
real(kind=dp):: x1,x2,x3,x4,y1,y2,y3,y4,dx,dy
k=1

select case(idir)

case(1)
  call cent_fac(bk,icell,jcell,idir,xfc,yfc)
  amx = blk(bk)%zhi(icell+1,jcell,k)%nx
  amy = blk(bk)%zhi(icell+1,jcell,k)%ny

  grad_met = dsqrt(amx**2 + amy**2)
  amx = amx/grad_met
  amy = amy/grad_met

  select case(iside)

  case(1)

    call centroid(bk,icell+two,jcell,idir,x1,y1)
    r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)

    call centroid(bk,icell+one,jcell,idir,x2,y2)
    r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)

    call centroid(bk,icell,jcell,idir,x3,y3)
    r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)

    if (icell .ne. 1) then
      call centroid(bk,icell-one,jcell,idir,x4,y4)
      r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
    else
      r4 = 3.d0*(r2+r3) - r1
    endif

  case(2)

    call centroid(bk,icell-one,jcell,idir,x1,y1)
    r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)

    call centroid(bk,icell,jcell,idir,x2,y2)
    r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)

    call centroid(bk,icell+one,jcell,idir,x3,y3)
    r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)

    if (icell .ne. cx-one) then
      call centroid(bk,icell+two,jcell,idir,x4,y4)
      r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
    else
      r4 = 3.d0*(r2+r3) - r1
    endif

  end select

case(2)

  call cent_fac(bk,icell,jcell,idir,xfc,yfc)
  amx = blk(bk)%eta(icell,jcell+one,k)%nx
  amy = blk(bk)%eta(icell,jcell+one,k)%ny

  grad_met = dsqrt(amx**2 + amy**2)
  amx = amx/grad_met
  amy = amy/grad_met

  select case(iside)

  case(1)

    call centroid(bk,icell,jcell+two,idir,x1,y1)
    r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)

    call centroid(bk,icell,jcell+one,idir,x2,y2)
    r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)

    call centroid(bk,icell,jcell,idir,x3,y3)
    r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)

    if(jcell .ne. 1) then
      call centroid(bk,icell,jcell-one,idir,x4,y4)
      r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
    else
      r4 = 3.d0*(r2+r3) - r1
    endif

  case(2)

    call centroid(bk,icell,jcell-one,idir,x1,y1)
    r1 = dabs((x1-xfc)*amx + (y1-yfc)*amy)

    call centroid(bk,icell,jcell,idir,x2,y2)
    r2 = dabs((x2-xfc)*amx + (y2-yfc)*amy)

    call centroid(bk,icell,jcell+one,idir,x3,y3)
    r3 = dabs((x3-xfc)*amx + (y3-yfc)*amy)

    if(jcell .ne. cy-one) then
      call centroid(bk,icell,jcell+two,idir,x4,y4)
      r4 = dabs((x4-xfc)*amx + (y4-yfc)*amy)
    else
      r4 = 3.d0*(r2+r3) - r1
    endif

  end select

end select

!========================
end subroutine cell_dist
!========================




!==================================================
subroutine centroid(bk,i,j,idir,xfc,yfc)
!==================================================

!simple centroid calculation

!===============================================================
!ndim  = space dimension (2)
!n     = total no. of cells in xi direction
!m     = total no. of cells in eta direction
!i     = xi index
!j     = eta index
!idir  = 1 for xi direction, 2 for eta direction
!xfc   = x coordinate for the centroid
!yfc   = y coordinate for the centroid
!x     = x coordinates
!y     = y coordinates
!===============================================================
use flux

implicit none
integer(kind=i2):: i,j,k,bk,idir
real(kind=dp):: xfc,yfc

integer(kind=i2):: ip1,jp1
real(kind=dp):: xa,xb,xc,xd,ya,yb,yc,yd
real(kind=dp):: xac,yac,xbd,ybd
k=1
ip1 = i+one
jp1 = j+one

xa  = blk(bk)%mesh(i  ,j,k)%x
xb  = blk(bk)%mesh(i  ,jp1,k)%x
xc  = blk(bk)%mesh(ip1,jp1,k)%x
xd  = blk(bk)%mesh(ip1,j,k)%x

ya  = blk(bk)%mesh(i  ,j,k)%y
yb  = blk(bk)%mesh(i  ,jp1,k)%y
yc  = blk(bk)%mesh(ip1,jp1,k)%y
yd  = blk(bk)%mesh(ip1,j,k)%y

xac = 0.5d0*(xa + xc)
yac = 0.5d0*(ya + yc)

xbd = 0.5d0*(xb + xd)
ybd = 0.5d0*(yb + yd)

xfc = 0.5d0*(xac + xbd)
yfc = 0.5d0*(yac + ybd)

!=======================
end subroutine centroid
!=======================
