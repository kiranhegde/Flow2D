!============================
subroutine write_result(iter,cfl,cl,cd,res)
!============================
use flux
use inf
use pri
use vis
use st_end
implicit none

integer(kind=i2):: i,j,k,bk,jterm,iter,iside  
integer(kind=i2),parameter::eight=8 
real(kind=dp)::pmin,mmin,rmin,umin,vmin,emin,aoa1,mumin 
real(kind=dp)::pmax,mmax,rmax,umax,vmax,emax,x1 ,mumax
real(kind=dp):: cfl,qtemp(2+ndim),mach,res,ent,p_p,cl,cd
real(kind=dp):: cp,cfx,cfy,cm,sx,sy,ds,pw,p1,var(8),mu
real(kind=dp):: amx,amy,gamma_m1,gterm2,vel_nor,u2,v2,c2,p2,term4,term5,term6

k=1
!pi = 4.0d0*datan(1.d0)
aoa1=aoa*pi/180.d0
pmin =  1.0d10
pmax = -1.0d10
mmin =  1.0d10
mmax = -1.0d10
rmin =  1.0d10
rmax = -1.0d10
umin =  1.0d10
umax = -1.0d10
vmin =  1.0d10
vmax = -1.0d10
emin =  1.0d10
emax = -1.0d10
mumin =  1.0d10
mumax = -1.0d10

if(ieuler .eq. 0) then
do bk = 1,nblk
   call start_end(bk)
do j=1,cy
  do i=1,cx
     qtemp(:)=blk(bk)%fv(i,j,k)%qvar(:)
     call q_conv(qtemp)
     call viscous(t,mu)

     !p_p   = 2*(p-znd)
     !p = p_p!-pinf
     mach=dsqrt(u*u+v*v)/a
     ent  = dlog10(p/rho**gamma/ent_inf)

     pmin=dmin1(pmin,p)  ; pmax=dmax1(pmax,p)
     mmin=dmin1(mmin,mach) ; mmax=dmax1(mmax,mach)
     rmin=dmin1(rmin,rho)  ; rmax=dmax1(rmax,rho)
     umin=dmin1(umin,u  )  ; umax=dmax1(umax,u  )
     vmin=dmin1(vmin,v  )  ; vmax=dmax1(vmax,v  )
     emin=dmin1(emin,ent)  ; emax=dmax1(emax,ent)
     mumin=dmin1(mumin,mu)  ; mumax=dmax1(mumax,mu)
     !emin=dmin1(emin,e)  ; emax=dmax1(emax,e)
  enddo
enddo
enddo
else
do bk = 1,nblk
   call start_end(bk)
do j=1,cy
  do i=1,cx
     qtemp(:)=blk(bk)%fv(i,j,k)%qvar(:)
     call q_conv(qtemp)

     !p_p   = 2*(p-znd)
     !p = p_p!-pinf
     mach=dsqrt(u*u+v*v)/a
     ent  = dlog10(p/rho**gamma/ent_inf)

     pmin=dmin1(pmin,p)  ; pmax=dmax1(pmax,p)
     mmin=dmin1(mmin,mach) ; mmax=dmax1(mmax,mach)
     rmin=dmin1(rmin,rho)  ; rmax=dmax1(rmax,rho)
     umin=dmin1(umin,u  )  ; umax=dmax1(umax,u  )
     vmin=dmin1(vmin,v  )  ; vmax=dmax1(vmax,v  )
     emin=dmin1(emin,ent)  ; emax=dmax1(emax,ent)
     mumin= 0.d0  ; mumax=0.d0
     !emin=dmin1(emin,e)  ; emax=dmax1(emax,e)

  enddo
enddo
enddo
endif

       M_max=mmax
       M_min=mmin

      write(*,9)('-',i=1,70)
9     format(70a)

      write(*,10)m_inf,aoa,re_l,cfl
10    format(' mach =',f6.3,', aoa =',f6.2, ', rey = ',e15.4,', cfl =',f8.2)
!      write(*,11)iflux,ilimit,gridfile
!11    format(' flux =',i2, ',     lim = ',i2,',    grid= ',a30)
      write(*,9)('-',i=1,70)
      write(*,'(" iterations        =",i12)')iter
!      write(*,'(" global dt         =",e16.6)') dtglobal
!      write(*,'(" l2 residue        =",e16.6)')dlog10(res)
      write(*,'(" l2 residue        =",e16.6)')res
!      write(*,'(" linf residue      =",e16.6)')fresi
!      write(*,'(" linf triangle     =",i12)') iresi
      write(*,*)
      write(*,'(" cl, cd            =",2f12.6)')cl,cd
      write(*,*)
      write(*,'(27x,"min",8x,"max")')          
      write(*,'(" density           =",2f12.6)')rmin, rmax
      write(*,'(" pressure          =",2f12.6)')pmin, pmax
      write(*,'(" mach number       =",2f12.6)')mmin, mmax
      write(*,'(" x velocity        =",2f12.6)')umin, umax
      write(*,'(" y velocity        =",2f12.6)')vmin, vmax
      write(*,'(" entropy           =",2f12.6)')emin, emax
      write(*,'(" Laminar Viscosity =",2f12.6)')mumin, mumax
!      write(*,'(" energy            =",2f12.6)')emin, emax
      write(*,9)('-',i=1,70)
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*)
!      call flush(6)

!.........write the soln
  
  do bk=1,nblk
     call blockinterface(bk)
     call blockinterface(bk)
  enddo
  open(unit=12,file='q_o.dat',form='unformatted')
  write(12) iter
  do bk = 1,nblk
     call start_end(bk)
     write(12) (((blk(bk)%fv(i,j,k)%qvar(jterm),i=1,cx),j=1,cy),jterm=1,ndim+2)
  enddo
  close(12)

iside = 1

  open(unit=12,file='p_tar')
  do bk = 1,nblk
     call start_end(bk)
     j=js
     cfx=0.d0
     cfy=0.d0
     cm=0.d0
  do i =is,ie
     !qtemp(:)=blk(bk)%fv(i,j,k)%qvar(:)
     !call q_conv(qtemp)
     !pw=p
     call pseudo_Laplacian(bk,i,j,var,eight)
     p1=var(4) 
     call pseudo_Laplacian(bk,i+one,j,var,eight)
     p2=var(4) 
     pw=0.5*(p1+p2)

     cp=(pw-pinf)/pdyna
     sx =-blk(bk)%eta(i  ,j,k)%nx
     sy =-blk(bk)%eta(i  ,j,k)%ny
     x1=blk(bk)%mesh(i,j,k)%x
 
     cfx=cfx+cp*sx
     cfy=cfy+cp*sy

     ds=0.5d0*(blk(bk)%mesh(i+one,j,k)%x+blk(bk)%mesh(i,j,k)%x) 
     cm=cm+cp*sx*(ds-0.25d0)
     write(12,100) i,x1,-cp,-cp*pdyna/(pmax-pinf)
  enddo
  enddo
  close(12)
  
  !cl =(-dsin(aoa1)*xf + dcos(aoa1)*yf)/pdyna
  !cd =( dcos(aoa1)*xf + dsin(aoa1)*yf)/pdyna
  cl= cfy*dcos(aoa1)-cfx*dsin(aoa1)
  cd= cfy*dsin(aoa1)+cfx*dcos(aoa1)
 
  print*,'Pst=',pdyna+pinf 

  write(*,102)cfx,cfy
  write(*,101)cl,cd,cm





  100 format(1x,i6,2x,f10.6,2x,f10.5,2x,f10.5,2x,f10.5)
  101 format(1x,'#',1x,3(f15.6,2x))
  102 format(1x,'#',1x,2(f15.6,2x))
!==========================
end subroutine write_result
!==========================

!.....calculate lift and drag coefficients
!=====================
subroutine clcd(cl,cd)
!=====================
use flux
use pri
use inf
use st_end
use vis
implicit none

integer(kind=i2):: i,j,k,bk,ib_m,jb_m,icell,jcell,idir,iside
integer(kind=i2),parameter::eight=8 
real(kind=dp):: cl,cd,clv,cdv, xf, yf,q_tmp(2+ndim),xn,yn,vrel 
real(kind=dp):: ds,sx,sy,aml(ndim),amr(ndim),qleft(2+ndim),qright(2+ndim) 
real(kind=dp):: tau_xx,tau_xy,tau_yx,tau_yy,com_fact,visc,sfx,sfy 
real(kind=dp):: p1,p2,p3,pw,un2,un,vn,vn2,r1,r2,r3,r4,var(8)
real(kind=dp):: dn,dn2,dudn, dvdn, dudx, dudy, dvdx, dvdy, sxx,syy,sxy
real(kind=dp):: tauw,tauw1,tauw2,tx,ty,sfc,ut,udotn,udotn2,ubarn,ubarn2,sgnm,x1,y1 
real(kind=dp):: amx,amy,gamma_m1,gterm2,vel_nor,u2,v2,c2,term4,term5,term6
real(kind=dp) :: u_x,u_y,v_x,v_y,t_x,t_y,aoa1

k=1
iside = 1

aoa1=aoa*pi/180.d0
cl=0.0d0
cd=0.0d0
xf = 0.0d0
yf = 0.0d0
do bk = 1,nblk
call start_end(bk)
j=js

do i=is,ie

   q_tmp(:) = blk(bk)%fv(i,j,k)%qvar(:)
   call q_conv(q_tmp)
   c2=a
   u2=u
   v2=v
   p2=p

   sx =-blk(bk)%eta(i  ,j,k)%nx
   sy =-blk(bk)%eta(i  ,j,k)%ny
   ds = dsqrt(sx*sx+sy*sy)
   amx= sx/ds
   amy= sy/ds
   
   gamma_m1 = gamma-1
   gterm2 = gamma/gamma_m1

!  calculation of the interior riemann state calculating the normal velocity
   vel_nor = u2*amx + v2*amy
   vel_nor = -vel_nor*(-1)**iside
!
!  calculates wall pressure by solving 1d riemann problem
!
   term4 = gamma_m1/2.d0*(vel_nor/c2)
   term5 = 1.d0 - term4
   term6 = term5**(2.d0*sngl(gterm2))
   pw    = p2*term6
   if(pw .lt. 0.d0) pw = p2
   !pw = p2
   !pw=2*(pw-pinf)
   !pw=pw-pinf
   !pw = 0.5d0*(3.0d0*p2-p3)
   !pw = p2

   call pseudo_Laplacian(bk,i,j,var,eight)
   p1=var(4) 
   call pseudo_Laplacian(bk,i+one,j,var,eight)
   p2=var(4) 
   pw=0.5*(p1+p2)
   !pw=2*(pw-znd)
   !write(52,210)blk(bk)%mesh(i,j,k)%x,p1,pw,p2 

   xf = xf + pw*sx
   yf = yf + pw*sy

enddo
enddo

!print 209,'x,y Pressure Forces',xf,yf
      print 211,xf,yf

sfx=0.0d0
sfy=0.0d0

if(ieuler == 0 ) then 
!if ( turb /= 0 .or. )  then
open(7,file='sfc.dat')
idir=2
do bk = 1,nblk
call start_end(bk)
j=js
do i=is,ie

   icell = i
   jcell = j 
   q_tmp(:) = blk(bk)%fv(i,j,k)%qvar(:)
   call q_conv(q_tmp)
   c2=a
   u2=u
   v2=v
   p2=p

   sx =blk(bk)%eta(i  ,j,k)%nx
   sy =blk(bk)%eta(i  ,j,k)%ny
   ds =dsqrt(sx*sx+sy*sy)
   amx=sx/ds
   amy=sy/ds
   
   gamma_m1 = gamma-1
   gterm2 = gamma/gamma_m1
!  calculation of the interior riemann state
!  calculating the normal velocity
!
   vel_nor = u2*amx + v2*amy
   vel_nor = -vel_nor*(-1)**iside
!
!  calculates wall pressure by solving 1d riemann problem
!
   term4 = gamma_m1/2.d0*(vel_nor/c2)
   term5 = 1.d0 - term4
   term6 = term5**(2.d0*sngl(gterm2))
   !pw = p2*term6
   !if(pw .lt. 0.d0) pw = p2


   !pw = p2
   !pw=2*(pw-znd)
   !call pseudo_Laplacian(bk,i,j,var,1)
   !p1=var(1) 
   !call pseudo_Laplacian(bk,i+one,j,var,1)
   !p2=var(1) 
   !pw=0.5*(p1+p2)
   

   aml(1)  = blk(bk)%eta(i  ,j,k)%nx
   aml(2)  = blk(bk)%eta(i  ,j,k)%ny
   amr(1)  = blk(bk)%eta(i,j+one,k)%nx
   amr(2)  = blk(bk)%eta(i,j+one,k)%ny

   qleft(:)   =blk(bk)%fv(i,j,k)%qvar(:) 
   qright(:)  =blk(bk)%fv(i,j+one,k)%qvar(:) 

   !call deriv_int_eta(bk,icell,jcell,idir,qleft,qright,amr,aml)
  
   u_x=blk(bk)%eta(icell,jcell,k)%u_x
   u_y=blk(bk)%eta(icell,jcell,k)%u_y
   v_x=blk(bk)%eta(icell,jcell,k)%v_x
   v_y=blk(bk)%eta(icell,jcell,k)%v_y
   t_x=blk(bk)%eta(icell,jcell,k)%t_x
   t_y=blk(bk)%eta(icell,jcell,k)%t_y
                     
   !sx = aml(1)
   !sy = aml(2)
   !ds = dsqrt(sx*sx+sy*sy)
   !tx = -sy/ds
   !ty = sx/ds

   !visc=0.5d0*(mu(i,j)+mu(i,j+1))
   visc=0.5d0*(blk(bk)%fv(i,j,k)%mu+blk(bk)%fv(i,j+1,k)%mu) 
 
 
   ds = dsqrt(aml(1)**2 + aml(2)**2)
   sx = -aml(1)
   sy = -aml(2)
   tx = sy
   ty =-sx

   tauw1 = visc*((u_x*tx + v_x*ty)*sx + (u_y*tx + v_y*ty)*sy)
 
   com_fact = 2.d0*(u_x + v_y)/3.d0  ! check this quantity
   tau_xx = visc*(2.d0*u_x - com_fact)
   tau_yy = visc*(2.d0*v_y - com_fact)
   tau_xy = visc*(u_y + v_x)
   tau_yx = tau_xy

   tauw2 = tau_xy*(tx*tx-ty*ty)-(tau_xx-tau_yy)*tx*ty
   !sfc = DSIGN(1.d0,ut)*tauw/pdyna
   !write(7,*)x(i,j),y(i,j),tauw1,tauw2

   write(7,*)blk(bk)%mesh(i,j,k)%x,tauw1,tauw2

   sfx = sfx+tau_xx*sx+tau_xy*sy
   sfy = sfy+tau_xy*sx+tau_yy*sy

enddo

enddo
      !print 209,'x,y Viscous Forces',sfx,sfy
      print 212,sfx,sfy

      close(7)
endif   

      xf = xf + sfx !tau_xx*sx-tau_xy*sy
      yf = yf + sfy !tau_xy*sx-tau_yy*sy


      cl =(-dsin(aoa1)*xf + dcos(aoa1)*yf)/pdyna
      cd =( dcos(aoa1)*xf + dsin(aoa1)*yf)/pdyna
      write(29,*) cl , cd

209 format(1x,2(f15.4,1x))
210 format(1x,4(f15.9,1x))
211 format(1x,'x,y Pressure Forces:',2x,2(f15.9,1x))
212 format(1x,'x,y Viscous  Forces:',2x,2(f15.9,1x))

!==================
end subroutine clcd
!==================

subroutine  interpol(bk,i,j,vari,mm)
use dims
use flux
use st_end
implicit none
integer(kind=i2):: i,j,k,bk,mm
real(kind=dp)::xi,yi,vari(mm),var(mm),r1,r2,r3,r4,rsum 
real(kind=dp)::cell_centroid_dist,x1,y1

k=1
      vari(:)=0.d0
      xi=blk(bk)%mesh(i,j,k)%x
      yi=blk(bk)%mesh(i,j,k)%y

      r1=cell_centroid_dist(bk,i-one,j-one,xi,yi)
      if(i==2.or.j==2) r1=0.0
      call primitive(bk,i-one,j-one,var,mm)
      vari(:)=var(:)*r1

      r2=cell_centroid_dist(bk,i-one,j  ,xi,yi)
      if(i==2.or.j==cy-1) r2=0.0
      call primitive(bk,i-one,j  ,var,mm)
      vari(:)=vari(:)+var(:)*r2

      r3=cell_centroid_dist(bk,i  ,j  ,xi,yi)
      if(i==cx-1.or.j==cy-1) r3=0.0
      call primitive(bk,i  ,j  ,var,mm)
      vari(:)=vari(:)+var(:)*r3

      r4=cell_centroid_dist(bk,i  ,j-one,xi,yi)
      if(i==cx-1.or.j==2) r4=0.0
      call primitive(bk,i  ,j-one,var,mm)
      vari(:)=vari(:)+var(:)*r4

      rsum=r1+r2+r3+r4
 
      vari(:)=vari(:)/rsum

end subroutine  interpol

real(kind=dp)  function cell_centroid_dist(bk,i,j,xi,yi)
use dims
use flux
implicit none

integer(kind=i2):: i,j,bk
real(kind=dp)::xi,yi,x1,y1,dist

     call cell_centroid(bk,i,j,x1,y1)
     dist=dsqrt( (xi-x1)**2+(yi-y1)**2)
     dist=1.0/dist
     cell_centroid_dist=dist

end function  cell_centroid_dist

subroutine  cell_centroid(bk,i,j,xi,yi)
use dims
use flux
implicit none

integer(kind=i2):: i,j,k,bk,ii
real(kind=dp)::xi,yi 
real(kind=dp)::x(5),y(5),area,tmp

k=1
     x(1)=blk(bk)%mesh(i,j,k)%x
     x(2)=blk(bk)%mesh(i+one,j,k)%x
     x(3)=blk(bk)%mesh(i+one,j+one,k)%x
     x(4)=blk(bk)%mesh(i,j+one,k)%x
     x(5)=x(1)
     y(1)=blk(bk)%mesh(i,j,k)%y
     y(2)=blk(bk)%mesh(i+one,j,k)%y
     y(3)=blk(bk)%mesh(i+one,j+one,k)%y
     y(4)=blk(bk)%mesh(i,j+one,k)%y
     y(5)=y(1)

     area=0.d0
     xi=0.d0
     yi=0.d0
     do ii=1,4
        tmp=x(ii)*y(ii+1)-x(ii+1)*y(ii)
        area=area+tmp
        xi=xi+tmp*(x(ii)+x(ii+1))
        yi=yi+tmp*(y(ii)+y(ii+1))
     enddo
     area=0.5d0*area
     xi=xi/6.d0/area
     yi=yi/6.d0/area
     !x1=0.25*(x(1)+x(2)+x(3)+ x(4))
     !y1=0.25*(y(1)+y(2)+y(3)+ y(4))

end subroutine cell_centroid


subroutine  pseudo_Laplacian(bk,i,j,vari,mm)
use dims
use st_end
use flux
implicit none

integer(kind=i2):: i,j,k,bk,ii,mm
real(kind=dp)::x(4),y(4),area,tmp,x1,y1,cell_centroid_dist
real(kind=dp)::xi,yi,vari(mm),var(mm),r1,r2,r3,r4,rsum 
real(kind=dp)::ixx,iyy,ixy,rx,ry,gx,gy,dx,dy,det 

k=1
vari(:)=0.d0
xi=blk(bk)%mesh(i,j,k)%x
yi=blk(bk)%mesh(i,j,k)%y
ixx=0.d0 ; iyy=0.d0
ixy=0.d0
rx=0.d0 ; ry=0.d0
gx=0.d0 ; gy=0.d0


   call cell_centroid(bk,i-one,j-one,x1,y1)
   dx=x1-xi ; dy=y1-yi
   rx=dx    ; ry=dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(1)=dx  ; y(1)=dy

   call cell_centroid(bk,i-one,j  ,x1,y1)
   dx=x1-xi ; dy=y1-yi
   rx=rx+dx    
   ry=ry+dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(2)=dx  ; y(2)=dy

   call cell_centroid(bk,i  ,j  ,x1,y1)
   dx=x1-xi ; dy=y1-yi 
   rx=rx+dx    
   ry=ry+dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(3)=dx  ; y(3)=dy

   call cell_centroid(bk,i  ,j-one,x1,y1)
   dx=x1-xi ; dy=y1-yi
   rx=rx+dx    
   ry=ry+dy   
   ixx=ixx+dx*dx
   iyy=iyy+dy*dy
   ixy=ixy+dx*dy
   x(4)=dx  ; y(4)=dy

   det=ixx*iyy-ixy*ixy 
   if(det==0.d0) then
      print*,'det=',det
   endif
   gx=(ixy*ry-iyy*rx)/det
   gy=(ixy*rx-ixx*ry)/det
 
   call start_end(bk)

   r1=1+gx*x(1)+gy*y(1)
   call primitive(bk,i-one,j-one,var,mm)
   vari(:)=var(:)*r1

   r2=1+gx*x(2)+gy*y(2)
   call primitive(bk,i-one,j  ,var,mm)
   vari(:)=vari(:)+var(:)*r2

   r3=1+gx*x(3)+gy*y(3)
   call primitive(bk,i  ,j  ,var,mm)
   vari(:)=vari(:)+var(:)*r3

   r4=1+gx*x(4)+gy*y(4)
   call primitive(bk,i  ,j-one,var,mm)
   vari(:)=vari(:)+var(:)*r4

   rsum=r1+r2+r3+r4
   vari(:)=vari(:)/rsum

end subroutine  pseudo_Laplacian

subroutine primitive(bk,i,j,var,mm)
use dims
use inf
use pri
use vis
use flux
implicit none
integer(kind=i2):: i,j,k,bk,mm
real(kind=dp)::var(mm),q_tmp(ndim+2),pr,pm

k=1
      q_tmp(1:ndim+2) = blk(bk)%fv(i,j,k)%qvar(1:ndim+2)
      call q_conv(q_tmp)
      !pr   = p/znd*p_inf
      !pr   = 2*(p-znd)
      pm  = dsqrt(u*u+v*v)/a
      var(1)=rho
      var(2)=u
      var(3)=v
      var(4)=p
      var(5)=t
      var(6)=pm
      if(ieuler == 0 .and. turb == 0 )  var(7)=blk(bk)%fv(i,j,k)%mu
      if(ieuler == 0 .and. turb == 1 )  var(8)=blk(bk)%fv(i,j,k)%mut

end subroutine primitive
