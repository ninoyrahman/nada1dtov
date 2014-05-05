  SUBROUTINE Schw_isotrop_solution
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
!  USE tmp_mdmetric_u
  USE tmp_mdderivatives
  USE MD_Subroutines
  USE MD_Boundary
  IMPLICIT NONE

  INTEGER :: i,k,m,index
  REAL(kind=double) :: dr,ggrr,betar,radius
  REAL(kind=double),DIMENSION(0:10000):: r_iso,Rtmp,alp_iso,ggrr_iso,betar_iso,psi_iso

  Rtmp(0) =1.495

!   nr at param
   dr = 0.01

! compute analytic solution by Baumgarte & Naculich 

   write(*,*) "Mbh",Mbh

 do i=1,10000

   Rtmp(i) = Rtmp(i-1)+dr

   
   r_iso(i) = ((2.d0*Rtmp(i)+Mbh+(4.d0*Rtmp(i)**2+4.d0*Mbh*Rtmp(i)+3.d0*Mbh**2)**0.5)*0.25d0)* & 
            & ( (4.d0+3.d0*sqrt(2.d0))*(2.d0*Rtmp(i)-3.d0*Mbh)/(8.d0*Rtmp(i)+6.d0*Mbh+3.d0* & 
            & (8.d0*Rtmp(i)**2+8.d0*Mbh*Rtmp(i)+6.d0*Mbh**2)**0.5))**(1.d0/sqrt(2.d0)) 
   
   psi_iso(i) = ((4.d0*Rtmp(i)/(2.d0*Rtmp(i)+Mbh+(4.d0*Rtmp(i)**2+4.d0*Mbh*Rtmp(i)+3.d0*Mbh**2)**0.5))**0.5)*& 
              & ( (8.d0*Rtmp(i)+6.d0*Mbh+3.d0*(8.d0*Rtmp(i)**2.d0+8.d0*Mbh*Rtmp(i)+6.d0*Mbh**2)**0.5)/ & 
              & ((4.d0+3.d0*sqrt(2.d0))*(2.d0*Rtmp(i)-3.d0*Mbh)))**(1.d0/(2.d0*sqrt(2.d0)))

   ggrr_iso(i) = psi_iso(i)**4

   alp_iso(i) = sqrt(1.d0 -2.d0*Mbh/Rtmp(i)+(27.d0*Mbh**4.d0)/(16.d0*Rtmp(i)**4.d0))

   betar_iso(i)= 3.d0*sqrt(3.d0)*Mbh**2*r_iso(i)/(4.d0*Rtmp(i)**3.d0)


  end do

   write(*,*) 'Rtmp,r_iso,psi_iso,alp_iso', Rtmp(1),r_iso(1),psi_iso(1),alp_iso(1)
   write(*,*) 'Rtmp,r_iso,psi_iso,1+M/2r,alp_iso', Rtmp(10000),r_iso(10000),psi_iso(10000)!,1.d0+1.d0/(2.d0*r_iso(10000)),alp_iso(10000)

! interpolation 


  DO k=kmin,kmax+ghzz
    DO i=imin,imax+ghzx
    
      radius = sqrt(xh(i)**2 + zh(k)**2)

         do m = 1,10000
           
           if (r_iso(m).gt.radius) then
             index = m
             exit
           end if
             
         end do


         psi(i,k) = psi_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
                    & psi_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))

!         write(*,*) psi(i,k),psi_iso(index-1),psi_iso(index)
!         write(*,*) radius,r_iso(index),index
!         stop 

         alp(i,k) = alp_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
                    & alp_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))

         ggrr = psi(i,k)**4

!         ggrr = ggrr_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
!                    & ggrr_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))
         
         gxx(i,k) = ggrr
         gyy(i,k) = ggrr
         gzz(i,k) = ggrr

         gxy(i,k) = 0.d0
         gxz(i,k) = 0.d0
         gyz(i,k) = 0.d0


         betar = betar_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
                    & betar_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))
         
         betaux(i,k) = xh(i)*betar/radius
         betauy(i,k) = 0.d0
         betauz(i,k) = zh(k)*betar/radius


! compute the covariant shift:

         betalx(i,k) = gxx(i,k)*betaux(i,k) !+ gxz(i,k)*betauz(i,k) 
         betaly(i,k) = 0.d0!gxy(i,k)*betaux(i,k) + gyy(i,k)*betauy(i,k) + gyz(i,k)*betauz(i,k) 
         betalz(i,k) = gzz(i,k)*betauz(i,k) !+ gxz(i,k)*betaux(i,k)


         Bx(i,k) = 0.d0
         By(i,k) = 0.d0
         Bz(i,k) = 0.d0
         

    ENDDO
  ENDDO


      phi = log(psi)


      xsi = psi**(-4.0)
      capw = EXP(-2.d0*phi)

 
  DO k=kmin,kmax+ghzz
    DO i=imin,imax+ghzx

      if ((capw(i,k).lt.0.001).or.(xsi(i,k).lt.1.d-4)) then 
!         write(*,*)  "lesss", i,k,
         capw(i,k) = 0.001d0
         phi(i,k)  = -log(capw(i,k))*0.5
 
        xsi(i,k) = 1.d-4!EXP(-4.d0*phi(i,k))

!         phi(i,k) = -log(xsi(i,k))/4.d0
         psi(i,k)  = EXP(phi(i,k))
         gxx(i,k) = (psi(i,k))**4
         gyy(i,k) = gxx(i,k)
         gzz(i,k) = gxx(i,k)

!         stop
      end if
 

    ENDDO
  ENDDO


  epsxsi = 0.5*(sqrt(xh(1)**2 + zh(1)**2)/(2.d0*Mbh))**4
  write(*,*) 'floor for xsi-->',epsxsi,xsi(0,0),xsi(1,1)


   call Boundary_xaxis_func(alp)
   call Boundary_zaxis_func(alp)

   call Boundary_xaxis_func(psi)
   call Boundary_zaxis_func(psi)

   call Boundary_xaxis_func(phi)
   call Boundary_zaxis_func(phi)

   call Boundary_xaxis_func(xsi)
   call Boundary_zaxis_func(xsi)

   call Boundary_xaxis_func(capw)
   call Boundary_zaxis_func(capw)

   call Boundary_xaxis_vector(betaux,betauy,betauz)
   call Boundary_zaxis_vector(betaux,betauy,betauz)
 
   call Boundary_xaxis_tensor(gxx,gxy,gxz,gyy,gyz,gzz)
   call Boundary_zaxis_tensor(gxx,gxy,gxz,gyy,gyz,gzz)
   

!  trK = 0.d0

! storing initial = analytic data

!!$        gaxx = gxx
!!$        gaxy = gxy
!!$        gaxz = gxz
!!$        gayy = gyy  
!!$        gayz = gyz
!!$        gazz = gzz
!!$        psia2 = psi

!!$    phi_o = 0.d0
!!$    alp_o = 1.d0
!!$    trK_o = 0.d0
!!$    xsi_o = 1.d0
!!$
!!$    betaux_o = 0.d0
!!$    betauy_o = 0.d0 
!!$    betauz_o = 0.d0
!!$
!!$    gphixx_o = 1.d0
!!$    gphixy_o = 0.d0
!!$    gphixz_o = 0.d0
!!$    gphiyy_o = 1.d0
!!$    gphiyz_o = 0.d0
!!$    gphizz_o = 1.d0
!!$
!!$    Aphixx_o = 0.d0
!!$    Aphixy_o = 0.d0
!!$    Aphixz_o = 0.d0
!!$    Aphiyy_o = 0.d0
!!$    Aphiyz_o = 0.d0
!!$    Aphizz_o = 0.d0
!!$
!!$    CCFx_o = 0.d0
!!$    CCFy_o = 0.d0
!!$    CCFz_o = 0.d0


      END SUBROUTINE Schw_isotrop_solution

















