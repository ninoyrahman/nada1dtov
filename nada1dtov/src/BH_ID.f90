
  Subroutine BH_ID
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdhydro
  USE tmp_mdmetric
  USE MD_Boundary

  implicit none

  REAL(KIND=double) :: Mbh,radius,dr
  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: phi_iso,xsi2_iso,error_xsi2
  REAL(KIND=double), DIMENSION(1:100000) :: psi_iso,r_iso,Rtmp,ggrr_iso
  REAL(KIND=double), DIMENSION(1:100000) :: alp_iso,betar_iso
  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: psi_iso2

  INTEGER :: i,index,m

  Mbh = 1.d0

100 FORMAT (2E22.11)
110 FORMAT (A9,E22.11)

  dr = 0.001

   do i=1,100000
 
        Rtmp(i) = Rtmp(i-1)+dr

   
   r_iso(i) = ((2.d0*Rtmp(i)+Mbh+(4.d0*Rtmp(i)**2+4.d0*Mbh*Rtmp(i)+3.d0*Mbh**2)**0.5)*0.25d0)* & 
            & ( (4.d0+3.d0*sqrt(2.d0))*(2.d0*Rtmp(i)-3.d0*Mbh)/(8.d0*Rtmp(i)+6.d0*Mbh+3.d0* & 
            & (8.d0*Rtmp(i)**2+8.d0*Mbh*Rtmp(i)+6.d0*Mbh**2)**0.5))**(1.d0/sqrt(2.d0)) 
   
   psi_iso(i) = ((4.d0*Rtmp(i)/(2.d0*Rtmp(i)+Mbh+(4.d0*Rtmp(i)**2+4.d0*Mbh*Rtmp(i)+3.d0*Mbh**2)**0.5))**0.5)*& 
              & ( (8.d0*Rtmp(i)+6.d0*Mbh+3.d0*(8.d0*Rtmp(i)**2.d0+8.d0*Mbh*Rtmp(i)+6.d0*Mbh**2)**0.5)/ & 
              & ((4.d0+3.d0*sqrt(2.d0))*(2.d0*Rtmp(i)-3.d0*Mbh)))**(1.d0/(2.d0*sqrt(2.d0)))


   alp_iso(i) = sqrt(1.d0 -2.d0*Mbh/Rtmp(i)+(27.d0*Mbh**4.d0)/(16.d0*Rtmp(i)**4.d0))

   betar_iso(i)= 3.d0*sqrt(3.d0)*Mbh**2*r_iso(i)/(4.d0*Rtmp(i)**3.d0)



!!$   phi_iso(i) = log(psi_iso(i))
!!$
!!$
!!$   xsi2_iso(i) = EXP(-2.d0*phi_iso(i))
!!$
!!$   error_xsi2(i) = abs(xsi2(i)-xsi2_iso(i))

   end do


  DO i=imin,imax+3
    
      radius = xh(i)

            do m = 1,100000
           
           if (r_iso(m).gt.radius) then
             index = m
             exit
           end if
             
         end do


    psi_iso2(i) = psi_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
        & psi_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))


    alp(i) = alp_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
        & alp_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))


    betaux(i) = betar_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + & 
        & betar_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))


   xsi(i) = log(psi_iso2(i))

   xsi2(i) = EXP(-2.d0*xsi(i))

   grr(i) = xsi(i)**4


      small_a(i) = 1.d0

      small_b(i) = 1.d0

      A_a(i) = 0.d0


      trK(i) = 0.d0

      delta_x(i) = 0.d0

      rho_matter(i) = 0.d0

      S_a(i) = 0.d0

      S_b(i) = 0.d0

      j_r(i)  = 0.d0

      capB(i) = 0.d0

      RHSdeltar(i) = 0.d0

   

  end DO



!!$  call Derivatives_func_12(betaux,dxdbetaux,dxxdbetaux)
!!$
!!$  A_a = (dxdbetaux+dxdbetaux-2.d0*dxdbetaux/3.d0)/(2.d0*alp) !OK



  call Boundary_zaxis_func(xsi2) 
  call Boundary_zaxis_func(xsi) 
  call Boundary_zaxis_func(small_a) 
  call Boundary_zaxis_func(small_b) 
  call Boundary_zaxis_func(trK) 
  call Boundary_zaxis_func(A_a) 
  call Boundary_zaxis_asymfunc(delta_x) 
  call Boundary_zaxis_func(alp)
  call Boundary_zaxis_asymfunc(betaux) 




  end Subroutine BH_ID
