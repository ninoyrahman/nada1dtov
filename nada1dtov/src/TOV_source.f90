  SUBROUTINE TOV_source
  USE tmp_mdparam 
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE tmp_mdhydro
  USE MD_Boundary

  IMPLICIT NONE

  INTEGER :: i,index,m


!   REAL(KIND=double), DIMENSION(1:4000) :: radius,press1,rho1,eps_r,alp1,g_rr1
   REAL(KIND=double), DIMENSION(1:80000) :: radius,press1,rho1,eps_r,alp1,g_rr1

   REAL(KIND=double) :: det

  
   kp = 100.d0
   gamma = 2.d0



!   OPEN(UNIT=300, FILE='TOV_r.dat',STATUS='UNKNOWN',ACTION= 'READ')

   OPEN(UNIT=300, FILE='TOV_r2.dat',STATUS='UNKNOWN',ACTION= 'READ') !stable star

!   OPEN(UNIT=300, FILE='TOV_r3.dat',STATUS='UNKNOWN',ACTION= 'READ') !unstable star collpase

!   OPEN(UNIT=300, FILE='TOV_r4.dat',STATUS='UNKNOWN',ACTION= 'READ') !unstable star collpase

   OPEN(UNIT=301, FILE='data/TOV_source1.dat',STATUS='UNKNOWN',ACTION= 'WRITE')

   OPEN(UNIT=302, FILE='data/TOV_data.dat',STATUS='UNKNOWN',ACTION= 'WRITE')


!    do index = 1,4000
    do index = 1,80000
     
     read(300,*) radius(index),press1(index),rho1(index),eps_r(index),alp1(index),g_rr1(index)
!     write(*,*) index, radius(index)

    end do

    close (unit=300)


   rho_atm = maxval(rho1)*1.d-10


      write(302,*) imax
 
    DO i=imin-ghzx,imax+ghzx


!         do m = 1,4000
         do m = 2,80000
           
           if (radius(m).gt.xh(i)) then
             index = m+1
             exit
           end if
             
         end do


          p(i) = press1(index-1)*(xh(i)-radius(index))/(radius(index-1)-radius(index)) + & 
                    & press1(index)*(xh(i)-radius(index-1))/(radius(index)-radius(index-1))


         rho(i) = rho1(index-1)*(xh(i)-radius(index))/(radius(index-1)-radius(index)) + & 
                    & rho1(index)*(xh(i)-radius(index-1))/(radius(index)-radius(index-1))

         alp(i) = alp1(index-1)*(xh(i)-radius(index))/(radius(index-1)-radius(index)) + & 
                    & alp1(index)*(xh(i)-radius(index-1))/(radius(index)-radius(index-1))


         eps(i) = eps_r(index-1)*(xh(i)-radius(index))/(radius(index-1)-radius(index)) + & 
                    & eps_r(index)*(xh(i)-radius(index-1))/(radius(index)-radius(index-1))

         grr(i) = g_rr1(index-1)*(xh(i)-radius(index))/(radius(index-1)-radius(index)) + & 
                    & g_rr1(index)*(xh(i)-radius(index-1))/(radius(index)-radius(index-1))


         gtt(i) = grr(i)*xh(i)**2

         gpp(i) = grr(i)*xh(i)**2

         det = grr(i)*gtt(i)*gpp(i)

         sqdetg(i)= sqrt(det)

         urr(i) = 1.d0/grr(i)
         utt(i) = 1.d0/gtt(i)
         upp(i) = 1.d0/gpp(i) 

         krr(i) = 0.d0
         ktt(i) = 0.d0
         kpp(i) = 0.d0

         vr(i) = 0.d0
         velr(i) = 0.d0

!      grr(i) = (1.d0 + 1.d0/(2.d0*xh(i)))**4


      small_a(i) = 1.d0

      small_b(i) = 1.d0

      betaux(i) = 0.d0

      A_a(i) = 0.d0

!      psi(i) = (1.d0 + 1.d0/(2.d0*xh(i)))

!      alp(i) = psi(i)**(-2)

      xsi(i) = log(grr(i)**0.25)
      xsi2(i) = exp(-2.d0*xsi(i))

      trK(i) = 0.d0

      delta_x(i) = 0.d0




!write(*,*)'introduce a perturbation in TOV, 0.0%'


       if(rho(i).gt.rho_atm) then 

!          rho(i) = rho(i) + (rho(i)*1.)/100.d0     
!          p(i) = p(i) + (p(i)*4.5)/100.d0     
!          eps(i)=p(i)/((gamma-1.d0)*rho(i)) 

       end if

 
      rho_matter(i) = rho(i)*(1. + eps(i))


!!$      S_a(i) = (1.d0/(xh(i)**4*exp(12.d0*xsi(i))))*p(i)*(g_rr(i)**2)
!!$
!!$      S_b(i) = (1.d0/(xh(i)**4*exp(12.d0*xsi(i))))*p(i)*(g_tt(i)**2)

!!$      S_a(i) = (1.d0/det)*p(i)*(g_rr(i)**2)
!!$
!!$      S_b(i) = (1.d0/det)*p(i)*(g_tt(i)**2)


      S_a(i) = p(i)

      S_b(i) = p(i)


      j_r(i) = 0.d0
!      SS(i) = (1.d0/(xh(i)**4*exp(12.d0*xsi(i))))*p(i)*(g_rr(i)**2+2.d0*g_tt(i)**2)

      SS(i) = 3.d0*p(i)

      lorentz(i) = 1.d0

      write(301,*) xh(i),S_a(i),3.d0*p(i),rho_matter(i)

    h(i) = 1.d0+eps(i)+p(i)/rho(i)

    dens(i) = rho(i)*sqdetg(i)

    tau(i) = sqdetg(i)*(rho(i)*h(i) -p(i))-dens(i)

    srval(i) = 0.d0
    


   ENDDO

    DO i=imin,imax

      write(302,*) xh(i),xsi(i),alp(i),p(i),rho_matter(i)

   ENDDO


    close (unit=302)

   do i = 1,imax

      if(rho(i).gt.1.1d0*rho_atm) then 
              atmpmask(i) = 1
 !             write(*,*)'i',i
            else 
              atmpmask(i) = 0
 !             write(*,*)'i0',i
            end if
            
    end do


      capB = 0.d0 

  call Boundary_zaxis_func(xsi2) 
  call Boundary_zaxis_func(xsi) 
  call Boundary_zaxis_func(small_a) 
  call Boundary_zaxis_func(small_b) 
  call Boundary_zaxis_func(trK) 
  call Boundary_zaxis_func(A_a) 
  call Boundary_zaxis_asymfunc(delta_x) 
  call Boundary_zaxis_func(alp)
  call Boundary_zaxis_asymfunc(betaux) 


        call Boundary_zaxis_func(rho) 
        call Boundary_zaxis_func(eps) 
        call Boundary_zaxis_func(p) 
        call Boundary_zaxis_func(lorentz) 
        call Boundary_zaxis_asymfunc(velr) 
        call Boundary_zaxis_asymfunc(vr) 

  write(*,*) 'INITIAL LAPSE------>>>>',alp(1)





END SUBROUTINE TOV_source

















