
  subroutine Con2Prim(xsitmp,urrtmp,utttmp,upptmp,denstmp,srvaltmp,tautmp, & 
                     & grrtmp,gtttmp,gpptmp,alptmp,betauxtmp,sqdetgtmp, & 
                     & rhotmp,velrtmp,vrtmp,epstmp,lorentztmp,htmp,ptmp)
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdhydro
  USE tmp_mdmetric
  implicit none

  integer i,recovery_iteration

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: xsitmp,alptmp,sqdetgtmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: urrtmp,utttmp,upptmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: grrtmp,gtttmp,gpptmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: betauxtmp

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(INOUT) :: denstmp,srvaltmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(INOUT) :: tautmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(INOUT) :: rhotmp,velrtmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(INOUT) :: vrtmp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(INOUT) :: epstmp,lorentztmp,htmp,ptmp

  REAL(KIND=double) :: tolgvf
  REAL(KIND=double):: udenst,usrt,utaut
  REAL(KIND=double):: WW,eps_atm,press_atm

  REAL(KIND=double):: pguess,hw,vrval,drhobydpress,depsbydpress,temp1 
  REAL(KIND=double):: velrval,vv,vv1

  REAL(KIND=double):: gamma_eos_tmp,S2,taupresdens,dpdeps,dpdrho
  REAL(KIND=double):: c_sound_squared,p_error,f_p,df_dp,cs


  tolgvf = 1.d-8

 
   press_atm= kp*rho_atm**(gamma)
   eps_atm   = press_atm/((gamma-1.d0)*rho_atm)


!!$OMP PARALLEL
!!$OMP DO PRIVATE(i,p_error,pguess,recovery_iteration,hw,velxval,velyval,velzval, &
!!$OMP vxval,vyval,vzval,vv,WW,c_sound_squared,gamma_eos_tmp,f_p,df_dp, udenst,usxt,usyt,uszt,utaut,taupresdens,S2,cs, &
!!$OMP  dpdeps,dpdrho)

!write(*,*)'inicio', denstmp(31),srvaltmp(31), tautmp(31)

   DO i=imin,imax

   p_error = 1.0d0

! initial guess for the pressure:
   pguess = p(i)

   recovery_iteration = 1 
! start Newton-Raphson:

! do recovery_iteration=1,100 

   udenst = denstmp(i)/sqdetgtmp(i)
   usrt = srvaltmp(i)/sqdetgtmp(i)
   utaut = tautmp(i)/sqdetgtmp(i)

   S2 = usrt*usrt*urrtmp(i)

   do while (p_error.gt.tolgvf.and.recovery_iteration.le.100)

            recovery_iteration = recovery_iteration + 1



   taupresdens = utaut+pguess+udenst

!    write(*,*)'con2prim',taupresdens,S2

   rhotmp(i) = (udenst*sqrt(taupresdens*taupresdens-S2))/taupresdens

!   write(*,*)'con2prim',rhotmp(i),udenst,sqdetgtmp(i)
!   stop

   WW =  taupresdens/(sqrt(taupresdens*taupresdens-S2))

   epstmp(i) = (sqrt((taupresdens)**2 - S2)-pguess*WW-udenst)/udenst


! here introduce if-then for EOS:
!------------------------------------------------------------------------
   
!   if (atmpmask(i).eq.1) then
 
   if (rhotmp(i).gt.1.1*rho_atm) then

      IF (eos_type.eq.0) then
 
          ptmp(i)   = kp*rhotmp(i)**(gamma)
          dpdeps = 0.d0
          dpdrho =kp*gamma*rhotmp(i)**(gamma-1.d0)

          htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)
 
          cs = sqrt((dpdrho+dpdeps*ptmp(i)/rhotmp(i)**2)/htmp(i))
          c_sound_squared = cs*cs

      else if  (eos_type.eq.1) THEN
         
          ptmp(i) = (gamma-1.d0)*rhotmp(i)*epstmp(i)         
          
          dpdeps = (gamma - 1.0d0) * rhotmp(i)
          dpdrho = (gamma - 1.0d0) * epstmp(i)

          htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)
        
          cs = sqrt((dpdrho+dpdeps*ptmp(i)/rhotmp(i)**2)/htmp(i))
          c_sound_squared = cs*cs
          
      end IF



!   else  if (atmpmask(i).eq.0) then

      else 

      p_error = tolgvf
          c_sound_squared = 0.d0

   end if


!------------------------------------------------------------------------



   vrval = usrt/(rhotmp(i)*htmp(i)*WW*WW)

   velrval = urrtmp(i)*vrval 

   vv =  vrval*velrval 


   if (vv.ge.1.0d0) then

    write(*,*)'Nonfatal problem in:1',i,xh(i)

      vv = 0.0d0
      WW = 1.d0

      rhotmp(i)= rho_atm 
      ptmp(i) = press_atm
      lorentztmp(i) = 1.d0
      epstmp(i) = eps_atm
      htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)

      velrtmp(i) = 0.d0
            
      p_error = tolgvf
!      stop

   end if


   
   if (c_sound_squared.lt.0.0d0) then
     print*
     print*, 'Nonfatal problem in:'
     print*, 'Equation of state calculation'
     print*
     print*, 'Diagnosis:'
     print*, 'c_s^2 < 0'
     print*
     print*, 'Action:'
     print*, 'Assigning this cell to atmosphere'
!!$     print*
     
     atmpmask(i) = 0
!    write(*,*)'Nonfatal problem in:2' 
     
     vv = 0.0d0
     WW = 1.d0
     
     rhotmp(i)= rho_atm 
     ptmp(i) = press_atm
     lorentztmp(i) = 1.d0
     epstmp(i) = eps_atm
     htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)
     
     velrtmp(i) = 0.d0

     c_sound_squared = 0.d0
     
     p_error = tolgvf
     
   endif
   
!   if (atmpmask(i).eq.1) then 
   if (rhotmp(i).gt.1.1*rho_atm) then   

   f_p = ptmp(i) - pguess
   df_dp = vv*c_sound_squared - 1.0d0


!!$   f_p =  pguess - ptmp(i) 
!!$   temp1 = (utaut+udenst+pguess)**2 - S2
!!$   drhobydpress = udenst * S2 / (sqrt(temp1)*(udenst+utaut+pguess)**2)
!!$   depsbydpress = pguess * S2 / (rhotmp(i) * (udenst + utaut + pguess) * temp1)
!!$   df_dp = 1.0d0 - dpdrho*drhobydpress - &
!!$          & dpdeps*depsbydpress
!!$

   
   ptmp(i) = pguess - f_p / df_dp
   
   p_error = dabs(ptmp(i)/pguess - 1.0d0)

   pguess = ptmp(i)

   end if
  
 
 end do

      if (recovery_iteration.ge.150) then
         print* 
         print*, 'Fatal problem in:'
         print*, 'con2prim'
         print*
         print*, 'Diagnosis:'
         print*, '   No convergence reached in recovery calculation'
         print*
         print*, 'Action:'
         print*, '   Stopping',p_error,i,xh(i),rhotmp(i)
         print*,  pguess,ptmp(i)

         stop
      endif


   taupresdens = utaut+ptmp(i)+udenst
 
   rhotmp(i) = udenst*sqrt(taupresdens*taupresdens-S2)/taupresdens

   WW =  taupresdens/(sqrt(taupresdens*taupresdens-S2))

   epstmp(i)=ptmp(i)/((gamma-1.d0)*rhotmp(i))

!   epstmp(i) = (sqrt( (taupresdens)**2 - S2)-ptmp(i)*WW-udenst)/udenst
   
   htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)


!!$   htmp(i) = taupresdens/(WW*WW*rhotmp(i))
!!$
!!$   epstmp(i) = htmp(i) - 1.d0 - pguess/rhotmp(i)

   vrtmp(i) = usrt/(rhotmp(i)*htmp(i)*WW*WW)

   velrtmp(i) = urrtmp(i)*vrtmp(i) 

   vv =  vrtmp(i)*velrtmp(i) 

   lorentztmp(i) = WW

   atmpmask(i) = 1

!!$

   if((rhotmp(i).gt.1.01*rho_atm).and.(epstmp(i).gt.1.d-20)) then


      else 



      rhotmp(i)= rho_atm 
      ptmp(i) = press_atm
      lorentztmp(i) = 1.d0
      epstmp(i) = eps_atm
      htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)

      velrtmp(i) = 0.d0
      vrtmp(i) = 0.d0
      
      denstmp(i)  = rhotmp(i)*sqdetgtmp(i)
      srvaltmp(i) = 0.d0
      tautmp(i)   = (rhotmp(i)*htmp(i)-ptmp(i))*sqdetgtmp(i)-denstmp(i)

      rho_matter(i) = 0.d0 
             
      S_a(i) = 0.d0 
      S_b(i) = 0.d0 
      
      j_r(i) = 0.d0 
      

!!$
!!$      if (i.lt.8) then 
!!$
!!$
!!$      rhotmp(i)=rhotmp(12) 
!!$      ptmp(i) = ptmp(12) 
!!$      lorentztmp(i) = lorentztmp(12)
!!$      epstmp(i) = epstmp(12)
!!$      htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)
!!$
!!$
!!$      velrtmp(i) = velrtmp(12)
!!$      vrtmp(i) = vrtmp(12) 
!!$      
!!$      denstmp(i)  = denstmp(12)
!!$      srvaltmp(i) = srvaltmp(12) 
!!$      tautmp(i)   = tautmp(12)
!!$
!!$
!!$
!!$!!      atmpmask(i) = 0
!!$
!!$!      write(*,*)'Nonfatal problem in:3',i
!!$
!!$     end if

      
    end if
!!$      




  if (I_EX.gt.0) then 

    if (i.lt.I_EX) then 
!    if (i.lt.7) then 

      rhotmp(i)= rho_atm 
      ptmp(i) = press_atm
      lorentztmp(i) = 1.d0
      epstmp(i) = eps_atm
      htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)

      velrtmp(i) = 0.d0
      vrtmp(i) = 0.d0
      
      denstmp(i)  = rhotmp(i)*sqdetgtmp(i)
      srvaltmp(i) = 0.d0
      tautmp(i)   = (rhotmp(i)*htmp(i)-ptmp(i))*sqdetgtmp(i)-denstmp(i)

      rho_matter(i) = 0.d0 
             
      S_a(i) = 0.d0 
      S_b(i) = 0.d0 
      
      j_r(i) = 0.d0 
      

!      else if (i.eq.7) then 
      else if (i.eq.I_EX) then 

      rhotmp(i)=rhotmp(i+1) 
      ptmp(i) = ptmp(i+1) 
      lorentztmp(i) = lorentztmp(i+1)
      epstmp(i) = epstmp(i+1)
      htmp(i) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)


      velrtmp(i) = velrtmp(i+1)
      vrtmp(i) = vrtmp(i+1) 
      
      denstmp(i)  = denstmp(i+1)
      srvaltmp(i) = srvaltmp(i+1) 
      tautmp(i)   = tautmp(i+1)


      rhotmp(i-1)=rhotmp(i+1) 
      ptmp(i-1) = ptmp(i+1) 
      lorentztmp(i-1) = lorentztmp(i+1)
      epstmp(i-1) = epstmp(i+1)
      htmp(i-1) = 1.d0 + epstmp(i) + ptmp(i)/rhotmp(i)


      velrtmp(i-1) = velrtmp(i+1)
      vrtmp(i-1) = vrtmp(i+1) 
      
      denstmp(i-1)  = denstmp(i+1)
      srvaltmp(i-1) = srvaltmp(i+1) 
      tautmp(i-1)   = tautmp(i+1)


     END if


    END if



  END DO


!!$OMP END DO
!!$OMP END PARALLEL

 
 END subroutine Con2Prim
