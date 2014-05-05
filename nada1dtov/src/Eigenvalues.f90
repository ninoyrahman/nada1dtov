
  subroutine Eigenvalues(rhotmp,velrtmp,epstmp,ptmp,vvtmp, & 
       & utmp,alptmp,betax,l1,l2,l3,l4,l5)

  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdhydro
  implicit none


  REAL(KIND=double),INTENT(IN) :: rhotmp,velrtmp,epstmp,ptmp,vvtmp
  REAL(KIND=double),INTENT(IN) :: utmp,alptmp
  REAL(KIND=double),INTENT(IN) :: betax
  REAL(KIND=double),INTENT(OUT):: l1, l2, l3, l4, l5
  
  integer i
  
  REAL(KIND=double)::  htmp

  REAL(KIND=double)::  cs,cs2
  REAL(KIND=double)::  dpdeps,kappa
  REAL(KIND=double)::  lbp,lbm,dpdrho
  REAL(KIND=double):: Vxcap,Vycap,Vzcap,Vuxcap,Vuycap,Vuzcap




       htmp = 1 + epstmp + ptmp/rhotmp  

       if (eos_type.eq.0) then 

          dpdeps = 0.d0
          dpdrho =kp*gamma*rhotmp**(gamma-1.d0)
          
       else if (eos_type.eq.1) then 
          
          dpdeps = (gamma - 1.0d0) * rhotmp
          dpdrho = (gamma - 1.0d0) * epstmp
          
       end if


        cs = sqrt((dpdrho+dpdeps*ptmp/rhotmp**2)/htmp)
        cs2 = cs*cs
        

!! eigenvalues: following M.Shibata's ordering:
!!
!!     l1 -> lambda +
!!     l2 -> lambda 0
!!     l3 -> lambda 0
!!     l4 -> lambda 0
!!     l5 -> lambda -

!!$        l1 = alptmp/(1.d0-vvtmp*cs2)
!!$        l1 = l1*(velrtmp*(1.d0-cs2)+cs*sqrt((1.d0-vvtmp)*(utmp*(1.d0-vvtmp*cs2)-velrtmp*velrtmp* & 
!!$           &  (1.d0-cs2) ) ) )-betax
!!$
!!$        l2 = alptmp*velrtmp-betax
!!$        l3 = alptmp*velrtmp-betax
!!$        l4 = alptmp*velrtmp-betax
!!$
!!$
!!$        l5 = alptmp/(1.d0-vvtmp*cs2)
!!$        l5 = l5*(velrtmp*(1.d0-cs2)-cs*sqrt((1.d0-vvtmp)*(utmp*(1.d0-vvtmp*cs2)-velrtmp*velrtmp* & 
!!$           &  (1.d0-cs2) ) ) )-betax



        l1 = alptmp/(1.d0-vvtmp*cs2)
        l1 = l1*(velrtmp*(1.d0-cs2)+cs*sqrt((1.d0-vvtmp)*(utmp*(1.d0-vvtmp*cs2)-velrtmp*velrtmp* & 
           &  (1.d0-cs2) ) ) )-betax

        l2 = alptmp*velrtmp-betax
        l3 = alptmp*velrtmp-betax
        l4 = alptmp*velrtmp-betax


        l5 = alptmp/(1.d0-vvtmp*cs2)
        l5 = l5*(velrtmp*(1.d0-cs2)-cs*sqrt((1.d0-vvtmp)*(utmp*(1.d0-vvtmp*cs2)-velrtmp*velrtmp* & 
           &  (1.d0-cs2) ) ) )-betax


!!$        if (abs(l1).gt.0d0) then 
!!$           else 
!!$           l1 = 0.d0
!!$        end if
!!$
!!$
!!$        if (abs(l2).gt.0d0) then 
!!$           else 
!!$           l2 = 0.d0
!!$        end if
!!$
!!$        if (abs(l3).gt.0d0) then 
!!$           else 
!!$           l3 = 0.d0
!!$        end if
!!$
!!$        if (abs(l4).gt.0d0) then 
!!$           else 
!!$           l4 = 0.d0
!!$        end if
!!$
!!$        if (abs(l5).gt.0d0) then 
!!$           else 
!!$           l5 = 0.d0
!!$        end if
!!$


        end subroutine Eigenvalues
        



