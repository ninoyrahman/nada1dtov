  SUBROUTINE Fluxespm(denstmp,srvaltmp,tautmp,vel,presstmp,alptmp,beta, & 
                     & sqdetgtmp,f1,f2,f5)
  USE tmp_mdparam
  implicit none
  
  REAL(KIND=double),INTENT(IN) :: denstmp,srvaltmp,tautmp,vel,sqdetgtmp
  REAL(KIND=double),INTENT(IN) :: alptmp,beta,presstmp 
  REAL(KIND=double),INTENT(OUT) :: f1,f2,f5

  REAL(KIND=double):: velnobeta

  velnobeta=vel-beta/alptmp

            f1 = denstmp*velnobeta
            f2 = srvaltmp*velnobeta + sqdetgtmp*presstmp
            f5 = tautmp*velnobeta + sqdetgtmp*presstmp*vel



   END SUBROUTINE Fluxespm

