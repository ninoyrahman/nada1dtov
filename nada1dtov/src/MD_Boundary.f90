
  MODULE MD_Boundary
  IMPLICIT NONE
  
  CONTAINS

    
!----------------------------------------------------------------------------------------------    

    SUBROUTINE Boundary_zaxis_func(func)
      USE tmp_mdparam
      USE tmp_mdgrid
      
      IMPLICIT NONE
      
      REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(INOUT) :: func
      REAL(KIND=double), PARAMETER :: onethird = 1.d0/3.d0
      INTEGER :: i,k,q 
      
      if(gridtype.eq.1) then 
         
          func(0)  = func(1)
          func(-1) = func(2)
          func(-2) = func(3)

      else if(gridtype.eq.0) then    

          func(0) = (4.d0*func(1)-func(2))*onethird
          func(-1) = func(1)
          func(-2) = func(2)
        
      end if
      
    END SUBROUTINE Boundary_zaxis_func


    SUBROUTINE Boundary_zaxis_asymfunc(func)
      USE tmp_mdparam
      USE tmp_mdgrid
      
      IMPLICIT NONE
      
      REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(INOUT) :: func
      REAL(KIND=double), PARAMETER :: onethird = 1.d0/3.d0
      INTEGER :: i,k,q 
      
      if(gridtype.eq.1) then 
         
          func(0)  = -func(1)
          func(-1) = -func(2)
          func(-2) = -func(3)

      else if(gridtype.eq.0) then    

          func(0) = 0.d0
          func(-1) = -func(1)
          func(-2) = -func(2)
        
      end if
      
    END SUBROUTINE Boundary_zaxis_asymfunc

        

! Radiative boundaries with RK4:


  SUBROUTINE RHSBCs(func,func_o,RHSBCtmp)
  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric

  IMPLICIT NONE

  REAL(KIND=double), DIMENSION(1:12,-ghzx+1:nx+ghzx),INTENT(IN) :: func
  REAL(KIND=double), DIMENSION(1:9),INTENT(IN) :: func_o
  REAL(KIND=double), DIMENSION(1:12,-ghzx+1:nx+ghzx),INTENT(INOUT) :: RHSBCtmp

  INTEGER :: i,k,t,q 
  REAL(KIND=double), DIMENSION(1:9)::velcc

        
        velcc=1.d0
        velcc(7)=sqrt(2.d0)
        velcc(6)=sqrt(2.d0)



        do q=1,9


          RHSBCtmp(q,imax+1) = ((func_o(q)-func(q,imax+1))/xh(imax+1) - & 
                           & ((func(q,imax+2)-func(q,imax))/(2.d0*dx)))*velcc(q) 

          RHSBCtmp(q,imax+2) = ((func_o(q)-func(q,imax+2))/xh(imax+2) - & 
                           & ((func(q,imax+3)-func(q,imax+1))/(2.d0*dx)))*velcc(q)
!!$
!!$
!!$
!!$
!!$          RHSBCtmp(q,imax+1) = ((func_o(q)-func(q,imax+1))/xh(imax+1) - & 
!!$                           & ((func(q,imax-1)-4.d0*func(q,imax)+ & 
!!$                           & 3.d0*func(q,imax+1))/(2.d0*dx)))*velcc(q)
!!$
!!$
!!$
!!$          RHSBCtmp(q,imax+2) = ((func_o(q)-func(q,imax+2))/xh(imax+2) - & 
!!$                           & ((func(q,imax)-4.d0*func(q,imax+1)+ & 
!!$                           & 3.d0*func(q,imax+2))/(2.d0*dx)))*velcc(q)
!!$
!!$



          RHSBCtmp(q,imax+3) = ((func_o(q)-func(q,imax+3))/xh(imax+3) - & 
                           & ((func(q,imax+1)-4.d0*func(q,imax+2)+ & 
                           & 3.d0*func(q,imax+3))/(2.d0*dx)))*velcc(q)





        end do



!!$
!!$          RHSBCtmp(6,imax+1) = ((func_o(6)-func(6,imax+1))/(xh(imax+1)**2) - & 
!!$                           & ((func(6,imax+2)-func(6,imax))/(2.d0*dx)))*velcc(6) 
!!$
!!$          RHSBCtmp(6,imax+2) = ((func_o(6)-func(6,imax+2))/(xh(imax+2)**2) - & 
!!$                           & ((func(6,imax+3)-func(6,imax+1))/(2.d0*dx)))*velcc(6)
!!$
!!$
!!$          RHSBCtmp(6,imax+3) = ((func_o(6)-func(6,imax+3))/(xh(imax+3)**2) - & 
!!$                           & ((func(6,imax+1)-4.d0*func(6,imax+2)+ & 
!!$                           & 3.d0*func(6,imax+3))/(2.d0*dx)))*velcc(6)
!!$
!!$
!!$
!!$          RHSBCtmp(8,imax+1) = ((func_o(8)-func(8,imax+1))/(xh(imax+1)**2) - & 
!!$                           & ((func(8,imax+2)-func(8,imax))/(2.d0*dx)))*velcc(8) 
!!$
!!$          RHSBCtmp(8,imax+2) = ((func_o(8)-func(8,imax+2))/(xh(imax+2)**2) - & 
!!$                           & ((func(8,imax+3)-func(8,imax+1))/(2.d0*dx)))*velcc(8)
!!$
!!$
!!$          RHSBCtmp(8,imax+3) = ((func_o(8)-func(8,imax+3))/(xh(imax+3)**2) - & 
!!$                           & ((func(8,imax+1)-4.d0*func(8,imax+2)+ & 
!!$                           & 3.d0*func(8,imax+3))/(2.d0*dx)))*velcc(8)
!!$

       

      END SUBROUTINE RHSBCs


END MODULE MD_Boundary



