  SUBROUTINE Runge_kutta_4th_BCs

  USE tmp_mdparam
  USE tmp_mdgrid
  USE tmp_mdmetric
  USE MD_Boundary

  IMPLICIT NONE


  REAL(KIND=double), DIMENSION(1:5,0:4) :: alphark4, betark4
  REAL(KIND=double), DIMENSION(:,:,:), ALLOCATABLE :: Urk4
  REAL(KIND=double), DIMENSION(:,:), ALLOCATABLE :: RHS,RHS3
  REAL(KIND=double), DIMENSION(1:8) :: urk_o


  INTEGER :: rkcount
  integer:: i,k,q,var,t
  integer :: coun1,coun2,coun_rate
  REAL(KIND=double) ::timelapse



  ALLOCATE ( Urk4(0:5,1:8,1-ghzx:nx+ghzx))
  ALLOCATE ( RHS(1:8,1-ghzx:nx+ghzx))
  ALLOCATE ( RHS3(1:8,1-ghzx:nx+ghzx))


    Urk4 = 0.d0

!----------------------------------------------------------------------------
! coefficients for RK4:
!----------------------------------------------------------------------------
    alphark4 = 0.d0

    alphark4(1,0)= 1.d0
    alphark4(2,0)= 0.44437049406734; alphark4(2,1)=0.55562950593266
    alphark4(3,0)= 0.62010185138540; alphark4(3,2)=0.37989814861460
    alphark4(4,0)= 0.17807995410773; alphark4(4,3)=0.82192004589227
    alphark4(5,0)= 0.00683325884039; alphark4(5,2)=0.51723167208978
    alphark4(5,3)= 0.12759831133288; alphark4(5,4)=0.34833675773694
    
    betark4 =0.d0
    
    betark4(1,0)= 0.39175222700392
    betark4(2,1)= 0.36841059262959
    betark4(3,2)= 0.25189177424738
    betark4(4,3)= 0.54497475021237
    betark4(5,3)= 0.08460416338212; betark4(5,4)=0.22600748319395

!----------------------------------------------------------------------------
! data at previous time step used for BCs:
!----------------------------------------------------------------------------

 !$OMP PARALLEL 
 !$OMP DO PRIVATE(i)

    do i=imin-ghzx,imax+ghzx

       Urk4(0,1,i) = xsi(i)
       Urk4(0,2,i) = small_a(i)
       Urk4(0,3,i) = small_b(i)
       Urk4(0,4,i) = trK(i)
       Urk4(0,5,i) = A_a(i)
       Urk4(0,6,i) = delta_x(i)
       Urk4(0,7,i) = alp(i)
       Urk4(0,8,i) = betaux(i)

    end do



 !$OMP END DO
 !$OMP END PARALLEL


     urk_o(1) = 0.d0 !xsi_o(:,:)
     urk_o(2) = 1.d0 !small_a_o(:,:)
     urk_o(3) = 1.d0 !small_b_o(:,:)
     urk_o(4) = 0.d0 !trK
     urk_o(5) = 0.d0
     urk_o(6) = 0.d0
     urk_o(7) = 1.d0 !alp
     urk_o(8) = 0.d0 !betaux


!----------------------------------------------------------------------------
! RK4 loop:
!----------------------------------------------------------------------------
rkcount=1


   do rkcount = 1,5


!************************************************************
!                  SPACETIME RHS                            !
!************************************************************

   call Einstein(RHS)


!************************************************************
!                      RK4 update:
!************************************************************


  IF (rkcount.eq.1) THEN

     call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)

     !$OMP PARALLEL 
     !$OMP DO PRIVATE(i,q)


        do i=imin,imax+3
           do q=1,8
              Urk4(1,q,i) = alphark4(1,0)*Urk4(0,q,i)+dt*betark4(1,0)*RHS(q,i)
           end do      


        xsi(i)     = Urk4(rkcount,1,i)
        small_a(i) = Urk4(rkcount,2,i)
        small_b(i) = Urk4(rkcount,3,i)
        trK(i)     = Urk4(rkcount,4,i)
        A_a(i)     = Urk4(rkcount,5,i)
        delta_x(i) = Urk4(rkcount,6,i)
        alp(i)     = Urk4(rkcount,7,i)
        betaux(i)  = Urk4(rkcount,8,i)


      end do
     !$OMP END DO 
     !$OMP END PARALLEL 
      
  else if (rkcount.eq.2) then 

     call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)


     !$OMP PARALLEL 
     !$OMP DO PRIVATE(i,q)

        do i=imin,imax+3

           do q=1,8
              Urk4(2,q,i) = alphark4(2,0)*Urk4(0,q,i)+alphark4(2,1)*Urk4(1,q,i)+dt*betark4(2,1)*RHS(q,i)
           end do      

        xsi(i)     = Urk4(rkcount,1,i)
        small_a(i) = Urk4(rkcount,2,i)
        small_b(i) = Urk4(rkcount,3,i)
        trK(i)     = Urk4(rkcount,4,i)
        A_a(i)     = Urk4(rkcount,5,i)
        delta_x(i) = Urk4(rkcount,6,i)
        alp(i)     = Urk4(rkcount,7,i)
        betaux(i)  = Urk4(rkcount,8,i)

        end do

     !$OMP END DO 
     !$OMP END PARALLEL 

      
  else if (rkcount.eq.3) then       

     call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)

     !$OMP PARALLEL 
     !$OMP DO PRIVATE(i,q)

        do i=imin,imax+3
           do q=1,8
              Urk4(3,q,i) = alphark4(3,0)*Urk4(0,q,i)+alphark4(3,2)*Urk4(2,q,i)+dt*betark4(3,2)*RHS(q,i)
           end do      

        xsi(i)     = Urk4(rkcount,1,i)
        small_a(i) = Urk4(rkcount,2,i)
        small_b(i) = Urk4(rkcount,3,i)
        trK(i)     = Urk4(rkcount,4,i)
        A_a(i)     = Urk4(rkcount,5,i)
        delta_x(i) = Urk4(rkcount,6,i)
        alp(i)     = Urk4(rkcount,7,i)
        betaux(i)  = Urk4(rkcount,8,i)


     end do
     !$OMP END DO 
     !$OMP END PARALLEL 
        
      
  else if (rkcount.eq.4) then       

     call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)

     !$OMP PARALLEL 
     !$OMP DO PRIVATE(i,q)

        do i=imin,imax+3
           do q=1,8
              Urk4(4,q,i) = alphark4(4,0)*Urk4(0,q,i)+alphark4(4,3)*Urk4(3,q,i)+dt*betark4(4,3)*RHS(q,i)
           end do      

        xsi(i)     = Urk4(rkcount,1,i)
        small_a(i) = Urk4(rkcount,2,i)
        small_b(i) = Urk4(rkcount,3,i)
        trK(i)     = Urk4(rkcount,4,i)
        A_a(i)     = Urk4(rkcount,5,i)
        delta_x(i) = Urk4(rkcount,6,i)
        alp(i)     = Urk4(rkcount,7,i)
        betaux(i)  = Urk4(rkcount,8,i)



     end do
     !$OMP END DO 

     !$OMP DO PRIVATE(i,q)
        do i=imin-3,imax+3
           do q=1,8
              RHS3(q,i) = RHS(q,i)
           end do
        end do
     !$OMP END DO 
     !$OMP END PARALLEL 
    
      
   else if (rkcount.eq.5) then 

     call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)


     !$OMP PARALLEL 
     !$OMP DO PRIVATE(i,q)

        do i=imin,imax+3

           do q=1,8
              Urk4(5,q,i) = alphark4(5,0)*Urk4(0,q,i)+alphark4(5,2)*Urk4(2,q,i)+alphark4(5,3)*Urk4(3,q,i)+ & 
             & alphark4(5,4)*Urk4(4,q,i)+dt*betark4(5,3)*RHS3(q,i)+dt*betark4(5,4)*RHS(q,i)
           end do      

        xsi(i)     = Urk4(rkcount,1,i)
        small_a(i) = Urk4(rkcount,2,i)
        small_b(i) = Urk4(rkcount,3,i)
        trK(i)     = Urk4(rkcount,4,i)
        A_a(i)     = Urk4(rkcount,5,i)
        delta_x(i) = Urk4(rkcount,6,i)
        alp(i)     = Urk4(rkcount,7,i)
        betaux(i)  = Urk4(rkcount,8,i)


     end do
     !$OMP END DO 
     !$OMP END PARALLEL 

            
  END IF


  call Boundary_zaxis_func(xsi) 
  call Boundary_zaxis_func(small_a) 
  call Boundary_zaxis_func(small_b) 
  call Boundary_zaxis_func(trK) 
  call Boundary_zaxis_func(A_a) 
  call Boundary_zaxis_asymfunc(delta_x) 
  call Boundary_zaxis_func(alp)
  call Boundary_zaxis_asymfunc(betaux) 
 

  
  call Derivatives
  


 end do
  
  DEALLOCATE (Urk4,RHS,RHS3)


END SUBROUTINE Runge_kutta_4th_BCs
    









