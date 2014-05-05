SUBROUTINE Runge_kutta_2nd_BCs

    USE tmp_mdparam
    USE tmp_mdgrid
    USE tmp_mdmetric
    USE tmp_mdhydro
    USE MD_HydroSubroutines
    USE MD_Boundary
    USE PPM_RECONSTRUCTION

    IMPLICIT NONE
 

    REAL(KIND=double), DIMENSION(:,:,:), ALLOCATABLE :: Urk4
    REAL(KIND=double), DIMENSION(:,:), ALLOCATABLE :: RHS,KO_disp1

    REAL(KIND=double), DIMENSION(1:9) :: urk_o


    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_umas,RHS2_umas
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_u1,RHS2_u1
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_un,RHS2_un


    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_umasT,RHS2_umasT
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_u1T,RHS2_u1T
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_unT,RHS2_unT


    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_umasD,RHS2_umasD
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_u1D,RHS2_u1D
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_unD,RHS2_unD

    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_umascapB,RHS2_umascapB
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_u1capB,RHS2_u1capB
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_uncapB,RHS2_uncapB


    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_umasbeta,RHS2_umasbeta
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_u1beta,RHS2_u1beta
    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3_unbeta,RHS2_unbeta

    REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: RHS3,RHS2

    REAL(KIND=double), DIMENSION(1:12,-ghzx+1:nx+ghzx) :: func_in


    REAL(KIND=double) ::  op1, op2, op3, op4, op5,alpi_right,alpk_right,alpi_left,alpk_left
    REAL(KIND=double), DIMENSION(1:nx) :: su1,su2,su5
    REAL(KIND=double), DIMENSION(0:nx+1) :: fna,fnr,fne




    INTEGER :: rkcount,hydro
    integer:: i,k,q,var,t
    integer :: coun1,coun2,coun_rate
    REAL(KIND=double) ::timelapse,Const



    ALLOCATE ( Urk4(0:2,1:12,1-ghzx:nx+ghzx))
    ALLOCATE ( RHS(1:12,1-ghzx:nx+ghzx))

    ALLOCATE ( KO_disp1(1:12,1-ghzx:nx+ghzx))


    Urk4 = 0.d0


    hydro = 1



    !----------------------------------------------------------------------------
    ! data at previous time step used for BCs:
    !----------------------------------------------------------------------------

    !!$OMP PARALLEL
    !!$OMP DO PRIVATE(i)

    do i=imin-ghzx,imax+ghzx

        !       Urk4(0,1,i) = xsi2(i)
        Urk4(0,1,i) = xsi(i)
        Urk4(0,2,i) = small_a(i)
        Urk4(0,3,i) = small_b(i)
        Urk4(0,4,i) = trK(i)
        Urk4(0,5,i) = A_a(i)
        Urk4(0,6,i) = delta_x(i)
        Urk4(0,7,i) = alp(i)
        Urk4(0,8,i) = betaux(i)

        Urk4(0,9,i) = capB(i)

        Urk4(0,10,i) = dens(i)
        Urk4(0,11,i) = srval(i)
        Urk4(0,12,i) = tau(i)
 


    end do



    !!$OMP END DO
    !!$OMP END PARALLEL


    !     urk_o(1) = 1.d0 !xsi2_o(:,:)

    urk_o(1) = 0.d0 !xsi_o(:,:)

    urk_o(2) = 1.d0 !small_a_o(:,:)
    urk_o(3) = 1.d0 !small_b_o(:,:)
    urk_o(4) = 0.d0 !trK
    urk_o(5) = 0.d0 !A_a
    urk_o(6) = 0.d0 !delta
    urk_o(7) = 1.d0 !alp
    urk_o(8) = 0.d0 !betaux

    urk_o(9) = 0.d0 !capB


    !----------------------------------------------------------------------------
    ! RK4 loop:
    !----------------------------------------------------------------------------
    rkcount=1


    !  eps= 0.00d0
       

   

    do rkcount = 1,2


        !************************************************************
        !                  SPACETIME RHS                            !
        !************************************************************


        !************************************************************
        !                      RK4 update:
        !************************************************************


        IF (rkcount.eq.1) THEN

            call Derivatives


            if (hydro.eq.1) then

                call matter_source(lorentz,rho,velr,eps,betaux,alp, &
                    & grr,urr)

            end if


            call Einstein(RHS)

            RHSdeltar(:) = RHS(6,:)

            !     call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)

            func_in(:,:) = Urk4(rkcount-1,:,:)

            call RHSBCs(func_in,urk_o,RHS)



            !************************************************************
            !************************************************************
            !                  HYDRO UPDATE                             !
            !************************************************************
            !************************************************************
            if (hydro.eq.1) then

                ! Reconstruction: subroutines in MD_HydroSubroutines.f90
!                call rlnlx(rho,rhom,rhop)
!                call rlnlx(velr,velrm,velrp)
!                call rlnlx(eps,epsm,epsp)
!                call rlnlx(p,pressm,pressp)

                call ppm(rho,rhom,rhop,velr,velrm,velrp,eps,epsm,epsp,p,pressm,pressp)

                !Riemann solver:
                call solverr_hlle(betaux,alp,grr,gtt,gpp, &
                    &  urr,utt,upp,xsi,fna,fnr,fne)


                ! compute source terms of hydro eqs.
                call source(lorentz,rho,velr,eps,betaux,alp, &
                    & krr,ktt,kpp,urr,utt,upp, &
                    & grr,gtt,gpp, &
                    & xsi,dxdalp,dxdbetaux, &
                    & dxdgrr,dxdgtt,dxdgpp, &
                    & su1,su2,su5)

  
                RHS(10,:)= 0.d0
                RHS(11,:)= 0.d0
                RHS(12,:)= 0.d0

                DO i=imin,imax

                    alpi_right = 0.5d0*(alp(i) + alp(i+1))
                    alpi_left  = 0.5d0*(alp(i) + alp(i-1))


                    op1 = (alpi_right*fna(i) - alpi_left*fna(i-1))*invdx


                    op2 = (alpi_right*fnr(i) - alpi_left*fnr(i-1))*invdx


                    op5 = (alpi_right*fne(i) - alpi_left*fne(i-1))*invdx


                    RHS(10,i)=-(op1)
                    RHS(11,i)=-(op2 - alp(i)*sqdetg(i)*su2(i))!RHSJx(i)
                    RHS(12,i)=-(op5 - alp(i)*sqdetg(i)*su5(i))!RHSE(i,k)



                END DO

     
            end if





            !!$OMP PARALLEL
            !!$OMP DO PRIVATE(i,q)


            do i=imin,imax+3

                do q=1,12

                    Urk4(1,q,i) = Urk4(0,q,i)+dt*RHS(q,i)

                end do

            end do

            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd(RHS2,RHS3)

            RHS3_un = RHS3
            RHS2_un = RHS2


            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_trK(RHS2,RHS3)

            RHS3_unT = RHS3
            RHS2_unT = RHS2


            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_deltax(RHS2,RHS3)

            RHS3_unD = RHS3
            RHS2_unD = RHS2


            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_capB(RHS2,RHS3)

            RHS3_uncapB = RHS3
            RHS2_uncapB = RHS2



            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_beta(RHS2,RHS3)

            RHS3_unbeta = RHS3
            RHS2_unbeta = RHS2


 

            do i=imin,imax+3


                !        xsi2(i)     = Urk4(rkcount,1,i)
                !        xsi(i) = -log(xsi2(i))/2.d0

                !!$        xsi(i)     = Urk4(rkcount,1,i)
                !!$
                !!$        small_a(i) = Urk4(rkcount,2,i)
                !!$        small_b(i) = Urk4(rkcount,3,i)
                !!$
                !!$        trK(i)     = Urk4(rkcount,4,i)
                !!$        A_a(i)     = Urk4(rkcount,5,i)
                !!$
                !!$        delta_x(i) = Urk4(rkcount,6,i)
                !!$
                !!$        alp(i)     = Urk4(rkcount,7,i)
                !!$
                !!$        betaux(i)  = Urk4(rkcount,8,i)
                !!$        capB(i)    = Urk4(rkcount,9,i)

                dens(i) = Urk4(rkcount,10,i)
                srval(i)= Urk4(rkcount,11,i)
                tau(i)  = Urk4(rkcount,12,i)


            end do
 
            call Boundary_zaxis_func(xsi2)
            call Boundary_zaxis_func(xsi)
            call Boundary_zaxis_func(small_a)
            call Boundary_zaxis_func(small_b)
            call Boundary_zaxis_func(trK)
            call Boundary_zaxis_func(A_a)
            call Boundary_zaxis_asymfunc(delta_x)
            call Boundary_zaxis_func(alp)
            call Boundary_zaxis_asymfunc(betaux)
            call Boundary_zaxis_asymfunc(capB)
 
            psi = exp(xsi)
 
            call Derivatives

 
            if (hydro.eq.1) then
     
                call Con2Prim(xsi,urr,utt,upp,dens,srval,tau, &
                    & grr,gtt,gpp,alp,betaux,sqdetg, &
                    & rho,velr,vr,eps,lorentz,h,p)


                call Boundary_zaxis_func(rho)
                call Boundary_zaxis_func(eps)
                call Boundary_zaxis_func(p)
                call Boundary_zaxis_func(lorentz)
                call Boundary_zaxis_asymfunc(velr)
                call Boundary_zaxis_asymfunc(vr)



                call matter_source(lorentz,rho,velr,eps,betaux,alp, &
                    & grr,urr)


            end if


            call Einstein_2nd(RHS2,RHS3)

            RHS2_u1 = RHS2
            RHS3_u1 = RHS3


            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_trK(RHS2,RHS3)


            RHS2_u1T = RHS2
            RHS3_u1T = RHS3


            !!$
            !!$      do i=imin,imax
            !!$
            !!$         !v^1 = v^n + Delta t * [0.5 * L_2 (u^n) + 0.5 * L_2 (u^1) + L_3 (u^n, v^n)- factor]
            !!$
            !!$         Urk4(1,5,i) =  Urk4(0,5,i)+dt*(0.5d0*(RHS2_un(i)+RHS2_u1(i)) + RHS3_un(i))
            !!$         A_a(i)     = Urk4(1,5,i)
            !!$
            !!$
            !!$        Urk4(1,4,i) = Urk4(0,4,i)+dt*(0.5d0*(RHS2_unT(i)+RHS2_u1T(i)) + RHS3_unT(i))
            !!$
            !!$         trK(i)     = Urk4(1,4,i)
            !!$
            !!$      end do
            !!$

            call Boundary_zaxis_func(A_a)
            call Boundary_zaxis_func(trK)

            call Derivatives
 
            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_deltax(RHS2,RHS3)


            RHS2_u1D = RHS2
            RHS3_u1D = RHS3

            !!$
            !!$      do i=imin,imax!+3
            !!$         !v^1 = v^n + Delta t * [0.5 * L_2 (u^n) + 0.5 * L_2 (u^1) + L_3 (u^n, v^n)- factor]
            !!$         Urk4(1,6,i) = Urk4(0,6,i)+dt*(0.5d0*(RHS2_unD(i)+RHS2_u1D(i)) + RHS3_unD(i))
            !!$         delta_x(i)     = Urk4(1,6,i)
            !!$      end do


            call Boundary_zaxis_asymfunc(delta_x)


            call Einstein(RHS)

            RHSdeltar(:) = RHS(6,:)

            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_capB(RHS2,RHS3)


            RHS2_u1capB = RHS2
            RHS3_u1capB = RHS3


            !!!$           RHS2 = 0.d0
            !!!$           RHS3 = 0.d0
            !!!$
            !!!$          call Einstein_2nd_beta(RHS2,RHS3)
            !!!$
            !!!$
            !!!$        RHS2_u1beta = RHS2
            !!!$        RHS3_u1beta = RHS3


            !!$      do i=imin,imax!+3
            !!$!v^1 = v^n + Delta t * [0.5 * L_2 (u^n) + 0.5 * L_2 (u^1) + L_3 (u^n, v^n)- factor]
            !!$         Urk4(1,9,i) = Urk4(0,9,i)+dt*(0.5d0*(RHS2_uncapB(i)+RHS2_u1capB(i)) + RHS3_uncapB(i))
            !!$         capB(i)     = Urk4(1,9,i)
            !!$      end do


            call Boundary_zaxis_asymfunc(capB)
            call Boundary_zaxis_asymfunc(betaux)

      
        else if (rkcount.eq.2) then

            call Derivatives

            if (hydro.eq.1) then

                call matter_source(lorentz,rho,velr,eps,betaux,alp, &
                    & grr,urr)

            end if


            call Einstein(RHS)

            RHSdeltar(:) = RHS(6,:)
  
            !   call RHSBCs(Urk4(rkcount-1,:,:),urk_o,RHS)

            func_in(:,:) = Urk4(rkcount-1,:,:)

            call RHSBCs(func_in,urk_o,RHS)



            !************************************************************
            !************************************************************
            !                  HYDRO UPDATE                             !
            !************************************************************
            !***********************************************************

            if (hydro.eq.1) then


!                call rlnlx(rho,rhom,rhop)
!                call rlnlx(velr,velrm,velrp)
!                call rlnlx(eps,epsm,epsp)
!                call rlnlx(p,pressm,pressp)
                call ppm(rho,rhom,rhop,velr,velrm,velrp,eps,epsm,epsp,p,pressm,pressp)


                call solverr_hlle(betaux,alp,grr,gtt,gpp, &
                    &  urr,utt,upp,xsi,fna,fnr,fne)


                call source(lorentz,rho,velr,eps,betaux,alp, &
                    & krr,ktt,kpp,urr,utt,upp, &
                    & grr,gtt,gpp, &
                    & xsi,dxdalp,dxdbetaux, &
                    & dxdgrr,dxdgtt,dxdgpp, &
                    & su1,su2,su5)

      
                RHS(10,:)= 0.d0
                RHS(11,:)= 0.d0
                RHS(12,:)= 0.d0


                DO i=imin,imax

                    !this is alp*xh:


                    alpi_right = 0.5d0*(alp(i) + alp(i+1))
                    alpi_left  = 0.5d0*(alp(i) + alp(i-1))


                    op1 = (alpi_right*fna(i) - alpi_left*fna(i-1))*invdx


                    op2 = (alpi_right*fnr(i) - alpi_left*fnr(i-1))*invdx


                    op5 = (alpi_right*fne(i) - alpi_left*fne(i-1))*invdx


                    RHS(10,i)=-(op1)
                    RHS(11,i)=-(op2 - alp(i)*sqdetg(i)*su2(i))
                    RHS(12,i)=-(op5 - alp(i)*sqdetg(i)*su5(i))



                END DO



            end if


            !!$OMP PARALLEL
            !!$OMP DO PRIVATE(i,q)

            do i=imin,imax+3

                do q=1,12


                    Urk4(2,q,i) = 0.5d0*(Urk4(0,q,i)+ Urk4(1,q,i)+dt*RHS(q,i))


                end do


                !        xsi2(i)     = Urk4(rkcount,1,i)
                !        xsi(i) = -log(xsi2(i))/2.d0

                !!$        xsi(i)     = Urk4(rkcount,1,i)
                !!$        small_a(i) = Urk4(rkcount,2,i)
                !!$        small_b(i) = Urk4(rkcount,3,i)
                !!$
                !!$        trK(i)     = Urk4(rkcount,4,i)
                !!$        A_a(i)     = Urk4(rkcount,5,i)
                !!$
                !!$        delta_x(i) = Urk4(rkcount,6,i)
                !!$
                !!$        alp(i)     = Urk4(rkcount,7,i)

                !!$        betaux(i)  = Urk4(rkcount,8,i)
                !!$        capB(i)    = Urk4(rkcount,9,i)

                dens(i) = Urk4(rkcount,10,i)
                srval(i)= Urk4(rkcount,11,i)
                tau(i)  = Urk4(rkcount,12,i)


            end do

            call Boundary_zaxis_func(xsi2)
            call Boundary_zaxis_func(xsi)
            call Boundary_zaxis_func(small_a)
            call Boundary_zaxis_func(small_b)
            call Boundary_zaxis_func(trK)
            call Boundary_zaxis_func(A_a)
            call Boundary_zaxis_asymfunc(delta_x)
            call Boundary_zaxis_func(alp)
            call Boundary_zaxis_asymfunc(betaux)
            call Boundary_zaxis_asymfunc(capB)
 

  
            call Derivatives

            psi = exp(xsi)
    
            if (hydro.eq.1) then

                call Con2Prim(xsi,urr,utt,upp,dens,srval,tau, &
                    & grr,gtt,gpp,alp,betaux,sqdetg, &
                    & rho,velr,vr,eps,lorentz,h,p)


                call Boundary_zaxis_func(rho)
                call Boundary_zaxis_func(eps)
                call Boundary_zaxis_func(p)
                call Boundary_zaxis_func(lorentz)
                call Boundary_zaxis_asymfunc(velr)
                call Boundary_zaxis_asymfunc(vr)


                call matter_source(lorentz,rho,velr,eps,betaux,alp, &
                    & grr,urr)

            end if
     

            call Einstein_2nd(RHS2,RHS3)

            RHS3_umas = RHS3
            RHS2_umas = RHS2


            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_trK(RHS2,RHS3)


            RHS3_umasT = RHS3
            RHS2_umasT = RHS2


            !!$      do i=imin,imax!+3
            !!$
            !!$         !v^(n+1) = v^n + 0.5 * Delta t * [L_2(u^n) + L_2(u^(n+1)+ L_3 (u^n,v^n) + L_3 (u^1, v^1)]
            !!$
            !!$         Urk4(2,5,i) = Urk4(0,5,i)+0.5*dt*(RHS2_un(i)+ RHS2_umas(i)+ RHS3_un(i)+RHS3_u1(i))
            !!$         A_a(i)     = Urk4(2,5,i)
            !!$
            !!$
            !!$         Urk4(2,4,i) = Urk4(0,4,i)+0.5*dt*(RHS2_unT(i)+ RHS2_umasT(i)+ RHS3_unT(i)+RHS3_u1T(i))
            !!$         trK(i)     = Urk4(2,4,i)
            !!$
            !!$      end do


            call Boundary_zaxis_func(A_a)
            call Boundary_zaxis_func(trK)


            call Derivatives

            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_deltax(RHS2,RHS3)


            RHS3_umasD = RHS3
            RHS2_umasD = RHS2


            !!$       do i=imin,imax!+3
            !!$
            !!$         !v^(n+1) = v^n + 0.5 * Delta t * [L_2(u^n) + L_2(u^(n+1)+ L_3 (u^n,v^n) + L_3 (u^1, v^1)]
            !!$
            !!$         Urk4(2,6,i) = Urk4(0,6,i)+0.5*dt*(RHS2_unD(i)+ RHS2_umasD(i)+ RHS3_unD(i)+RHS3_u1D(i))
            !!$         delta_x(i)     = Urk4(2,6,i)
            !!$
            !!$      end do


            call Boundary_zaxis_asymfunc(delta_x)

            call Einstein(RHS)
            RHSdeltar(:) = RHS(6,:)


            RHS2 = 0.d0
            RHS3 = 0.d0

            call Einstein_2nd_capB(RHS2,RHS3)


            RHS3_umascapB = RHS3
            RHS2_umascapB = RHS2



            !!!$           RHS2 = 0.d0
            !!!$           RHS3 = 0.d0
            !!!$
            !!!$           call Einstein_2nd_beta(RHS2,RHS3)
            !!!$
            !!!$
            !!!$           RHS3_umasbeta = RHS3
            !!!$           RHS2_umasbeta = RHS2

 
 
            !!$      do i=imin,imax!+3
            !!$         Urk4(2,9,i) = Urk4(0,9,i)+0.5*dt*(RHS2_uncapB(i)+ RHS2_umascapB(i)+ RHS3_uncapB(i)+RHS3_u1capB(i))
            !!$
            !!$         capB(i)     = Urk4(2,9,i)
            !!$
            !!$
            !!$      end do



            call Boundary_zaxis_asymfunc(capB)
            call Boundary_zaxis_asymfunc(betaux)



        END IF




    end do
  
    DEALLOCATE (Urk4,RHS)


END SUBROUTINE Runge_kutta_2nd_BCs
    









