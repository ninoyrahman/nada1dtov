Subroutine Shock_Tube
    USE tmp_mdparam
    USE tmp_mdgrid
    USE tmp_mdhydro
    USE tmp_mdmetric
    USE MD_Boundary

    implicit none

    REAL(KIND=double) :: Mbh,radius,det

    INTEGER :: i,index,m



100 FORMAT (2E22.11)
110 FORMAT (A9,E22.11)



    gamma = 5.d0/3.d0

    DO i=imin,imax+3


        alp(i) = 1.d0

        betaux(i) = 0.d0

        xsi(i) = log(1.d0)

        xsi2(i) = EXP(-2.d0*xsi(i))

        grr(i) = 1.d0

        gtt(i) = grr(i)*xh(i)**2

        gpp(i) = grr(i)*xh(i)**2

        det = grr(i)*gtt(i)*gpp(i)

        sqdetg(i)= sqrt(det)

        urr(i) = 1.d0/grr(i)
        utt(i) = 1.d0/gtt(i)
        upp(i) = 1.d0/gpp(i)


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


        if(i .le. ihalf) then
!            rho(i) = 10.0d0
!            p(i) = 13.3

            rho(i) = 1.0d0
            p(i) = 10.0e+03

        else
!            rho(i) = 1.0d0
!            p(i) = 0.66e-06

            rho(i) = 1.0d0
            p(i) = 10.0e-02

        endif

        velr(i) = 0.0d0

        lorentz(i)= 1.d0/sqrt(1.d0-velr(i)**2)

        eps(i) = p(i)/(gamma-1.d0)*rho(i)

        h(i) = 1 + eps(i) + p(i)/rho(i)


        dens(i) = rho(i)*sqdetg(i)*lorentz(i)

        tau(i) = sqdetg(i)*(lorentz(i)**2*rho(i)*h(i) -p(i))-dens(i)

        srval(i) = sqdetg(i)*lorentz(i)**2*rho(i)*h(i)*velr(i)


    end DO



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

    rho_atm = maxval(rho)*1.d-10

end Subroutine Shock_Tube
