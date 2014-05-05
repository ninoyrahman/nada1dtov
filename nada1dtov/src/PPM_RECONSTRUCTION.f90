! PPM_RECONSTRUCTION.f90
!
! Created on: May 5, 2014
! Author: ninoy

MODULE PPM_RECONSTRUCTION
    IMPLICIT NONE

CONTAINS


    !!!quadratic interpolation for ppm
    SUBROUTINE ppm_quadratic_interpolation(w, wm, wp)

        USE tmp_mdparam
        USE tmp_mdgrid
        implicit none

        integer i
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: w
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(OUT) :: wm, wp
        REAL(KIND=double) :: dmwi, dmwip1, dwi, dwip1, wi12

        wm = 0.d0
        wp = 0.d0

        call OMP_SET_NUM_THREADS(4)

        !$OMP PARALLEL SHARED(w, wm, wp) PRIVATE(i, dwi, dwip1, dmwi, dmwip1, wi12)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-2, imax+1
            dwi = (1.0/2.0)*(w(i+1) - w(i-1))
            dwip1 = (1.0/2.0)*(w(i+2) - w(i))
            if((w(i+1) - w(i))*(w(i) - w(i-1)) .gt. 0.0) then
                dmwi = (dwi/abs(dwi))*min(abs(dwi), 2.0*min(abs(w(i) - w(i-1)), abs(w(i+1) - w(i))))
            else
                dmwi = 0.0
            endif
            if((w(i+2) - w(i+1))*(w(i+1) - w(i)) .gt. 0.0) then
                dmwip1 = (dwip1/abs(dwip1))*min(abs(dwip1), 2.0*min(abs(w(i+1) - w(i)), abs(w(i+2) - w(i+1))))
            else
                dmwip1 = 0.0
            endif
            wi12 = w(i) + 0.5*(w(i+1) - w(i)) + (1.0/(4.0*dx))*(-(2.0/3.0)*dx*dmwip1 + (2.0/3.0)*dx*dmwi)
            wp(i) = wi12
            wm(i) = wi12
        enddo

        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE ppm_quadratic_interpolation

    !!! Monotonicity Preserving Enforcement
    SUBROUTINE ppm_monotonocity_enforcement(w, wtm, wtp, wm, wp)
        USE tmp_mdparam
        USE tmp_mdgrid
        implicit none

        integer i
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: w, wtm, wtp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(OUT) :: wm, wp

        call OMP_SET_NUM_THREADS(4)

        !$OMP PARALLEL SHARED(w, wtm, wtp, wm, wp) PRIVATE(i)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-1, imax
            if((wtm(i) - w(i))*(w(i) - wtp(i-1)) .le. 0.0) then
                wm(i) = w(i)
            else if(-(wtm(i) - wtp(i-1))*(w(i) - 0.5*(wtp(i-1) + wtm(i))) &
                .gt. (wtm(i) - wtp(i-1))*(wtm(i) - wtp(i-1))/6.0) then
                wm(i) = 3.0*w(i) - 2.0*wtp(i-1)
            else
                wm(i) = wtm(i)
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL SHARED(w, wtm, wtp, wm, wp) PRIVATE(i)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-1, imax
            if((wtm(i+1) - w(i+1))*(w(i+1) - wtp(i)) .le. 0.0) then
                wp(i) = w(i+1)
            else if((wtm(i+1) - wtp(i))*(w(i+1) - 0.5*(wtp(i) + wtm(i+1))) &
                .gt. (wtm(i+1) - wtp(i))*(wtm(i+1) - wtp(i))/6.0) then
                wp(i) = 3.0*w(i+1) - 2.0*wtm(i+1)
            else
                wp(i) = wtp(i)
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE ppm_monotonocity_enforcement

    !!! Flattening at Shock
    SUBROUTINE ppm_flattening(w, wm, wp, velr, p)

        USE tmp_mdparam
        USE tmp_mdgrid
        implicit none

        integer i
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: w, velr, p
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(OUT) :: wm, wp
        REAL(KIND=double) :: omega0, omega1, omega2
        REAL(KIND=double) :: fi, fip1, fbari, fbarip1, fbarip2, fbarim1

        omega0 = 0.33
        omega1 = 0.75
        omega2 = 10

        call OMP_SET_NUM_THREADS(4)

        !$OMP PARALLEL SHARED(w, p, velr, wm) PRIVATE(i, fbarim1, fbari, fbarip1, fi)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-1, imax
            if(abs(p(i+1)- p(i-1)) > omega0*min(p(i+1), p(i-1)) .and. velr(i-1) > velr(i+1)) then
                fbari = max(0.0, 1.0 - max(0.0, omega2*((p(i+1) - p(i-1))/(p(i+2) - p(i-2)) - omega1)))
                if((p(i-1) - p(i+1)) > 0) then
                    fbarip1 = max(0.0, 1.0 - max(0.0, omega2*((p(i+2) - p(i))/(p(i+3) - p(i-1)) - omega1)))
                    fi = max(fbari, fbarip1)
                else if((p(i-1) - p(i+1)) < 0) then
                    fbarim1 = max(0.0, 1.0 - max(0.0, omega2*((p(i) - p(i-2))/(p(i+1) - p(i-3)) - omega1)))
                    fi = max(fbari, fbarim1)
                else
                    fi = fbari
                endif
                wm(i) = w(i)*(1.0 - fi) + wm(i)*fi
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL SHARED(w, p, velr, wp) PRIVATE(i, fbari, fbarip1, fbarip2, fip1)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-1, imax
            if(abs(p(i+2)- p(i)) > omega0*min(p(i+2), p(i)) .and. velr(i) > velr(i+2)) then
                fbarip1 = max(0.0, 1.0 - max(0.0, omega2*((p(i+2) - p(i))/(p(i+3) - p(i-1)) - omega1)))
                if((p(i) - p(i+2)) > 0) then
                    fbarip2 = max(0.0, 1.0 - max(0.0, omega2*((p(i+3) - p(i+1))/(p(i+4) - p(i)) - omega1)))
                    fip1 = max(fbarip1, fbarip2)
                else if((p(i) - p(i+2)) < 0) then
                    fbari = max(0.0, 1.0 - max(0.0, omega2*((p(i+1) - p(i-1))/(p(i+2) - p(i-2)) - omega1)))
                    fip1 = max(fbarip1, fbari)
                else
                    fip1 = fbarip1
                endif
                wp(i) = w(i+1)*(1.0 - fip1) + wp(i)*fip1
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE ppm_flattening


    !!! Steepening at Contact Discontinuity
    SUBROUTINE ppm_steeping(rho, rhom, rhop, p)

        USE tmp_mdparam
        USE tmp_mdgrid
        implicit none

        integer i
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: rho, p
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(OUT) :: rhom, rhop
        REAL(KIND=double) :: dpi, dpip1, drhoi, drhoip1
        REAL(KIND=double) :: d2rhoim1, d2rhoip1, d2rhoip2, d2rhoi
        REAL(KIND=double) :: etai, etaip1, etabari, etabarip1
        REAL(KIND=double) :: dmrhoi, dmrhoip1
        REAL(KIND=double) :: K, e, eta1, eta2

        e = 0.01
        eta1 = 20
        eta2 = 0.05
        K = 0.1*gamma

        call OMP_SET_NUM_THREADS(4)

        !$OMP PARALLEL SHARED(rho,p,K,e,eta1,eta2,rhom) PRIVATE(i,dpi,drhoi,drhoip1,d2rhoip1, &
        !$OMP& d2rhoim1,etabari,etai,dmrhoip1)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-1, imax
            drhoi = 0.5*(rho(i+1) - rho(i-1))
            drhoip1 = 0.5*(rho(i+2) - rho(i))
            dpi = 0.5*(p(i+1) - p(i-1))
            if(K*abs(drhoi)/min(rho(i+1), rho(i-1)) >= abs(dpi)/min(p(i+1), p(i-1))) then
                d2rhoip1 = (rho(i+2) -2.0*rho(i+1) + rho(i))/(dx*dx)
                d2rhoim1 = (rho(i) -2.0*rho(i-1) + rho(i-2))/(dx*dx)
                if(d2rhoip1*d2rhoim1 < 0.0 .and. abs(rho(i+1) - rho(i-1)) &
                 .gt. e*min(abs(rho(i+1)), abs(rho(i-1)))) then
                    etabari = (rho(i-2) - rho(i+2) + 4.0*drhoi)/(12.0*drhoi)
                else
                    etabari = 0.0
                endif
                etai = max(0.0, min(eta1*(etabari - eta2), 1.0))
                if((rho(i+2) - rho(i+1))*(rho(i+1) - rho(i)) > 0) then
                    dmrhoip1 = (drhoip1/abs(drhoip1))*min(abs(drhoip1), 2.0*min(abs(rho(i+1) - rho(i)), &
                    abs(rho(i+2) - rho(i+1))))
                else
                    dmrhoip1 = 0.0
                endif
                rhom(i) = rhom(i)*(1.0 - etai) + (rho(i+1) - 0.5*dmrhoip1)
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL SHARED(rho,p,K,e,eta1,eta2,rhop) PRIVATE(i,dpip1,drhoi,drhoip1,d2rhoip2, &
        !$OMP& d2rhoi,etabarip1,etaip1,dmrhoi)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-1, imax
            drhoi = 0.5*(rho(i+1) - rho(i-1))
            drhoip1 = 0.5*(rho(i+2) - rho(i))
            dpip1 = 0.5*(p(i+2) - p(i))
            if(K*abs(drhoip1)/min(rho(i+2), rho(i)) >= abs(dpip1)/min(p(i+2), p(i))) then
                d2rhoip2 = (rho(i+3) -2.0*rho(i+2) + rho(i+1))/(dx*dx)
                d2rhoi = (rho(i+1) -2.0*rho(i) + rho(i-1))/(dx*dx)
                if(d2rhoip2*d2rhoi < 0.0 .and. abs(rho(i+2) - rho(i)) > e*min(abs(rho(i+2)), &
                abs(rho(i)))) then
                    etabarip1 = (rho(i-1) - rho(i+3) + 4.0*drhoip1)/(12.0*drhoip1)
                else
                    etabarip1 = 0.0
                endif
                etaip1 = max(0.0, min(eta1*(etabarip1 - eta2), 1.0))
                if((rho(i+1) - rho(i))*(rho(i) - rho(i-1)) > 0) then
                    dmrhoi = (drhoi/abs(drhoi))*min(abs(drhoi), 2.0*min(abs(rho(i) - rho(i-1)), &
                    abs(rho(i+1) - rho(i))))
                else
                    dmrhoi = 0.0
                endif
                rhop(i) = rhop(i)*(1.0 - etaip1) + (rho(i) + 0.5*dmrhoi)
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE ppm_steeping

    !PPM reconstruction
    SUBROUTINE ppm(rho,rhom,rhop,velr,velrm,velrp,eps,epsm,epsp,p,pressm,pressp)
        USE tmp_mdparam
        USE tmp_mdgrid
        implicit none

        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: rho, velr, eps, p
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(OUT) :: rhom, rhop, velrm, velrp, &
         epsm, epsp, pressm, pressp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx):: rhotm, rhotp, velrtm, velrtp, epstm, &
        epstp, presstm, presstp

        call ppm_quadratic_interpolation(rho, rhotm, rhotp)
        call ppm_quadratic_interpolation(velr, velrtm, velrtp)
        call ppm_quadratic_interpolation(eps, epstm, epstp)
        call ppm_quadratic_interpolation(p, presstm, presstp)

        call ppm_steeping(rho, rhom, rhop, p)

        call ppm_flattening(rho, rhotm, rhotp, velr, p)
        call ppm_flattening(velr, velrtm, velrtp, velr, p)
        call ppm_flattening(eps, epstm, epstp, velr, p)
        call ppm_flattening(p, presstm, presstp, velr, p)

        call ppm_monotonocity_enforcement(rho, rhotm, rhotp, rhom, rhop)
        call ppm_monotonocity_enforcement(velr, velrtm, velrtp, velrm, velrp)
        call ppm_monotonocity_enforcement(eps, epstm, epstp, epsm, epsp)
        call ppm_monotonocity_enforcement(p, presstm, presstp, pressm, pressp)

    END SUBROUTINE ppm

    END MODULE PPM_RECONSTRUCTION

