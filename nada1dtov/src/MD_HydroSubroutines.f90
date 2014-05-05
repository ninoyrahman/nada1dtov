MODULE MD_HydroSubroutines
    IMPLICIT NONE

CONTAINS



    !----------------------------------------------------------------------------
    !!     This subroutine computes the source terms for thee hydro equations
    !----------------------------------------------------------------------------


    subroutine source(lorentztmp,rhotmp,velrtmp,epstmp,betauxtmp,alptmp, &
        & krrtmp,ktttmp,kpptmp,urrtmp,utttmp,upptmp, &
        & grrtmp,gtttmp,gpptmp,& 
        & xsitmp,dxdalptmp,dxdbetauxtmp, & 
        & dxdgrrtmp,dxdgtttmp,dxdgpptmp, & 
        & su1,su2,su5)

        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        USE tmp_mdmetric

        implicit none

        integer i
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: rhotmp,epstmp,velrtmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: lorentztmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: betauxtmp,alptmp,xsitmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: krrtmp,ktttmp,kpptmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: urrtmp,utttmp,upptmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: grrtmp,gtttmp,gpptmp

        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: dxdalptmp

        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: dxdbetauxtmp

        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx), INTENT(IN) :: dxdgrrtmp,dxdgtttmp,dxdgpptmp

        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: vrtmp

        REAL(KIND=double), DIMENSION(1:nx), INTENT(OUT) :: su1,su2,su5

 
        REAL(KIND=double):: w,rhoenthalpy,vel1rtmp
        REAL(KIND=double) :: ptmp,t0_r,t00,t0r,t0t,t0p
        REAL(KIND=double) :: trr,ttt,tpp


        su1 = 0.0
        su2 = 0.0
        su5 = 0.0


        !!$OMP PARALLEL
        !!$OMP DO PRIVATE(i,k,w,ptmp,rhoenthalpy,vel1xtmp,vel1ytmp,vel1ztmp,t00,t0x,t0y,t0z, &
        !!$OMP & txx,txy,txz,tyy,tyz,tzz,t0_x,t0_y,t0_z)


        do i = imin, imax
             
            w = lorentztmp(i)
            ptmp = p(i)
            !             enthalpy = (1. + epstmp(i) + ptmp/rhotmp(i))

            rhoenthalpy = rhotmp(i) * (1. + epstmp(i)) + ptmp

            !!$ T^{ij}
            !-----------------------------------------------------------------------------------------
            t00 = (rhoenthalpy*w**2 - ptmp)/(alptmp(i)**2) !OK

            vel1rtmp = alptmp(i)*velrtmp(i)-betauxtmp(i)


            t0r = (rhoenthalpy*w**2*vel1rtmp + ptmp*betauxtmp(i))/ (alptmp(i)**2)!OK
            t0t = 0.d0
            t0p = 0.d0

            trr = rhoenthalpy*w*w*vel1rtmp*vel1rtmp/(alptmp(i)**2) + &
                ptmp * (urrtmp(i) - betauxtmp(i)*betauxtmp(i)/(alptmp(i)**2))

            ttt = ptmp*utttmp(i)

            tpp = ptmp*upptmp(i)


            vrtmp(i) = grrtmp(i)*velrtmp(i)

 
            t0_r = rhoenthalpy*w*w*vrtmp(i)/alptmp(i)

            !-----------------------------------------------------------------------------------------
 

            su1(i) = 0.d0


            su2(i) = t00*(0.5d0*betauxtmp(i)*betauxtmp(i)*dxdgrrtmp(i)-alptmp(i)*dxdalptmp(i))+ &
                & t0r*betauxtmp(i)*dxdgrrtmp(i) + t0_r*dxdbetauxtmp(i) + &
                & 0.5d0*trr*dxdgrrtmp(i)+0.5d0*ttt*dxdgtttmp(i)+0.5d0*tpp*dxdgpptmp(i)

            !!$      if (su2(i).ne.0.d0) then
            !!$      write(*,*) su2(i),i
            !!$      stop
            !!$      end if

            su5(i) = t00*(betauxtmp(i)*betauxtmp(i)*krrtmp(i)- &
                & betauxtmp(i)*dxdalptmp(i)) + &
                & t0r*(-dxdalptmp(i)+2.d0*betauxtmp(i)*krrtmp(i)) + &
                & trr*krrtmp(i)+ttt*ktttmp(i)+tpp*kpptmp(i)


        !!$      if (i.eq.31) then
        !!$
        !!$           write(*,*) ''
        !!$           write(*,*) 'source', su5(i),t00,t0r,trr,betauxtmp(i)
        !!$           write(*,*) ''
        !!$
        !!$           end if


        !      j_r(i) = alptmp(i)*grrtmp(i)*t0r


        end do


    !!$OMP END DO
    !!$OMP END PARALLEL


    end subroutine source

    !reconstrcution:

    SUBROUTINE rlnlx(w, wm, wp)

        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        implicit none

        integer i,j,k
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: w
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(OUT) :: wm,wp
        REAL(KIND=double) :: slope, ss, ssp,c


        ss= 0.d0
        ssp = 0.d0
        wm = 0.d0
        wp = 0.d0
        slope = 0.d0
        c =  0.0d0

        irec = 2

        if (irec .eq. 1) then !doing minmod


            !!$ computing the slope using upwind and downwind method:
            !!$          slope-limiter ---> slope =minmod(ss,ssp) where the minmod function is
            !!$
            !!$                     | ss    if |ss|<|ssp| and ss*ssp>0
            !!$             minmod =| ssp   if |ssp|<|ss| and ss*ssp>0
            !!$                     | 0     if ss*ssp<=0
            !!$

            !$@OMP PARALLEL
            !$@OMP DO PRIVATE(i,k,ss,ssp,slope,wm,wp)
            do i = imin-1, imax

                ss = (w(i)-w(i-1))/dx
                ssp = (w(i+1)-w(i))/dx

                if ((ss*ssp) .gt. 0.d0) then
                    slope=min(abs(ss),abs(ssp))*ss/abs(ss)
                else
                    slope=0.d0
                endif
                !this is the values at the left of the i+1/2 interface
                wm(i)=w(i)+slope*(dx/2.)!(xh_i(i) - xh(i))
            enddo

            !$@OMP END DO
            !$@OMP END PARALLEL

            !    ss= 0.d0
            !    ssp = 0.d0
            !    slope = 0.d0

            !$@OMP PARALLEL
            !$@OMP DO PRIVATE(i,k,ss,ssp,slope,wm,wp)

            do i = imin-1, imax

                ss = (w(i+1)-w(i))/dx
                ssp =(w(i+2)-w(i+1))/dx

                if ((ss*ssp) .gt. 0.d0) then
                    slope=min(abs(ss),abs(ssp))*ss/abs(ss)
                else
                    slope=0.d0
                endif

                !this is the values at the right of the i+1/2 interface

                wp(i)=w(i+1)+slope*(-dx/2.)!(xh_i(i) - xh(i+1))

            enddo


        !$@OMP END DO
        !$@OMP END PARALLEL


        else if (irec.eq.2) then !doing MC slope limiter


            !!$      MC slope-limiter ---> slope =MC(ss,ssp) where the minmod function is
            !!$
            !!$                     | 2ss    if |ss| <|ssp| and 2|ss|<|c|  and ss*ssp>0
            !!$             minmod =| 2ssp   if |ssp|<|ss|  and 2|ssp|<|c| and ss*ssp>0
            !!$                     | c      if |c|  <2|ss| and |c|<2|ssp| and ss*ssp>0
            !!$                     | 0      if ss*ssp<=0
            !!$
            !!$               where c = (ss+ssp)/2


            !$@OMP PARALLEL
            !$@OMP DO PRIVATE(i,k,ss,ssp,slope,c,wm,wp)

            do i = imin-1, imax

                ss  = (w(i)-w(i-1))/dx
                ssp = (w(i+1)-w(i))/dx
                c   = (ss+ssp)*0.5

                if ((ss*ssp) .gt. 0.d0) then

                    slope=min(2.d0*abs(ss),2.d0*abs(ssp),abs(c))*ss/abs(ss)
                !!$             if ((abs(ss).lt.abs(ssp)).and.((2.*abs(ss)).lt.abs(c))) then
                !!$               slope = 2.*ss
                !!$             else if ((abs(ssp).lt.abs(ss)).and.((2.*abs(ssp)).lt.abs(c))) then
                !!$               slope = 2.*ssp
                !!$             else if ((abs(c).lt.(2.*abs(ss))).and.(abs(c).lt.(2.*abs(ssp)))) then
                !!$               slope = c
                !!$              end if
                else
                    slope=0.d0
                endif
                !this is the values at the left of the i+1/2 interface
                wm(i)=w(i)+slope*(dx/2.)!(xh_i(i) - xh(i))

            enddo


            !$@OMP END DO
            !$@OMP END PARALLEL


            ss= 0.d0
            ssp = 0.d0
            slope = 0.d0
            c =  0.0d0


            !$@OMP PARALLEL
            !$@OMP DO PRIVATE(i,k,ss,ssp,slope,c,wm,wp)

            do i = imin-1, imax

                ss = (w(i+1)-w(i))/dx
                ssp =(w(i+2)-w(i+1))/dx
                c   = (ss+ssp)*0.5

                if ((ss*ssp) .gt. 0.d0) then
                    slope=min(2.d0*abs(ss),2.d0*abs(ssp),abs(c))*ss/abs(ss)
                else
                    slope=0.d0
                endif
                !this is the values at the right of the i+1/2 interface
                wp(i)=w(i+1)+slope*(-dx/2.)!(xh_i(i) - xh(i+1))

            enddo


        !$@OMP END DO
        !$@OMP END PARALLEL


        endif



    END SUBROUTINE rlnlx


    SUBROUTINE solverr_hlle(betauxtmp,alptmp,grrtmp,gtttmp,gpptmp, &
        &  urrtmp,utttmp,upptmp,xsitmp,fna,fnr,fne)

        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        USE tmp_mdmetric
  
        implicit none
        integer i,j,k,dir

        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: betauxtmp,alptmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: grrtmp,gtttmp,gpptmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: urrtmp,utttmp,upptmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: xsitmp

        REAL(KIND=double), DIMENSION(0:nx+1),INTENT(OUT) :: fna,fnr,fne

        real(KIND=double):: sqdetgmed,l1m,l2m,l3m,l4m,l5m
        real(KIND=double):: l1p,l2p,l3p,l4p,l5p,Evalmax,Evalmin,Evalmaxmin
        real(KIND=double):: pm, pp, hm, hp, wm, wp
        real(KIND=double):: umedrr,umedtt,umedpp,alpmed, vvm, vvp
        real(KIND=double):: gmedrr,gmedtt,gmedpp,xsimed
        real(KIND=double):: betauxmed
        real(KIND=double):: f1m,f2m,f5m
        real(KIND=double):: f1p,f2p,f5p

        real(KIND=double):: densp,srvalp,taup
        real(KIND=double):: densm,srvalm,taum
        real(KIND=double):: ddens,dsrval,dtau
        real(KIND=double):: velrval
        real(KIND=double):: vrval




        !!$OMP DO PRIVATE(i,k,sqdetgmed,l1m,l2m,l3m,l4m,l5m,l1p,l2p,l3p,l4p,l5p,Evalmax,Evalmin, &
        !!$OMP  Evalmaxmin,pm,pp,hm,hp,wm,wp,umedxx,umedxy,umedxz,umedyy,umedyz,umedzz,alpmed,vvm, &
        !!$OMP  vvp,gmedxx,gmedxy,gmedxz,gmedyy,gmedyz,gmedzz,phimed,betauxmed,betauymed,betauzmed, &
        !!$OMP  f1m,f2m,f3m,f4m,f5m,f1p,f2p,f3p,f4p,f5p,densp,sxvalp,syvalp,szvalp,taup,   &
        !!$OMP  densm,sxvalm,syvalm,szvalm,taum,ddens,dsxval,dsyval,dszval,dtau,velxval, &
        !!$OMP  velyval,velzval,vxval,vyval,vzval,dir)


        DO i=imin-1,imax

            umedrr = 0.5*(urrtmp(i)+urrtmp(i+1))
            umedtt = 0.5*(utttmp(i)+utttmp(i+1))
            umedpp = 0.5*(upptmp(i)+upptmp(i+1))

            gmedrr = 0.5*(grrtmp(i)+grrtmp(i+1))
            gmedtt = 0.5*(gtttmp(i)+gtttmp(i+1))
            gmedpp = 0.5*(gpptmp(i)+gpptmp(i+1))


            sqdetgmed = sqrt(gmedrr*gmedtt*gmedpp)

            alpmed = 0.5*(alptmp(i)+alptmp(i+1))

            betauxmed = 0.5*(betauxtmp(i)+betauxtmp(i+1))

            xsimed = 0.5*(xsitmp(i)+xsitmp(i+1))

            ! u_ij*veli*vj

            velrval = velrm(i)!+betauxmed)/alpmed


            vrval = gmedrr*velrval

            vvm =  vrval*velrval

            wm = 1.d0/sqrt(1.-vvm)
              

            IF(eos_type.eq.0) then !polytropic
                pm = kp*rhom(i)**gamma
            !           pm = pressm(i)
            else if (eos_type.eq.1) then !ideal gas
                pm = (gamma-1.d0)*rhom(i)*epsm(i)
            end IF


            hm = 1 + epsm(i) + pm / rhom(i)

            !Primitive to conservative (Valencia formulation)

            densm  = sqdetgmed*rhom(i)*wm
            srvalm = sqdetgmed*rhom(i)*wm*wm*hm*vrval
            taum   = sqdetgmed*(rhom(i)*wm*wm*hm-pm)-densm

            call Eigenvalues(rhom(i),velrm(i),epsm(i),pm,vvm,&
                & umedrr,alpmed,betauxmed,l1m,l2m,l3m,l4m,l5m)

            velrval = velrp(i)!+betauxmed)/alpmed

            vrval = gmedrr*velrval

            vvp =  vrval*velrval

            wp  = 1.d0/sqrt(1.-vvp)

            IF(eos_type.eq.0) then !polytropic
                pp = kp*rhop(i)**gamma
            !           pp = pressp(i)
            else if (eos_type.eq.1) then !ideal gas
                pp = (gamma-1.d0)*rhop(i)*epsp(i)
            end IF


            hp = 1 + epsp(i) + pp / rhop(i)

            !Primitive to conservative (Valencia formulation)

            densp  = sqdetgmed*rhop(i)*wp
            srvalp = sqdetgmed*rhop(i)*wp*wp*hp*vrval
            taup   = sqdetgmed*(rhop(i)*wp*wp*hp-pp)-densp



            call Eigenvalues(rhop(i),velrp(i),epsp(i),pp,vvp, &
                & umedrr,alpmed,betauxmed,l1p,l2p,l3p,l4p,l5p)

 
            ! Maximum and minimum eigenvalues:
            Evalmin = min(0.d0,l1m,l2m,l3m,l4m,l5m,l1p,l2p,l3p,l4p,l5p)
            Evalmax = max(0.d0,l1m,l2m,l3m,l4m,l5m,l1p,l2p,l3p,l4p,l5p)

            Evalmaxmin = Evalmax - Evalmin
             
            ! friedrichs variables (state_right - state_left)

            ddens  = (densp-densm)
            dsrval = (srvalp-srvalm)
            dtau   = (taup - taum)
 
            ! numerical fluxes at the left of the interface
            call Fluxespm(densm,srvalm,taum,velrm(i),pm,alpmed,betauxmed,sqdetgmed,f1m,f2m,f5m)

            ! numerical fluxes at the right of the interface
            call Fluxespm(densp,srvalp,taup,velrp(i),pp,alpmed,betauxmed,sqdetgmed,f1p,f2p,f5p)

            fna(i) = (Evalmax*f1m - Evalmin*f1p + Evalmax*Evalmin*ddens)/Evalmaxmin
            fnr(i) = (Evalmax*f2m - Evalmin*f2p + Evalmax*Evalmin*dsrval)/Evalmaxmin
            fne(i) = (Evalmax*f5m - Evalmin*f5p + Evalmax*Evalmin*dtau)/Evalmaxmin

        enddo


    !!$OMP END DO


 
    END SUBROUTINE solverr_hlle




    subroutine matter_source(lorentztmp,rhotmp,velrtmp,epstmp,betauxtmp,alptmp, &
        & grrtmp,urrtmp)

        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        USE tmp_mdmetric

        implicit none

        integer i
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: rhotmp,epstmp,velrtmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: lorentztmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: betauxtmp,alptmp
        REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx),INTENT(IN) :: grrtmp,urrtmp


        REAL(KIND=double):: w,rhoenthalpy,vel1rtmp
        REAL(KIND=double) :: ptmp,t0r,ur,utheta



    
        do i = imin, imax

            if (rhotmp(i).gt.1.01*rho_atm) then

                w = lorentztmp(i)
                ptmp = p(i)

                rhoenthalpy = rhotmp(i) * (1. + epstmp(i)) + ptmp


                vel1rtmp = alptmp(i)*velrtmp(i)-betauxtmp(i)



                t0r = (rhoenthalpy*w**2*vel1rtmp + ptmp*betauxtmp(i))/ (alptmp(i)**2)!OK


                !              vel1rtmp = velrtmp(i)-betauxtmp(i)/alptmp(i)

                !             t0r = rhoenthalpy*w**2*vel1rtmp/alptmp(i) + ptmp*(betauxtmp(i)**2/alptmp(i)**2)

                ur = w*(velrtmp(i)-betauxtmp(i)/alptmp(i))

                !             rho_matter(i) = rhoenthalpy*w*w - ptmp !rhotmp(i) *(1. + epstmp(i))
             
                !             S_a(i) = grrtmp(i)*(rhoenthalpy*ur*ur+ptmp*(urrtmp(i)-betauxtmp(i)**2/alptmp(i)**2))!ptmp

                !OK

                rho_matter(i) = rhoenthalpy*w*w - ptmp
                S_a(i) = grrtmp(i)*(rhoenthalpy*velrtmp(i)*velrtmp(i)*w*w +ptmp*urrtmp(i))
                S_b(i) = ptmp

                j_r(i) = rhoenthalpy*w*w*velrtmp(i)*grrtmp(i)

            else


                rho_matter(i) = 0.d0
             
                S_a(i) = 0.d0
                S_b(i) = 0.d0

                j_r(i) = 0.d0


            end if

        end do




    end subroutine matter_source


    Subroutine Constraints
        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        USE tmp_mdmetric

        implicit none

        REAL(KIND=double) :: fac
        REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: R,Hamiltonian

        INTEGER :: i

        do i=imin,imax

            fac = -xsi2(i)**2/(small_a(i))

            R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- &
                (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* &
                (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ &
                4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))) - (1.d0/small_a(i))*8.d0*(dxxdxsi_x2(i)+dxdxsi(i)*dxdxsi_x2(i)) + &
                (1.d0/small_a(i))*8.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i))

       
            Hamiltonian(i)  = R(i) - (A_a(i)**2+0.5*A_a(i)**2)+2.d0*trK(i)**2/3.d0-16.d0*3.14159265358979323d0*rho_matter(i)


        end do

        OPEN(UNIT=21, FILE='data/Hamiltonian.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')

100     FORMAT (2E22.11)
110     FORMAT (A9,E22.11)


        write(21,*)''
        write(21,110)'"time = ',time
        DO i=imin,imax
            write(21,100) xh(i),Hamiltonian(i)
        ENDDO

        close(unit=21)



    end Subroutine Constraints


    Subroutine l2norm_ham
        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        USE tmp_mdmetric

        implicit none

        REAL(KIND=double) :: fac,l2norm,l2norm2
        REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: R,Hamiltonian

        INTEGER :: i

        l2norm = 0.d0
  
        do i=imin,imax

            fac = -xsi2(i)**2/(small_a(i))
            !!$
            !!$       R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- &
            !!$            & (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* &
            !!$            & (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ &
            !!$            & 4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))) - (1.d0/small_a(i))*8.d0*(dxxdxsi_x2(i)+dxdxsi(i)*dxdxsi_x2(i)) + &
            !!$            & (1.d0/small_a(i))*8.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i))



            R(i) = fac*(dxxdsmall_a(i)/(2.d0*small_a(i))+dxxdsmall_b(i)/small_b(i)-small_a(i)*dxddelta_x(i)- &
                (dxdsmall_a(i)/small_a(i))*(dxdsmall_a(i)/small_a(i))+0.5d0*(dxdsmall_b(i)/small_b(i))* &
                (dxdsmall_b(i)/small_b(i))+2.d0/(xh(i)*small_b(i))*(3.d0-small_a(i)/small_b(i))*dxdsmall_b(i)+ &
                4.d0/(xh(i)*xh(i))*(1.d0-small_a(i)/small_b(i))) - (1.d0/small_a(i))*8.d0*(dxxdxsi_x2(i)+dxdxsi(i)*dxdxsi_x2(i)) + &
                (1.d0/small_a(i))*8.d0*dxdxsi_x2(i)*(dxdsmall_a(i)/(2.d0*small_a(i))-dxdsmall_b(i)/(small_b(i))-2.d0/xh(i))
       

            Hamiltonian(i)  = R(i) - (A_a(i)**2+0.5*A_a(i)**2)+2.d0*trK(i)**2/3.d0-16.d0*3.14159265358979323d0*rho_matter(i)

            l2norm = l2norm + Hamiltonian(i)**2

        end do


        OPEN(UNIT=22, FILE='data/l2norm_Hamiltonian.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
    
        write(22,*) time,sqrt(l2norm/imax)
        !    write(22,*) time,sqrt(l2norm/490)
    
        close(unit=22)

        !!$
        !!$    If (I_EX.gt.0) then
        !!$
        !!$       write(*,*) I_EX
        !!$
        !!$       if (I_EX.gt.6) then
        !!$          I_EX = I_EX+10
        !!$          else
        !!$          I_EX = I_EX
        !!$       end if
        !!$
        !!$       l2norm2 = 0.d0
        !!$
        !!$       do i=I_EX,imax-10
        !!$          l2norm2 = l2norm2 + abs(Hamiltonian(i))**2
        !!$       end do
        !!$


        l2norm2 = 0.d0
        I_EX = 0
        do i=1,imax-50

            If (xh(i).gt.0.7) then
                l2norm2 = l2norm2 + abs(Hamiltonian(i))**2
                I_EX = I_EX+1
            end If

        end do

        OPEN(UNIT=23, FILE='data/l2norm_Hamiltonian2.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')
        write(23,*) time,sqrt(l2norm2/(I_EX))
        close(unit=23)


    !!$    end if


    end Subroutine l2norm_ham







    Subroutine AH_Finder
        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdmetric

        implicit none

        REAL(KIND=double) :: AHrad1d,AHarea1d,Mirr1d,PI_D,dgttdr
        REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: R,AHrad,psitmp,dpsidr,CapitalB
        integer:: i

        i_EX = 0.d0
        !-------------------------------------------------------------
        ! 1-D AH finder:
        !-------------------------------------------------------------
        PI_D = 3.141592653589793238462d0

        !         psitmp = EXP(xsi)

        psitmp = exp(4.d0*xsi)*small_b

        do i=imin+2,imax-10!400

            !        dgttdr= (gtt(i+1)-gtt(i-1))/(2.*dx)

            !         dgttdr= (gtt(i-2)-8.d0*gtt(i-1)+ 8.d0*gtt(i+1)-gtt(i+2))*inv12dx

            !         AHrad(i)=dgttdr/(gtt(i)*sqrt(grr(i)))-2.d0*ktt(i)/gtt(i)


            dpsidr(i) = (psitmp(i-2)-8.d0*psitmp(i-1)+ 8.d0*psitmp(i+1)-psitmp(i+2))*inv12dx

            AHrad(i)= (2.d0/xh(i) + dpsidr(i)/psitmp(i))/sqrt(small_a(i)*exp(4.d0*xsi(i))) - 2.d0*utt(i)*ktt(i)


        end do

        do i=imin+2,imax-12!398

            if ((AHrad(i).lt.0.d0).and.(AHrad(i+1).gt.0.d0)) then

                AHrad1d=xh(i)

                AHarea1d=4.d0*PI_D*AHrad1d**2*EXP(xsi(i))**4

                Mirr1d=sqrt((4.d0*PI_D*AHrad1d**2*EXP(xsi(i))**4)/(16.d0*PI_D))

 
                if (i.gt.6) then

                    I_EX = i-10

                else
             
                    I_EX = i

                end if

            !          goto 20

            !!$          write(*,*) ''
            !!$          write(*,*) '1D AH: '
            !!$          write(*,*) 'radius AH ', AHrad1d
            !!$          write(*,*) 'Area AH ', AHarea1d
            !!$          write(*,*) 'Irreducible Mass AH ',Mirr1d
            !!$          write(*,*) ''

            end if

        end do

        !20   continue


        OPEN(UNIT=12, FILE='data/AH.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')


        write(12,*) time,I_EX,AHrad1d,AHarea1d,Mirr1d


        close(unit=12)


    !    stop

    end Subroutine AH_Finder










    Subroutine BH
        USE tmp_mdparam
        USE tmp_mdgrid
        USE tmp_mdhydro
        USE tmp_mdmetric

        implicit none

        REAL(KIND=double) :: Mbh,radius,dr
        REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: phi_iso,xsi2_iso,error_xsi2
        REAL(KIND=double), DIMENSION(1:100000) :: psi_iso,r_iso,Rtmp
        REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: psi_iso2

        INTEGER :: i,index,m

        Mbh = 1.d0

100     FORMAT (2E22.11)
110     FORMAT (A9,E22.11)

        dr = 0.001

        do i=1,100000
 
            Rtmp(i) = Rtmp(i-1)+dr

   
            r_iso(i) = ((2.d0*Rtmp(i)+Mbh+(4.d0*Rtmp(i)**2+4.d0*Mbh*Rtmp(i)+3.d0*Mbh**2)**0.5)*0.25d0)* &
                & ( (4.d0+3.d0*sqrt(2.d0))*(2.d0*Rtmp(i)-3.d0*Mbh)/(8.d0*Rtmp(i)+6.d0*Mbh+3.d0* &
                & (8.d0*Rtmp(i)**2+8.d0*Mbh*Rtmp(i)+6.d0*Mbh**2)**0.5))**(1.d0/sqrt(2.d0))
   
            psi_iso(i) = ((4.d0*Rtmp(i)/(2.d0*Rtmp(i)+Mbh+(4.d0*Rtmp(i)**2+4.d0*Mbh*Rtmp(i)+3.d0*Mbh**2)**0.5))**0.5)*&
                & ( (8.d0*Rtmp(i)+6.d0*Mbh+3.d0*(8.d0*Rtmp(i)**2.d0+8.d0*Mbh*Rtmp(i)+6.d0*Mbh**2)**0.5)/ &
                & ((4.d0+3.d0*sqrt(2.d0))*(2.d0*Rtmp(i)-3.d0*Mbh)))**(1.d0/(2.d0*sqrt(2.d0)))




        !!$   phi_iso(i) = log(psi_iso(i))
        !!$
        !!$
        !!$   xsi2_iso(i) = EXP(-2.d0*phi_iso(i))
        !!$
        !!$   error_xsi2(i) = abs(xsi2(i)-xsi2_iso(i))

        end do


        DO i=imin,imax
    
            radius = xh(i)

            do m = 1,100000
           
                if (r_iso(m).gt.radius) then
                    index = m
                    exit
                end if
             
            end do


            psi_iso2(i) = psi_iso(index-1)*(radius-r_iso(index))/(r_iso(index-1)-r_iso(index)) + &
                & psi_iso(index)*(radius-r_iso(index-1))/(r_iso(index)-r_iso(index-1))


            phi_iso(i) = log(psi_iso2(i))


            xsi2_iso(i) = EXP(-2.d0*phi_iso(i))

            error_xsi2(i) = abs(xsi2(i)-xsi2_iso(i))


        end DO



        OPEN(UNIT=15, FILE='data/Error_xsi2.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')


        write(15,*)''
        write(15,110)'"time = ',time
        DO i=imin-ghzx,imax+ghzx
            write(15,*) xh(i),error_xsi2(i)
        ENDDO

        OPEN(UNIT=16, FILE='data/xsi2_sol.dat',STATUS='UNKNOWN',ACTION= 'WRITE',POSITION='append')


        write(16,*)''
        write(16,110)'"time = ',time
        DO i=imin-ghzx,imax+ghzx
            write(16,*) xh(i),xsi2_iso(i)
        ENDDO



        close(unit=15)
        close(unit=16)


    end Subroutine BH



END MODULE MD_HydroSubroutines


