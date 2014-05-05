  MODULE tmp_mdhydro
  USE tmp_mdparam


  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: rho,vr,eps,lorentz,h
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: velr,p,atmpmask
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: dens,srval,tau

 

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: rhom,vrm,epsm
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: rhop,vrp,epsp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: velrm,pressp,pressm
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: velrp,tm,tp







ENDMODULE tmp_mdhydro

















