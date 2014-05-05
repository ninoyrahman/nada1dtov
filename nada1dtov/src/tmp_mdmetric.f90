  MODULE tmp_mdmetric
  USE tmp_mdparam


   REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: psi

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: xsi,dxdxsi,dxxdxsi,capB, RHSdeltar

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: dxdxsi_x2,dxxdxsi_x2

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: xsi2,dxdxsi2,dxxdxsi2

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: alp,dxdalp,dxxdalp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: betaux,dxdbetaux,dxxdbetaux,dxdbetaux_x

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: small_a,dxdsmall_a,dxxdsmall_a
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: small_b,dxdsmall_b,dxxdsmall_b

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: delta_x,dxddelta_x
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: trK,dxdtrK,dxdab2,A_a

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: laplacian_alp,Confdiv_beta,dxdConfdiv_beta

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: Cov_alp,dxdA_a
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: rho_matter, S_a, S_b,SS,j_r

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: grr,gtt,gpp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: urr,utt,upp
  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: krr,ktt,kpp

  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: dxdgrr,dxdgtt,dxdgpp,sqdetg

!  REAL(KIND=double), DIMENSION(1-ghzx:nx+ghzx) :: 

!  REAL(KIND=double), DIMENSION(-ghzx+1:nx+ghzx) :: 

  integer :: I_EX


ENDMODULE tmp_mdmetric

















