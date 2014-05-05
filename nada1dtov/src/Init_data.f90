  SUBROUTINE Init_data
  USE tmp_mdparam 
  USE tmp_mdgrid 
  USE tmp_mdmetric
  USE tmp_mdhydro
  USE MD_Boundary
  USE MD_Subroutines

  IMPLICIT NONE


      WRITE(*,*) ''
      WRITE(*,*) '   Initial data     '
      write(*,*) ' *****************************************'
      WRITE(*,*) ''


      IF (InitSpc.EQ.0) THEN
  
      WRITE(*,*) ''
      WRITE(*,*) 'Setting Minkowski data'
      WRITE(*,*) ''
      call Metric_minkowski

      ELSE IF (InitSpc.EQ.1) THEN


 
      write(*,*)'*************************'
      write(*,*)'    TOV Initial Data     '
      write(*,*)'*************************'
      write(*,*)''

!!$
!!$      call TOV_solver
!!$
!!$      call Detg(gxx,gxy,gxz,gyy,gyz,gzz,det_g)
!!$      sqdetg =sqrt(det_g)
!!$
!!$      psi = det_g**twelfth
!!$      phi = log(psi)
!!$      xsi = 1.d0/psi**4 !EXP(-4*phi)
!!$
!!$      gphixx = EXP(-4*phi)*gxx 
!!$      gphiyy = EXP(-4*phi)*gyy 
!!$      gphizz = EXP(-4*phi)*gzz 
!!$      gphixy = EXP(-4*phi)*gxy 
!!$      gphixz = EXP(-4*phi)*gxz
!!$      gphiyz = EXP(-4*phi)*gyz 
!!$


      write(*,*)'**************************************'
      write(*,*)' Schwarzschild isotropic coordinates  '
      write(*,*)'**************************************'
      write(*,*)''


       call Schwarzschild_isotropic

       

    ELSE IF (InitSpc.EQ.10) THEN


      write(*,*)'**************************************'
      write(*,*)' Schwarzschild isotropic coordinates  '
      write(*,*)' analytic time independent  solution  '
      write(*,*)'**************************************'
      write(*,*)''

!      call Schw_isotrop_solution

!      xsi = EXP(-4.d0*phi)!1.d0/psi**4 !EXP(-4*phi)


     ELSE IF (InitSpc.EQ.11) THEN

      write(*,*)'**************************************'
      write(*,*)' Using checkpoint data                '
      write(*,*)'**************************************'
      write(*,*)''

!      call Read_Checkpoint

     ELSE IF (InitSpc.EQ.12) THEN

      write(*,*)'**************************************'
      write(*,*)'ScalarField                           '
      write(*,*)'**************************************'
      write(*,*)''

      call ScalarField


     ELSE IF (InitSpc.EQ.13) THEN

      write(*,*)'**************************************'
      write(*,*)'Shock                           '
      write(*,*)'**************************************'
      write(*,*)''

      call Shock

   ELSE 

      write(*,*) 'wrong initial data'
      stop

  ENDIF
 
END SUBROUTINE Init_data
