MODULE tmp_mdparam
    IMPLICIT NONE
    SAVE

    !******************************************************
    ! Pedro Montero - 27/11/11
    !------------------------------------------------------
    ! This module defines different parameters that some
    ! of them will go into a paramater file
    !******************************************************
  
    INTEGER, PARAMETER :: single = 4
    !  INTEGER, PARAMETER :: double = 8
    INTEGER, PARAMETER :: double = KIND(1.0D0)
 
    REAL(KIND=double), PARAMETER :: half = 0.5d0
    REAL(KIND=double), PARAMETER :: forth = 0.25d0
    REAL(KIND=double), PARAMETER :: twelfth = 1.d0/12.d0
    REAL(KIND=double), PARAMETER :: twothird = 2.d0/3.d0
    REAL(KIND=double), PARAMETER :: untercio = 1.d0/3.d0
    REAL(KIND=double), PARAMETER :: threehalf = 3.d0/2.d0
    REAL(KIND=double), PARAMETER :: fourthird = 4.d0/3.d0
  
    INTEGER, PARAMETER :: zero = 0
    INTEGER, PARAMETER :: one = 1

!    INTEGER, PARAMETER :: eos_type =0
    INTEGER, PARAMETER :: eos_type = 1

    INTEGER, PARAMETER ::  eps1= -0.001d0

    INTEGER, PARAMETER :: nx = 400

    INTEGER, PARAMETER :: ghzx = 3

    INTEGER, PARAMETER :: imin = 1
    INTEGER, PARAMETER :: imax = nx
    INTEGER, PARAMETER :: ihalf = nx/2

    CHARACTER (len=10),PARAMETER :: phi_bd= 'phi_bd'
    CHARACTER (len=10),PARAMETER :: sphi_bd= 'sphi_bd'
    CHARACTER (len=10),PARAMETER :: xsi_bd= 'xsi_bd'
    CHARACTER (len=10),PARAMETER :: trK_bd= 'trK_bd'
    CHARACTER (len=10),PARAMETER :: gphi_bd= 'gphi_bd'
    CHARACTER (len=10),PARAMETER :: Aphi_bd= 'Aphi_bd'
    CHARACTER (len=10),PARAMETER :: CCF_bd= 'CCF_bd'
    CHARACTER (len=10),PARAMETER :: alp_bd= 'alp_bd'
    CHARACTER (len=10),PARAMETER :: beta_bd= 'beta_bd'


    REAL(KIND=double), PARAMETER :: g_grav=6.673d-8 !# gravitational constant in cm^3 g^-1 s^-2
    REAL(KIND=double), PARAMETER :: c_light=2.99792458d10 !# speed of light in cm s^-1



!d_light_year=9.46053d17           # light year in cm
!d_parsec=3.2616                   # parsec in light years
!d_detect=1.0d7                    # gravitational wave detection distance in parsec
!m_solar=1.989d33                  # solar mass in gm
!unit_charge=1.60217733d-19        # unit charge in C
!k_boltzmann=1.380658d-16          # Boltzmann's constant in erg K^-1
!h_planck=6.626176d-27             # Planck's constant in erg s
!n_avogadro=6.0221415d23           # Avogadro's number
!m_u=1.6605402d-24                 # atomic mass unit in g
!rho_nuc_lorene=1.66d14            # nuclear matter density of the Lorene code in cgs units
!rho_trapping=2.0d12               # neutrino trapping density threshold



ENDMODULE tmp_mdparam
