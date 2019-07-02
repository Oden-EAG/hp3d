! Optimally blended quadrature rules for wave propagation (see Ainsworth, Wajid, 2009)

MODULE ainsworth_quadrature_module
  IMPLICIT NONE
  LOGICAL :: INITIALIZED = .FALSE.
  DOUBLE PRECISION :: XI_AINSWORTH_SCALED(10,10)
  DOUBLE PRECISION :: W_AINSWORTH_SCALED(10,10)
  DOUBLE PRECISION, PARAMETER :: XI_AINSWORTH(10,10) = reshape( &
       (/ &
       0.d0, &
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       0.816496580927726d0, -0.816496580927726d0, & ! p = 1
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       0.0d0, 0.930949336251263d0, -0.930949336251263d0, & ! p = 2
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       -0.964335275879562d0, 0.964335275879562d0, -0.429352058315787d0, &
       0.429352058315787d0, & ! p = 3
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       0.000000000000000d0, -0.978315678013417d0, 0.978315678013417d0, &
       -0.638731398345590d0, 0.638731398345590d0, & ! p = 4
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       -0.985446820998315d0, -0.752558388054789d0, -0.280556681820821d0, &
       0.280556681820821d0, 0.752558388054789d0, 0.985446820998315d0, & ! p = 5
       0.d0, 0.d0, 0.d0, 0.d0, &
       -0.989564901331163d0, -0.820496210793208d0, -0.463335883847022d0, &
       0.000000000000000d0, 0.463335883847022d0, 0.820496210793208d0, &
       0.989564901331163d0, & ! p = 6
       0.d0, 0.d0, 0.d0, &
       -0.992154829409481d0, -0.864059339845500d0, -0.586467949432683d0, &
       -0.207447135295099d0, 0.207447135295099d0, 0.586467949432683d0, &
       0.864059339845500d0, 0.992154829409481d0, & ! p = 7
       0.d0, 0.d0, &
       -0.993888391435683d0, -0.893581190017803d0, -0.672520467240063d0, &
       -0.360613635820052d0, 0.000000000000000d0, 0.360613635820052d0, &
       0.672520467240063d0, 0.893581190017803d0, 0.993888391435683d0, & ! p = 8
       0.d0, &
       -0.995105205867138d0, -0.914477642987090d0, -0.734696655470703d0, &
       -0.475285190219825d0, -0.164365837601352d0, 0.164365837601352d0, &
       0.475285190219825d0, 0.734696655470703d0, 0.914477642987090d0, &
       0.995105205867138d0 &
       /), (/10,10/) )

  DOUBLE PRECISION, PARAMETER :: W_AINSWORTH(10,10) = reshape( &
       (/ &
       0.d0, &
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       1.d0, 1.d0, & ! p = 1
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       1.230769230769231d0, 0.384615384615385d0, 0.384615384615385d0, & ! p = 2
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       0.199826014447922d0, 0.199826014447922d0, 0.800173985552078d0, &
       0.800173985552078d0, & ! p = 3
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       0.693766937669377d0, 0.121787277062268d0, 0.121787277062268d0, &
       0.531329254103044d0, 0.531329254103044d0, & ! p = 4
       0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
       0.081835174559638d0, 0.372395751222672d0, 0.545769074217690d0, &
       0.545769074217690d0, 0.372395751222672d0, 0.081835174559638d0, & ! p = 5
       0.d0, 0.d0, 0.d0, 0.d0, &
       0.058719277436163d0, 0.273663854856191d0, 0.426675691237058d0, &
       0.481882352941176d0, 0.426675691237058d0, 0.273663854856191d0, &
       0.058719277436163d0, & ! p = 6
       0.d0, 0.d0, 0.d0, &
       0.044164752444346d0, 0.208912408709088d0, 0.338113373846498d0, &
       0.408809465000068d0, 0.408809465000068d0, 0.338113373846498d0, &
       0.208912408709088d0, 0.044164752444346d0, & ! p = 7
       0.d0, 0.d0, &
       0.034415599386995d0, 0.164411441988856d0, 0.272653718847312d0, &
       0.344040703534754d0, 0.368957072484166d0, 0.344040703534754d0, &
       0.272653718847312d0, 0.164411441988856d0, 0.034415599386995d0, & ! p = 8
       0.d0, &
       0.027569114782485d0, 0.132615801678385d0, 0.223654050135981d0, &
       0.290430742781218d0, 0.325730290621930d0, 0.325730290621930d0, &
       0.290430742781218d0, 0.223654050135981d0, 0.132615801678385d0, &
       0.027569114782485d0 &
       /), (/10,10/) )
END MODULE ainsworth_quadrature_module

SUBROUTINE init_ainsworth_quadrature
  USE ainsworth_quadrature_module
  IMPLICIT NONE
  INTEGER :: i, j

  DO i = 2, 10
     DO j = 1, i
        XI_AINSWORTH_SCALED(j,i) = (1.0 + XI_AINSWORTH(j,i)) / 2.0
        W_AINSWORTH_SCALED(j,i) = W_AINSWORTH(j,i) / 2.0
     ENDDO
  ENDDO
  INITIALIZED = .TRUE.
END SUBROUTINE init_ainsworth_quadrature

SUBROUTINE ainsworth_quadrature_3D(norder, n_xi, xi, w)
  ! For hexahedral elements only.
  USE ainsworth_quadrature_module
  USE parameters, ONLY: MAXbrickH
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: norder(19) ! order of hexa nodes
  INTEGER, INTENT(OUT) :: n_xi ! total number of integration points
  DOUBLE PRECISION, INTENT(OUT) :: xi(3,MAXbrickH) ! integration points
  DOUBLE PRECISION, INTENT(OUT) :: w(MAXbrickH) ! integration weights
  INTEGER :: nordx, nordy, nordz, nordxy, n, i, j, k

  IF (.NOT. INITIALIZED) THEN
     CALL init_ainsworth_quadrature
  ENDIF
  CALL decode(norder(19), nordxy, nordz)
  CALL decode(nordxy, nordx, nordy)
  n_xi = (nordx+1)*(nordy+1)*(nordz+1)
  n = 0
  DO i = 1, nordx+1
     DO j = 1, nordy+1
        DO k = 1, nordz+1
           n = n + 1
           xi(1, n) = XI_AINSWORTH_SCALED(i, nordx+1)
           xi(2, n) = XI_AINSWORTH_SCALED(j, nordy+1)
           xi(3, n) = XI_AINSWORTH_SCALED(k, nordz+1)
           w(n) = W_AINSWORTH_SCALED(i,nordx+1) * W_AINSWORTH_SCALED(j,nordy+1) * W_AINSWORTH_SCALED(k,nordz+1)
        ENDDO
     ENDDO
  ENDDO
END SUBROUTINE ainsworth_quadrature_3D
