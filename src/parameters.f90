!===============================================================================!
MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE InitializeParameters
  MODULE PROCEDURE InitializeParameters
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: InitializeParameters
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
CONTAINS
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE InitializeParameters()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: TEnd
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: KappaP1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: Reconstruction
USE MOD_FiniteVolume2D_vars,ONLY: ReconstructionFix
USE MOD_FiniteVolume2D_vars,ONLY: WhichOutput
USE MOD_FiniteVolume2D_vars,ONLY: nOutputFiles
USE MOD_FiniteVolume2D_vars,ONLY: InitialCondition
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: VarNameVisu
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

InitialCondition = 211

PrimRefState1 = (/1.4,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0/)
PrimRefState2 = (/1.4,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0/)
PrimRefState3 = (/1.4,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0/)
PrimRefState4 = (/1.4,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0/)

SELECT CASE(InitialCondition)
  CASE(200) ! Constant State
    TEnd    = 1.0
    Kappa   = 1.4
    nElemsX = 100
    nElemsY = 100
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE(211) ! Magnetic Field Loop Advection
    TEnd    = 2.0
    Kappa   = 1.66666666666667
    nElemsX = 800
    nElemsY = 400
    MESH_X0 = (/-1.0,-0.5/)
    MESH_X1 = (/+1.0,+0.5/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE(213) ! Orszag-Tang Vortex
    TEnd    = 1.0
    Kappa   = 1.66666666666667
    nElemsX = 600
    nElemsY = 600
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE(215) ! Rotor Problem
    TEnd    = 0.25
    Kappa   = 1.4
    nElemsX = 600
    nElemsY = 600
    MESH_X0 = (/0.0,0.0/)
    MESH_X1 = (/1.0,1.0/)
    BoundaryConditionsType = (/2,2,2,2/)
  CASE(217) ! Kelvin-Helmholtz Instability
    TEnd    = 10.0
    Kappa   = 1.4
    nElemsX = 600
    nElemsY = 600
    MESH_X0 = (/-0.5,-0.5/)
    MESH_X1 = (/+0.5,+0.5/)
    BoundaryConditionsType = (/1,1,1,1/)
  CASE DEFAULT
    ErrorMessage = "Initial condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

CFL      = 0.95
KappaM1  = Kappa-1.0
KappaP1  = Kappa+1.0
sKappaM1 = 1.0/KappaM1

Reconstruction    = 4
ReconstructionFix = 3

WhichOutput  = 2
nOutputFiles = 100

VarNameVisu(1) = 'Density'
VarNameVisu(2) = 'VelocityX'
VarNameVisu(3) = 'VelocityY'
VarNameVisu(4) = 'VelocityZ'
VarNameVisu(5) = 'Pressure'
VarNameVisu(6) = 'MagneticFieldX'
VarNameVisu(7) = 'MagneticFieldY'
VarNameVisu(8) = 'MagneticFieldZ'
VarNameVisu(9) = 'Psi'

!-------------------------------------------------------------------------------!
END SUBROUTINE InitializeParameters
!===============================================================================!
!
!
!
!===============================================================================!
END MODULE MOD_Parameters
!-------------------------------------------------------------------------------!
