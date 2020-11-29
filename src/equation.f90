!===============================================================================!
MODULE MOD_Equation
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
PRIVATE
!-------------------------------------------------------------------------------!
INTERFACE ExactFunction
  MODULE PROCEDURE ExactFunction
END INTERFACE

INTERFACE SourceTerms
  MODULE PROCEDURE SourceTerms
END INTERFACE

INTERFACE BoundaryConditions
  MODULE PROCEDURE BoundaryConditions
END INTERFACE

INTERFACE TimeStep
  MODULE PROCEDURE TimeStep
END INTERFACE

INTERFACE RiemannSolver
  MODULE PROCEDURE RiemannSolver
END INTERFACE

INTERFACE EvaluateFlux1D
  MODULE PROCEDURE EvaluateFlux1D
END INTERFACE

INTERFACE ConsToPrim
  MODULE PROCEDURE ConsToPrim
END INTERFACE

INTERFACE PrimToCons
  MODULE PROCEDURE PrimToCons
END INTERFACE

INTERFACE MixedGLM
  MODULE PROCEDURE MixedGLM
END INTERFACE
!-------------------------------------------------------------------------------!
PUBLIC :: ExactFunction
PUBLIC :: SourceTerms
PUBLIC :: BoundaryConditions
PUBLIC :: TimeStep
PUBLIC :: RiemannSolver
PUBLIC :: EvaluateFlux1D
PUBLIC :: ConsToPrim
PUBLIC :: PrimToCons
PUBLIC :: MixedGLM
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
SUBROUTINE ExactFunction(WhichInitialCondition,t,x,Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nDims
USE MOD_FiniteVolume2D_vars,ONLY: PI
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MESH_SX
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: KappaP1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
INTEGER,INTENT(IN) :: WhichInitialCondition
REAL,INTENT(IN)    :: t
REAL,INTENT(IN)    :: x(1:nDims)
REAL,INTENT(OUT)   :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                :: Prim(1:nVar)
REAL                :: Cons_In(1:nVar),Cons_Out(1:nVar)
REAL                :: r, r0, rho0, p0, B0, vs, rho1, rho2, am, om 
REAL                :: xc(2), xm, ym
REAL                :: rho_in,rho_out,v_in,v_out,p_in,p_out,B_in,B_out
REAL                :: delta_rho, delta_vx, delta_vy, delta_p, delta_bx, delta_by
CHARACTER(LEN=255)  :: ErrorMessage
!-------------------------------------------------------------------------------!

Cons = 0.0
Prim = 0.0
SELECT CASE(WhichInitialCondition)
  !----------------------------------------------------------------------!
  ! [200] Constant State                                                 !
  !----------------------------------------------------------------------!
  CASE(200)
    Prim(1:nVar) = PrimRefState1(1:nVar)

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [211] Magnetic Field Loop Advection                                  !
  !----------------------------------------------------------------------!
  CASE(211)
    xm    = MESH_X0(1)+0.5*MESH_SX(1)
    ym    = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1) = x(1)-xm
    xc(2) = x(2)-ym
    r     = SQRT(xc(1)**2 + xc(2)**2)

    rho0  = 1.0
    B0    = 1.0E-03
    p0    = 1.0

    Prim(1) = rho0
    Prim(2) = 2.00
    Prim(3) = 1.00
    Prim(4) = 1.00
    Prim(5) = p0
    Prim(6) = 0.0
    Prim(7) = 0.0
    Prim(8) = 0.0

    IF (r .LT. 0.3) THEN
      Prim(6) = -B0*x(2)/r
      Prim(7) =  B0*x(1)/r
      Prim(8) = 0.0
    END IF

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [213] Orszag-Tang Vortex                                             !
  !----------------------------------------------------------------------!
  CASE(213)
    Prim(1)= Kappa**2
    Prim(2)= -SIN(2.0*PI*x(2))
    Prim(3)= +SIN(2.0*PI*x(1))
    Prim(4)= 0.0
    Prim(5)= Kappa
    Prim(6)= -SIN(2.0*PI*x(2))
    Prim(7)= +SIN(4.0*PI*x(1))
    Prim(8)= 0.0

    CALL PrimToCons(Prim,Cons)
  !----------------------------------------------------------------------!
  ! [215] Rotor Problem                                                  !
  !----------------------------------------------------------------------!
  CASE(215)
    Cons_In  = 0.0
    Cons_Out = 0.0

    xm      = MESH_X0(1)+0.5*MESH_SX(1)
    ym      = MESH_X0(2)+0.5*MESH_SX(2)
    xc(1)   = x(1)-xm
    xc(2)   = x(2)-ym
    r       = SQRT(xc(1)**2 + xc(2)**2)

    r0      = 1.0E-1
    om      = 1.0E+0

    rho_in  = 1.0E+1
    rho_out = 1.0E+0
    v_in    = 0.0E+0
    v_out   = 0.0E+0
    p_in    = 1.0E+0
    p_out   = 1.0E+0
    B_in    = 0.5*SQRT(2.0)
    B_out   = 0.5*SQRT(2.0)

    Prim(1) = rho_in
    Prim(2) = -om*xc(2)/r0
    Prim(3) = +om*xc(1)/r0
    Prim(4) = v_in
    Prim(5) = p_in
    Prim(6) = B_in
    Prim(7) = 0.0
    Prim(8) = 0.0

    CALL PrimToCons(Prim,Cons_In)

    Prim(1) = rho_out
    Prim(2) = v_out
    Prim(3) = v_out
    Prim(4) = v_out
    Prim(5) = p_out
    Prim(6) = B_out
    Prim(7) = 0.0
    Prim(8) = 0.0
    CALL PrimToCons(Prim,Cons_Out)
    Cons    = -0.5*(Cons_In-Cons_Out)*TANH(80.0*(r-r0)) + Cons_Out + 0.5*(Cons_In-Cons_Out)
  !----------------------------------------------------------------------!
  ! [217] Kelvin-Helmholtz Instability                                   !
  !----------------------------------------------------------------------!
  CASE(217)
    vs = 0.50
    am = 0.01

    rho1 = 1.0
    rho2 = 2.0
    p0   = 2.5
    B0   = 0.2

    IF (ABS(x(2)) .GE. 0.25) THEN
      delta_rho = rho1
      delta_vx  = +(vs+am*SIN(2.0*PI*x(1)))
    ELSE
      delta_rho = rho2
      delta_vx  = -(vs+am*SIN(2.0*PI*x(1)))
    END IF
    IF (x(2) .GE. 0.00) THEN
      delta_vy  = +am*SIN(2.0*PI*x(1))
    ELSE
      delta_vy  = +am*SIN(2.0*PI*x(1))
    END IF

    Prim(1) = delta_rho
    Prim(2) = delta_vx
    Prim(3) = delta_vy
    Prim(4) = 0.0
    Prim(5) = p0
    Prim(6) = B0
    Prim(7) = 0.0
    Prim(8) = 0.0
    CALL PrimToCons(Prim,Cons)
  CASE DEFAULT
    ErrorMessage = "Exact function not specified"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!-------------------------------------------------------------------------------!
END SUBROUTINE ExactFunction
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE SourceTerms(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: S
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER          :: ii, jj
!-------------------------------------------------------------------------------!

S = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE SourceTerms
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE BoundaryConditions(t)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: nGhosts
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X0
USE MOD_FiniteVolume2D_vars,ONLY: MESH_X1
USE MOD_FiniteVolume2D_vars,ONLY: MeshBary
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState1
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState2
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState3
USE MOD_FiniteVolume2D_vars,ONLY: PrimRefState4
USE MOD_FiniteVolume2D_vars,ONLY: BoundaryConditionsType
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: V
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)    :: t
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
INTEGER            :: ii, jj
INTEGER            :: idx_vx, idx_vy
REAL               :: x0, xc, xt
REAL               :: Prim_in(1:nVar), Prim_out(1:nVar)
REAL               :: Cons_in(1:nVar), Cons_out(1:nVar)
REAL               :: ConsRefState1(1:nVar), ConsRefState2(1:nVar)
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

idx_vx = 2
idx_vy = 3

!------------------------------!
! Left Boundary Conditions     !
!------------------------------!
SELECT CASE(BoundaryConditionsType(4))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nElemsX-nGhosts+ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState4(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,-nGhosts+ii,jj) = U(1:nVar,nGhosts-ii+1,jj)
        U(idx_vx,-nGhosts+ii,jj) =-U(idx_vx,nGhosts-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Right Boundary Conditions    !
!------------------------------!
SELECT CASE(BoundaryConditionsType(2))
  CASE(1) ! Periodic
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO jj=1,nElemsY    
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState2(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO jj=1,nElemsY
      DO ii=1,nGhosts
        U(1:nVar,nElemsX+ii,jj) = U(1:nVar,nElemsX-ii+1,jj)
        U(idx_vx,nElemsX+ii,jj) =-U(idx_vx,nElemsX-ii+1,jj)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Top Boundary Conditions      !
!------------------------------!
SELECT CASE(BoundaryConditionsType(3))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState3(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,nElemsY+jj) = U(1:nVar,ii,nElemsY-jj+1)
        U(idx_vy,ii,nElemsY+jj) =-U(idx_vy,ii,nElemsY-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

!------------------------------!
! Bottom Boundary Conditions   !
!------------------------------!
SELECT CASE(BoundaryConditionsType(1))
  CASE(1) ! Periodic
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nElemsY-nGhosts+jj)
      END DO
    END DO
  CASE(2) ! Transmissive
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE(3) ! Inflow
    Prim_in(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_in(1:nVar),Cons_in(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_in(1:nVar)
      END DO
    END DO
  CASE(4) ! Outflow
    Prim_out(1:nVar) = PrimRefState1(1:nVar)
    CALL PrimToCons(Prim_out(1:nVar),Cons_out(1:nVar))
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = Cons_out(1:nVar)
      END DO
    END DO
  CASE(5) ! Reflecting
    DO ii=1,nElemsX
      DO jj=1,nGhosts
        U(1:nVar,ii,-nGhosts+jj) = U(1:nVar,ii,nGhosts-jj+1)
        U(idx_vy,ii,-nGhosts+jj) =-U(idx_vy,ii,nGhosts-jj+1)
      END DO
    END DO
  CASE DEFAULT
    ErrorMessage = "Boundary condition not implemented"
    WRITE(*,*) ErrorMessage
    STOP
END SELECT

DO jj=-nGhosts,nElemsY+nGhosts+1
  DO ii=-nGhosts,nElemsX+nGhosts+1
    CALL ConsToPrim(U(1:nVar,ii,jj),V(1:nVar,ii,jj))
  END DO
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE BoundaryConditions
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION TimeStep()
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: MESH_DX
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nElemsX
USE MOD_FiniteVolume2D_vars,ONLY: nElemsY
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxX
USE MOD_FiniteVolume2D_vars,ONLY: LambdaMaxY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_TIMESTEP
USE MOD_FiniteVolume2D_vars,ONLY: GLM_ch
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL    :: TimeStep
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL    :: FastestWaveX, FastestWaveY
REAL    :: Prim(1:nVar)
INTEGER :: ii, jj
!-------------------------------------------------------------------------------!

GLM_ch = 0.0
LambdaMaxX = 0.0
LambdaMaxY = 0.0
TimeStep = HUGE(1.0)

DO jj=1,nElemsY
  DO ii=1,nElemsX
    CALL ConsToPrim(U(1:nVar,ii,jj),Prim(1:nVar))
    CALL WaveSpeeds2D(Prim(1:nVar),FastestWaveX,FastestWaveY)
    LambdaMaxX = MAX(LambdaMaxX,ABS(FastestWaveX))
    LambdaMaxY = MAX(LambdaMaxY,ABS(FastestWaveY))
    TimeStep  = MIN(TimeStep,MESH_DX(1)/LambdaMaxX,MESH_DX(2)/LambdaMaxY)
    GLM_ch = MAX(GLM_ch,LambdaMaxX,LambdaMaxY)
  END DO
END DO

TimeStep = CFL*TimeStep

IF (TimeStep .LT. MIN_TIMESTEP) THEN
  TimeStep = MIN_TIMESTEP
END IF

!-------------------------------------------------------------------------------!
END FUNCTION TimeStep
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds1D(Prim,slowest,fastest)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)           :: Prim(1:nVar)
REAL,INTENT(OUT),OPTIONAL :: slowest
REAL,INTENT(OUT),OPTIONAL :: fastest
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL                      :: rho, vx, vy, vz, p, Bx, By, Bz
REAL                      :: c, v2, B2, c2
REAL                      :: cf, cs, ca2, cf2
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)

c  = EOS_SoundSpeed(rho,p)
c2 = c**2
v2 = vx*vx + vy*vy + vz*vz
B2 = Bx*Bx + By*By + Bz*Bz

ca2 = Bx**2/rho
cf2 = 0.5*((c2+B2/rho) + SQRT((c2+B2/rho)**2.0 - 4.0*c2*ca2))
cf  = SQRT(cf2)

!--------------------!
! Slowest wave       !
!--------------------!
IF (PRESENT(slowest)) THEN
  slowest = vx - cf
END IF

!--------------------!
! Fastest wave       !
!--------------------!
IF (PRESENT(fastest)) THEN
  fastest = vx + cf
END IF

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE WaveSpeeds2D(Prim,fastestx,fastesty)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: fastestx
REAL,INTENT(OUT) :: fastesty
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: c, v2, B2, c2
REAL             :: cfx, cfy, cax2, cay2, cfx2, cfy2
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)

c  = EOS_SoundSpeed(rho,p)
c2 = c**2
v2 = vx*vx + vy*vy + vz*vz
B2 = Bx*Bx + By*By + Bz*Bz

cax2 = Bx**2/rho
cfx2 = 0.5*((c2+B2/rho) + SQRT((c2+B2/rho)**2.0 - 4.0*c2*cax2))
cfx  = SQRT(cfx2)

cay2 = By**2/rho
cfy2 = 0.5*((c2+B2/rho) + SQRT((c2+B2/rho)**2.0 - 4.0*c2*cay2))
cfy  = SQRT(cfy2)

fastestx = vx + cfx
fastesty = vy + cfy

!-------------------------------------------------------------------------------!
END SUBROUTINE WaveSpeeds2D
!===============================================================================!
!
!
!
!===============================================================================!
FUNCTION EOS_SoundSpeed(rho,p)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN) :: rho, p
REAL            :: EOS_SoundSpeed
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
CHARACTER(LEN=255) :: ErrorMessage
!-------------------------------------------------------------------------------!

EOS_SoundSpeed = ABS(Kappa*p/rho)
IF (EOS_SoundSpeed .LT. 0.0) THEN
  ErrorMessage = "Negative speed of sound"
  WRITE(*,*) ErrorMessage
  STOP
END IF
EOS_SoundSpeed = SQRT(EOS_SoundSpeed)

!-------------------------------------------------------------------------------!
END FUNCTION EOS_SoundSpeed
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE ConsToPrim(Cons, Prim)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_ENERGY
USE MOD_FiniteVolume2D_vars,ONLY: MIN_MOMENTUM, MIN_MAGNETIC
USE MOD_FiniteVolume2D_vars,ONLY: MIN_PSI
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Cons(1:nVar)
REAL,INTENT(OUT) :: Prim(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, Sx, Sy, Sz, E, Bx, By, Bz
REAL             :: S2, B2
REAL             :: Psi
!-------------------------------------------------------------------------------!

rho = Cons(1)
Sx  = Cons(2)
Sy  = Cons(3)
Sz  = Cons(4)
E   = Cons(5)
Bx  = Cons(6)
By  = Cons(7)
Bz  = Cons(8)
Psi = Cons(9)

S2 = Sx**2 + Sy**2 + Sz**2
B2 = Bx**2 + By**2 + Bz**2

IF (rho .LT. MIN_DENSITY) THEN
  rho = MIN_DENSITY
END IF
IF (E .LT. MIN_ENERGY) THEN
  E = MIN_ENERGY
END IF
IF (ABS(Sx) .LT. MIN_MOMENTUM) THEN
  Sx = 0.0
END IF
IF (ABS(Sy) .LT. MIN_MOMENTUM) THEN
  Sy = 0.0
END IF
IF (ABS(Sz) .LT. MIN_MOMENTUM) THEN
  Sz = 0.0
END IF
IF (ABS(Bx) .LT. MIN_MAGNETIC) THEN
  Bx = 0.0
END IF
IF (ABS(By) .LT. MIN_MAGNETIC) THEN
  By = 0.0
END IF
IF (ABS(Bz) .LT. MIN_MAGNETIC) THEN
  Bz = 0.0
END IF
IF (ABS(Psi) .LT. MIN_MAGNETIC) THEN
  Psi = 0.0
END IF

Prim(1) = rho
Prim(2) = Sx/rho
Prim(3) = Sy/rho
Prim(4) = Sz/rho
Prim(5) = KappaM1*(E - 0.5*S2/rho - 0.5*B2)
Prim(6) = Bx
Prim(7) = By
Prim(8) = Bz
Prim(9) = Psi

!-------------------------------------------------------------------------------!
END SUBROUTINE ConsToPrim
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE PrimToCons(Prim, Cons)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: KappaM1
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
USE MOD_FiniteVolume2D_vars,ONLY: MIN_DENSITY, MIN_PRESSURE
USE MOD_FiniteVolume2D_vars,ONLY: MIN_SPEED, MIN_MAGNETIC
USE MOD_FiniteVolume2D_vars,ONLY: MIN_PSI
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Cons(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: v2, B2
REAL             :: Psi
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)
Psi = Prim(9)

v2 = vx**2 + vy**2 + vz**2
B2 = Bx**2 + By**2 + Bz**2

IF (rho .LT. MIN_DENSITY) THEN
  rho = MIN_DENSITY
END IF
IF (p .LT. MIN_PRESSURE) THEN
  p = MIN_PRESSURE
END IF
IF (ABS(vx) .LT. MIN_SPEED) THEN
  vx = 0.0
END IF
IF (ABS(vy) .LT. MIN_SPEED) THEN
  vy = 0.0
END IF
IF (ABS(vz) .LT. MIN_SPEED) THEN
  vz = 0.0
END IF
IF (ABS(Bx) .LT. MIN_MAGNETIC) THEN
  Bx = 0.0
END IF
IF (ABS(By) .LT. MIN_MAGNETIC) THEN
  By = 0.0
END IF
IF (ABS(Bz) .LT. MIN_MAGNETIC) THEN
  Bz = 0.0
END IF
IF (ABS(Psi) .LT. MIN_PSI) THEN
  Psi = 0.0
END IF

Cons(1) = rho
Cons(2) = rho*vx
Cons(3) = rho*vy
Cons(4) = rho*vz
Cons(5) = sKappaM1*p + 0.5*rho*v2 + 0.5*B2
Cons(6) = Bx
Cons(7) = By
Cons(8) = Bz
Cons(9) = Psi

!-------------------------------------------------------------------------------!
END SUBROUTINE PrimToCons
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE EvaluateFlux1D(Prim,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: Kappa
USE MOD_FiniteVolume2D_vars,ONLY: sKappaM1
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: Prim(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: rho, vx, vy, vz, p, Bx, By, Bz
REAL             :: v2, B2, pt
!-------------------------------------------------------------------------------!

rho = Prim(1)
vx  = Prim(2)
vy  = Prim(3)
vz  = Prim(4)
p   = Prim(5)
Bx  = Prim(6)
By  = Prim(7)
Bz  = Prim(8)

v2 = vx**2 + vy**2 + vz**2
B2 = Bx**2 + By**2 + Bz**2
pt = (p + 0.5*B2) - Bx*Bx

Flux(1) = rho*vx
Flux(2) = rho*vx*vx + pt
Flux(3) = rho*vy*vx - By*Bx
Flux(4) = rho*vz*vx - Bz*Bx
Flux(5) = (sKappaM1*p + 0.5*rho*v2 + 0.5*B2)*vx + pt*vx - Bx*By*vy - Bx*Bz*vz
Flux(6) = 0.0
Flux(7) = By*vx - Bx*vy
Flux(8) = Bz*vx - Bx*vz
Flux(9) = 0.0

!-------------------------------------------------------------------------------!
END SUBROUTINE EvaluateFlux1D
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolver(PrimL,PrimR,NormVect,TangVect,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: nGPs
USE MOD_FiniteVolume2D_vars,ONLY: nDims
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: PrimL(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: PrimR(1:nVar,1:nGPs)
REAL,INTENT(IN)  :: NormVect(1:nDims,1:nGPs)
REAL,INTENT(IN)  :: TangVect(1:nDims,1:nGPs)
REAL,INTENT(OUT) :: Flux(1:nVar,1:nGPs)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: PrimLL(1:nVar,1:nGPs), PrimRR(1:nVar,1:nGPs)
REAL             :: ConsLL(1:nVar,1:nGPs), ConsRR(1:nVar,1:nGPs)
REAL             :: Flux_GLM(2,1:nGPs)
INTEGER          :: iGP
!-------------------------------------------------------------------------------!

DO iGP=1,nGPs
  ! Rotating the vector quantities       !
  PrimLL(1,iGP) = PrimL(1,iGP)
  PrimLL(2,iGP) = NormVect(1,iGP)*PrimL(2,iGP) + NormVect(2,iGP)*PrimL(3,iGP)
  PrimLL(3,iGP) = TangVect(1,iGP)*PrimL(2,iGP) + TangVect(2,iGP)*PrimL(3,iGP)
  PrimLL(4,iGP) = PrimL(4,iGP)
  PrimLL(5,iGP) = PrimL(5,iGP)
  PrimLL(6,iGP) = NormVect(1,iGP)*PrimL(6,iGP) + NormVect(2,iGP)*PrimL(7,iGP)
  PrimLL(7,iGP) = TangVect(1,iGP)*PrimL(6,iGP) + TangVect(2,iGP)*PrimL(7,iGP)
  PrimLL(8,iGP) = PrimL(8,iGP)
  PrimLL(9,iGP) = PrimL(9,iGP)

  PrimRR(1,iGP) = PrimR(1,iGP)
  PrimRR(2,iGP) = NormVect(1,iGP)*PrimR(2,iGP) + NormVect(2,iGP)*PrimR(3,iGP)
  PrimRR(3,iGP) = TangVect(1,iGP)*PrimR(2,iGP) + TangVect(2,iGP)*PrimR(3,iGP)
  PrimRR(4,iGP) = PrimR(4,iGP)
  PrimRR(5,iGP) = PrimR(5,iGP)
  PrimRR(6,iGP) = NormVect(1,iGP)*PrimR(6,iGP) + NormVect(2,iGP)*PrimR(7,iGP)
  PrimRR(7,iGP) = TangVect(1,iGP)*PrimR(6,iGP) + TangVect(2,iGP)*PrimR(7,iGP)
  PrimRR(8,iGP) = PrimR(8,iGP)
  PrimRR(9,iGP) = PrimR(9,iGP)
  
  CALL GLMFlux(PrimLL(:,iGP),PrimRR(:,iGP),Flux_GLM(:,iGP))

  CALL PrimToCons(PrimLL(1:nVar,iGP),ConsLL(1:nVar,iGP))
  CALL PrimToCons(PrimRR(1:nVar,iGP),ConsRR(1:nVar,iGP))  

  CALL RiemannSolverByRusanov(&
    ConsLL(1:nVar,iGP),ConsRR(1:nVar,iGP),&
    PrimLL(1:nVar,iGP),PrimRR(1:nVar,iGP),Flux(1:nVar,iGP))

  Flux(6,iGP) = Flux(6,iGP) + Flux_GLM(1,iGP)
  Flux(9,iGP) = Flux(9,iGP) + Flux_GLM(2,iGP)

  ! Rotating back the vector quantities  !
  Flux(2:3,iGP) = NormVect(1:nDims,iGP)*Flux(2,iGP) &
                + TangVect(1:nDims,iGP)*Flux(3,iGP)
  Flux(6:7,iGP) = NormVect(1:nDims,iGP)*Flux(6,iGP) &
                + TangVect(1:nDims,iGP)*Flux(7,iGP)
END DO

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolver
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE RiemannSolverByRusanov(ConsL,ConsR,PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(IN)  :: ConsL(1:nVar), ConsR(1:nVar)
REAL,INTENT(IN)  :: PrimL(1:nVar), PrimR(1:nVar)
REAL,INTENT(OUT) :: Flux(1:nVar)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL             :: FluxL(1:nVar), FluxR(1:nVar)
REAL             :: LambdaMax, fastestL, fastestR
!-------------------------------------------------------------------------------!

CALL EvaluateFlux1D(PrimL,FluxL)
CALL EvaluateFlux1D(PrimR,FluxR)
CALL WaveSpeeds1D(PrimL,fastest=fastestL)
CALL WaveSpeeds1D(PrimR,fastest=fastestR)

LambdaMax = MAX(ABS(fastestL),ABS(fastestR))

Flux = 0.5*((FluxL + FluxR) - LambdaMax*(ConsR - ConsL))

!-------------------------------------------------------------------------------!
END SUBROUTINE RiemannSolverByRusanov
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE GLMFlux(PrimL,PrimR,Flux)
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: GLM_ch
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
REAL,INTENT(INOUT) :: PrimL(1:nVar)
REAL,INTENT(INOUT) :: PrimR(1:nVar)
REAL,INTENT(OUT)   :: Flux(1:2)
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!
REAL               :: Bx_L, Bx_R, Psi_L, Psi_R
REAL               :: Bx_m, Psi_m
!-------------------------------------------------------------------------------!

Bx_L  = PrimL(6)
Bx_R  = PrimR(6)
Psi_L = PrimL(9)
Psi_R = PrimR(9)

Bx_m  = PrimL(6) + (0.5*(Bx_R-Bx_L)   - 0.5/GLM_ch*(Psi_R-Psi_L))
Psi_m = PrimL(9) + (0.5*(Psi_R-Psi_L) - 0.5*GLM_ch*(Bx_R-Bx_L))

PrimL(6) = Bx_m
PrimL(9) = Psi_m

PrimR(6) = Bx_m
PrimR(9) = Psi_m

Flux(1) = Psi_m
Flux(2) = GLM_ch*GLM_ch*Bx_m

!-------------------------------------------------------------------------------!
END SUBROUTINE GLMFlux
!===============================================================================!
!
!
!
!===============================================================================!
SUBROUTINE MixedGLM
!-------------------------------------------------------------------------------!
! Divergence Cleaning with mixed GLM
!-------------------------------------------------------------------------------!
USE MOD_FiniteVolume2D_vars,ONLY: U
USE MOD_FiniteVolume2D_vars,ONLY: CFL
USE MOD_FiniteVolume2D_vars,ONLY: nVar
USE MOD_FiniteVolume2D_vars,ONLY: GLM_alpha
!-------------------------------------------------------------------------------!
IMPLICIT NONE
!-------------------------------------------------------------------------------!
! >> FORMAL ARGUMENTS                                                           !
!-------------------------------------------------------------------------------!
! >> LOCAL VARIABLES                                                            !
!-------------------------------------------------------------------------------!

U(nVar,:,:) = U(nVar,:,:)*EXP(-CFL*GLM_alpha)

!===============================================================================!
END SUBROUTINE MixedGLM
!-------------------------------------------------------------------------------!
!
!
!
!===============================================================================!
END MODULE MOD_Equation
!-------------------------------------------------------------------------------!
