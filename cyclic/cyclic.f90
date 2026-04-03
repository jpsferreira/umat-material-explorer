!C234567890123456789012345678901234567890123456789012345678901234567890
PROGRAM TEST
use,intrinsic :: ISO_Fortran_env
INCLUDE 'aba_param.inc'
INCLUDE 'param_umat.inc'

!C     ADD COMMON BLOCKS HERE IF NEEDED ()
!C      COMMON /KBLOCK/KBLOCK

PARAMETER(NTENS = 6, NSTATEV = NSDV, NPROPS = 13, NDI=3, NSHR=3)
PARAMETER(NOEL = 1, NPT = 8)

CHARACTER*8 CMNAME
DIMENSION STRESS(NTENS),STATEV(NSTATEV),STATEVP(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

CHARACTER*8 filename
integer :: un = 200
!
real*8,ALLOCATABLE::curvepoints(:,:)
real*8::maxstresstime,maxstress,maxgammatime,maxgamma,minstresstime,minstress
real*8:: maxstretch, maxstretchtime
real*8:: delta,lossmodulus

DO I1=1,NTENS
    DO J1=1,NTENS
        DDSDDE(I1,J1)=0.0D0
    ENDDO
    STRESS(I1)=0.0D0
ENDDO
!!!
time(1)=ZERO
time(2)=ZERO
call UEXTERNALDB(0,0,time,ZERO,0,0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MATERIAL PROPERTIES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k PENALTY PARAMETER
PROPS(1)=2.d0/1000.000d0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ELASTIC PARAM
PROPS(2)=1.00d0  ! C10
PROPS(3)=0.00d0 ! C01
PROPS(4)=1.00d0 !K1
PROPS(5)=1.d0 !K2
PROPS(6)=0.1d0 !kdisp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!viscous parameters - maxwell
! v - number of dashpots
PROPS(7)=1
!tau1 %
PROPS(8)=2.0d0
!teta1
PROPS(9)=0.835d0
!tau2 %
PROPS(10)=1.2d0
!teta2
PROPS(11)=7.0d0
!tau3 %
PROPS(12)=12.d0
!teta3
PROPS(13)=2.0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
STATEV=ZERO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
erf=ZERO
RHO=ZERO
PI=FOUR*ATAN(ONE)
!
!################################################################################################!
NSTEPS=100 !NUMBER OF POINTS PER CYCLE
!
!################################################################################################!
!CYLCIC RELATED VARIABLES
NTESTSF=61 !NUMBER OF TESTS
NTESTSA=31 !NUMBER OF TESTS
NTESTSAS=21 !NUMBER OF TESTS
!
NCYCLES=3 !NUMBER OF CYCLES PER TEST
!
ALLOCATE(CURVEPOINTS(NTESTSF,3))
!
!FREQUENCY SWEEP
FREQMIN=0.003d0
DFREQ=0.01d0
!
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS - FREQ VALUES
DO I1=1,NTESTSF
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,1)=FREQMIN*10.0d0**COUNTER
ENDDO
!
!AMPLITUDE SWEEP
FREQ0 = 1.d0
AMP_STRETCH=0.05D0
PRE_STRETCH=1.1D0

AMP_MIN = 0.01D0
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS - AMP STRETCH VALUES
DO I1=1,NTESTSA
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,2)=AMP_MIN*10.0d0**COUNTER
ENDDO
COUNTER=-0.05d0
!UNIFORM DISTRIBUTITED LOG-SCALE POINTS - AMP GAMMA VALUES
DO I1=1,NTESTSAS
  COUNTER=COUNTER+0.05d0
  CURVEPOINTS(I1,3)=AMP_MIN*10.0d0**COUNTER
ENDDO

AMP_GAMMA=0.3D0
PRE_GAMMA=0.1D0

! Create output directories
CALL SYSTEM('mkdir -p stress_curves/freq_sweep/uniaxial')
CALL SYSTEM('mkdir -p stress_curves/freq_sweep/equibiaxial')
CALL SYSTEM('mkdir -p stress_curves/freq_sweep/shear')
CALL SYSTEM('mkdir -p stress_curves/freq_sweep/sshear')
CALL SYSTEM('mkdir -p stress_curves/amp_sweep/uniaxial')
CALL SYSTEM('mkdir -p stress_curves/amp_sweep/equibiaxial')
CALL SYSTEM('mkdir -p stress_curves/amp_sweep/shear')
CALL SYSTEM('mkdir -p stress_curves/amp_sweep/sshear')

!########################################################################################!
!############################## |||   FREQUENCY SWEEP   ||| #############################!
!########################################################################################!
!
!################################## |||   UNIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=150, FILE='stress_curves/freq_sweep/uniaxial.out', STATUS='UNKNOWN')
rewind(150)
!
DO I1=1,NTESTSF
  CALL RUN_CYCLIC_TEST(I1, CURVEPOINTS(I1,1), NSTEPS, NCYCLES, &
    1, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/freq_sweep/uniaxial/')
  write(150,*) CURVEPOINTS(I1,1),storagemodulus,lossmodulus,tandelta
ENDDO
close(150)

!################################## |||   BIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=151, FILE='stress_curves/freq_sweep/biaxial.out', STATUS='UNKNOWN')
rewind(151)
!
DO I1=1,NTESTSF
  CALL RUN_CYCLIC_TEST(I1, CURVEPOINTS(I1,1), NSTEPS, NCYCLES, &
    2, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/freq_sweep/equibiaxial/')
  write(151,*) CURVEPOINTS(I1,1),storagemodulus,tandelta,lossmodulus
ENDDO
close(151)

!################################## |||   SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=152, FILE='stress_curves/freq_sweep/shear.out', STATUS='UNKNOWN')
rewind(152)
!
DO I1=1,NTESTSF
  CALL RUN_CYCLIC_TEST(I1, CURVEPOINTS(I1,1), NSTEPS, NCYCLES, &
    3, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/freq_sweep/shear/')
  write(152,*) CURVEPOINTS(I1,1),storagemodulus,tandelta,lossmodulus
ENDDO
close(152)

!################################## |||   SIMPLE SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=153, FILE='stress_curves/freq_sweep/sshear.out', STATUS='UNKNOWN')
rewind(153)
!
DO I1=1,NTESTSF
  CALL RUN_CYCLIC_TEST(I1, CURVEPOINTS(I1,1), NSTEPS, NCYCLES, &
    4, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/freq_sweep/sshear/')
  write(153,*) CURVEPOINTS(I1,1),storagemodulus,tandelta,lossmodulus
ENDDO
close(153)


!########################################################################################!
!############################## |||   AMPLITUDE SWEEP   ||| #############################!
!########################################################################################!

!################################## |||   UNIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=150, FILE='stress_curves/amp_sweep/uniaxial.out', STATUS='UNKNOWN')
rewind(150)
!
DO I1=1,NTESTSA
  AMP_STRETCH=CURVEPOINTS(I1,2)
  CALL RUN_CYCLIC_TEST(I1, FREQ0, NSTEPS, NCYCLES, &
    1, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/amp_sweep/uniaxial/')
  write(150,*) AMP_STRETCH,storagemodulus,lossmodulus,tandelta
ENDDO
close(150)

!################################## |||   BIAXIAL   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=150, FILE='stress_curves/amp_sweep/biaxial.out', STATUS='UNKNOWN')
rewind(150)
!
DO I1=1,NTESTSA
  AMP_STRETCH=CURVEPOINTS(I1,2)
  CALL RUN_CYCLIC_TEST(I1, FREQ0, NSTEPS, NCYCLES, &
    2, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/amp_sweep/equibiaxial/')
  write(150,*) AMP_STRETCH,storagemodulus,lossmodulus,tandelta
ENDDO
close(150)

!################################## |||   SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=152, FILE='stress_curves/amp_sweep/shear.out', STATUS='UNKNOWN')
rewind(152)
!
DO I1=1,NTESTSAS
  AMP_GAMMA=CURVEPOINTS(I1,3)
  CALL RUN_CYCLIC_TEST(I1, FREQ0, NSTEPS, NCYCLES, &
    3, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/amp_sweep/shear/')
  write(152,*) AMP_GAMMA,storagemodulus,tandelta,lossmodulus
ENDDO
close(152)

!################################## |||   SIMPLE SHEAR   ||| ################################!
CALL RESETDFGRD(DFGRD1,NDI)
OPEN (UNIT=153, FILE='stress_curves/amp_sweep/sshear.out', STATUS='UNKNOWN')
rewind(153)
!
DO I1=1,NTESTSAS
  AMP_GAMMA=CURVEPOINTS(I1,3)
  CALL RUN_CYCLIC_TEST(I1, FREQ0, NSTEPS, NCYCLES, &
    4, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
    STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
    STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
    NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
    CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
    maxstress, minstress, maxstresstime, minstresstime, &
    maxstretch, maxstretchtime, maxgamma, maxgammatime, &
    storagemodulus, lossmodulus, tandelta, delta, un, &
    'stress_curves/amp_sweep/sshear/')
  write(153,*) AMP_GAMMA,storagemodulus,tandelta,lossmodulus
ENDDO
close(153)
!
! !################################################################################################!
CALL SYSTEM('gnuplot -p plots/data_cyclicfreq.plt')
! !################################################################################################!
!
! !################################################################################################!
CALL SYSTEM('gnuplot -p plots/data_cyclicamp.plt')
! !################################################################################################!
END PROGRAM

!========================================================================================!
! Subroutine to run a single cyclic test (one frequency/amplitude point)
!
! deform_type: 1=uniaxial, 2=biaxial, 3=pure shear, 4=simple shear
!========================================================================================!
SUBROUTINE RUN_CYCLIC_TEST(ITEST, FREQ, NSTEPS, NCYCLES, &
  deform_type, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA, &
  STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
  STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
  NDI, NSHR, NTENS, NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
  CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KINC, PI, &
  maxstress, minstress, maxstresstime, minstresstime, &
  maxstretch, maxstretchtime, maxgamma, maxgammatime, &
  storagemodulus, lossmodulus, tandelta, delta, un, outpath)

INCLUDE 'aba_param.inc'
INCLUDE 'param_umat.inc'

INTEGER ITEST, NSTEPS, NCYCLES, deform_type, un
CHARACTER*8 CMNAME
CHARACTER*(*) outpath
CHARACTER*8 filename
INTEGER NDI, NSHR, NTENS, NSTATEV, NPROPS, NOEL, NPT, LAYER, KSPT, KINC
DIMENSION STRESS(NTENS),STATEV(NSTATEV),DDSDDE(NTENS,NTENS),DDSDDT(NTENS),      &
DRPLDE(NTENS),STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),            &
PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
real*8 FREQ, AMP_STRETCH, PRE_STRETCH, AMP_GAMMA, PRE_GAMMA
real*8 maxstress, minstress, maxstresstime, minstresstime
real*8 maxstretch, maxstretchtime, maxgamma, maxgammatime
real*8 storagemodulus, lossmodulus, tandelta, delta
real*8 F, CYCLETIME, STRETCH, GAMMA, stressamp, amp_val, PI
real*8 SSE, SPD, SCD, RPL, DRPLDT, DTIME, TEMP, DTEMP, PNEWDT, CELENT
INTEGER stress_idx

! Determine which stress component to track
IF (deform_type .LE. 2) THEN
  stress_idx = 1
ELSE
  stress_idx = 4
ENDIF

TIME(1) = ZERO
F = FREQ * TWO * PI

WRITE(filename, fmt='(i0,a)') ITEST, '.out'
un = un + 2
OPEN(unit=un, file=outpath//filename, status='UNKNOWN')

CYCLETIME = ONE / FREQ
DTIME = CYCLETIME / NSTEPS

! Initialize tracking variables
IF (deform_type .LE. 2) THEN
  maxstress = -1.0d30
  maxstresstime = -1.0d30
  minstress = 1.0d30
  minstresstime = 1.0d30
  maxstretch = -1.0d30
  maxstretchtime = -1.0d30
ELSE
  maxstress = -1.0d30
  maxstresstime = -1.0d30
  minstress = 1.0d30
  minstresstime = 1.0d30
  maxgamma = -1.0d30
  maxgammatime = -1.0d30
ENDIF

DO KK = 1, NCYCLES
  DO KSTEP = 1, NSTEPS
    ! Apply deformation
    SELECT CASE (deform_type)
    CASE (1) ! Uniaxial
      STRETCH = AMP_STRETCH * DSIN(F * TIME(1)) + PRE_STRETCH
      DFGRD1(1,1) = STRETCH
      DFGRD1(2,2) = ONE / SQRT(STRETCH)
      DFGRD1(3,3) = ONE / SQRT(STRETCH)
    CASE (2) ! Biaxial
      STRETCH = AMP_STRETCH * DSIN(F * TIME(1)) + PRE_STRETCH
      DFGRD1(1,1) = STRETCH
      DFGRD1(2,2) = STRETCH
      DFGRD1(3,3) = ONE / (STRETCH * STRETCH)
    CASE (3) ! Pure shear
      GAMMA = AMP_GAMMA * DSIN(F * TIME(1))
      DFGRD1(1,2) = GAMMA
      DFGRD1(2,1) = GAMMA
    CASE (4) ! Simple shear
      GAMMA = AMP_GAMMA * DSIN(F * TIME(1)) + PRE_GAMMA
      DFGRD1(1,2) = GAMMA
    END SELECT

    CALL UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, STRAN, &
      DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATEV, PROPS, &
      NPROPS, COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)

    ! Track extrema on last cycle
    IF (KK .EQ. NCYCLES) THEN
      IF (STRESS(stress_idx) .GT. maxstress) THEN
        maxstress = STRESS(stress_idx)
        maxstresstime = TIME(1)
      ENDIF
      IF (STRESS(stress_idx) .LT. minstress) THEN
        minstress = STRESS(stress_idx)
        minstresstime = TIME(1)
      ENDIF
      IF (deform_type .LE. 2) THEN
        IF (STRETCH .GT. maxstretch) THEN
          maxstretch = STRETCH
          maxstretchtime = TIME(1)
        ENDIF
      ELSE
        IF (GAMMA .GT. maxgamma) THEN
          maxgamma = GAMMA
          maxgammatime = TIME(1)
        ENDIF
      ENDIF
    ENDIF

    TIME(1) = TIME(1) + DTIME
    IF (deform_type .LE. 2) THEN
      write(un, *) TIME(1), STRETCH, STRESS(stress_idx)
    ELSE
      write(un, *) TIME(1), GAMMA, STRESS(stress_idx)
    ENDIF
  ENDDO
ENDDO

! Compute viscoelastic moduli
stressamp = maxstress - minstress
IF (deform_type .LE. 2) THEN
  delta = maxstretchtime - maxstresstime
  amp_val = AMP_STRETCH
ELSE
  delta = maxgammatime - maxstresstime
  amp_val = AMP_GAMMA
ENDIF
delta = delta * TWO * PI / CYCLETIME

storagemodulus = (maxstress / amp_val) * DCOS(delta)
lossmodulus = (maxstress / amp_val) * DSIN(delta)
tandelta = lossmodulus / storagemodulus

close(unit=un)

END SUBROUTINE RUN_CYCLIC_TEST
