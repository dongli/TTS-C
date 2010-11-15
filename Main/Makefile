PROJECT_NAME = "TTS Test Case Drivers"
DEBUG = yes
OPTIONS = #-DLAGRANGE_TO_EULER
TARGET = DriveSolidRotation DriveMovingVortices DriveDeformation \
         GenSearchList DumpSearchList CalcError_1 CalcError_2 \
         CalcExactSolidRotation CalcExactMovingVortices
OBJECT = MsgManager.o CommonTypes.o RunManager.o FloatingPoint.o \
         NFWrap.o SphereService.o MeshManager.o FlowManager.o \
         TracerManager.o TTS.o LagrangeToEuler.o SolidRotationTestbed.o \
         MovingVorticesTestbed.o DeformationTestbed.o

include ../Makefile.basic
