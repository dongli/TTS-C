PROJECT_NAME = "TTS Unit Tests"
DEBUG = yes
LIBRARY = -L$${NCARG_ROOT}/lib -lngmath
TARGET = TestMsgManager TestMeshManager TestPolarStereo TestNFWrap \
         TestVoronoi
OBJECT = MsgManager.o CommonTypes.o RunManager.o FloatingPoint.o \
         NFWrap.o SphereService.o MeshManager.o FlowManager.o \
         TracerManager.o TTS.o LagrangeToEuler.o SolidRotationTestbed.o \
         MovingVorticesTestbed.o DeformationTestbed.o

include ../Makefile.basic
