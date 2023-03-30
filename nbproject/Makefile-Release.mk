#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=Cygwin_1-Windows
CND_DLIB_EXT=dll
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/DEM.o \
	${OBJECTDIR}/IO.o \
	${OBJECTDIR}/LB.o \
	${OBJECTDIR}/elmt.o \
	${OBJECTDIR}/hybird.o \
	${OBJECTDIR}/myvector.o \
	${OBJECTDIR}/node.o \
	${OBJECTDIR}/utils.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=-fopenmp -w -O3 -std=c++11
CXXFLAGS=-fopenmp -w -O3 -std=c++11

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hybird.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hybird.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	g++ -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hybird.exe ${OBJECTFILES} ${LDLIBSOPTIONS} -fopenmp

${OBJECTDIR}/DEM.o: DEM.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/DEM.o DEM.cpp

${OBJECTDIR}/IO.o: IO.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/IO.o IO.cpp

${OBJECTDIR}/LB.o: LB.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/LB.o LB.cpp

${OBJECTDIR}/elmt.o: elmt.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/elmt.o elmt.cpp

${OBJECTDIR}/hybird.o: hybird.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/hybird.o hybird.cpp

${OBJECTDIR}/myvector.o: myvector.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/myvector.o myvector.cpp

${OBJECTDIR}/node.o: node.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/node.o node.cpp

${OBJECTDIR}/utils.o: utils.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O3 -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/utils.o utils.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
