TARGET_NAME=qtr_serial

#-----------------------------------------------------------

#C compiler settings

CC=icc
CXX=icpc
LD=ld
AR=ar cru
RANLIB=ranlib

ifndef DEBUG
        NDEBUG=1
        OPT=-O3 -qopenmp -xCORE-AVX512 -ipo
endif

CXXFLAGS +=  -Wall -Wno-sign-compare -Wno-unused-function -Wfatal-errors $(OPT) -std=c++11
LDFLAGS += -qopenmp -xCORE-AVX512 -ipo

#-----------------------------------------------------------

#MPI settings

ifdef QTRMPI
        CXX=mpicxx
        CXXFLAGS += -DQTRMPI
        TARGET_NAME=qtr_mpi
endif

#-----------------------------------------------------------

# Finite difference order

#POT settings

ifdef NADBPOT_TN1
        CXXFLAGS += -DNADBPOT_TN1
endif

ifdef NADBPOT_TN2
        CXXFLAGS += -DNADBPOT_TN2
endif

ifdef NADBPOT_GM3
        CXXFLAGS += -DNADBPOT_GM3
endif

ifdef NADBPOT_MM4
        CXXFLAGS += -DNADBPOT_MM4
endif

ifdef NADBPOT_GM5
        CXXFLAGS += -DNADBPOT_GM5
endif

ifdef NADBPOT_GM6
        CXXFLAGS += -DNADBPOT_GM6
endif

ifdef NADBPOT_PBC
        CXXFLAGS += -DNADBPOT_PBC
endif

#-----------------------------------------------------------

#Ensures all code is statically linked on a Linux machine

UNAME := $(shell uname)
MACHINE := $(shell uname -m)

ifeq ($(UNAME), Linux)

        LINUX=1

        ifeq ($(MACHINE), x86_64)
                LINUX64=1
        else
                LINUX32=1
        endif

#MPI often doesn't mix well with static linking...

        ifndef QTRMPI
                #LDFLAGS += -static
        endif
endif

ifeq ($(UNAME), Darwin)
        OSX=1
endif

ifeq ($(UNAME), CYGWIN_NT-5.1)
        WIN32=1
endif

ifdef LINUX
        CXXFLAGS += -DLINUX
endif

ifdef OSX
        CXXFLAGS += -DOSX
endif

ifdef WIN32
        CXXFLAGS += -DWIN32
endif

include Rules.mk
