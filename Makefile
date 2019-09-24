# When we use ifort,
# -heap-arrays option is needed for built in transpose function with
#  large dimension matrix
#--------------------------------------------------
TARGET=NdSpline
MODDIR = main
HOST=$(strip $(Host))
DEBUG_MODE=on
HOST=other

OS = Linux
ifneq (,$(findstring arwin,$(shell uname)))
  OS = OSX
endif
$(info Debug mode is $(DEBUG_MODE).)

#--------------------------------------------------
# Default Parameters
#--------------------------------------------------
ifeq ($(strip $(HOST)),other)
  FDEP=makedepf90
  FC=gfortran
  ifeq ($(MPI), on)
    FC=mpif90 -DMPI
  endif
  LFLAGS+= -I/usr/local/include -L/usr/local/lib
  LFLAGS+= -lblas -llapack -lgsl -lz
  FFLAGS=-O3
  FFLAGS+= -fopenmp
  ifeq ($(DEBUG_MODE),on)
    DFLAGS+=-Wall -pedantic -fbounds-check -O -Wuninitialized -fbacktrace
  endif
endif


ifeq ($(DEBUG_MODE),on)
  DFLAGS+=
endif



#--------------------------------------------------
# Source Files
#--------------------------------------------------

SRCDIR = src
OBJDIR = obj
DEPDIR = .

SRCS=
OBJS=
MODS=

SRCF95:=$(wildcard $(SRCDIR)/*.F90)
OBJF95:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95))))
SRCS= $(SRCF95)
OBJS= $(OBJF95)

SRCDIR_MAIN = main

SRCF95_MAIN:=$(wildcard $(SRCDIR_MAIN)/*.F90)
OBJF95_MAIN:=$(addprefix $(OBJDIR)/, $(patsubst %F90, %o, $(notdir $(SRCF95_MAIN))))
SRCS_MAIN= $(SRCF95_MAIN)
OBJS_MAIN= $(OBJF95_MAIN)

SRCS_ALL = $(SRCS) $(SRCS_MAIN)
OBJS_ALL = $(OBJS) $(OBJS_MAIN)

MODOUT=
ifeq ($(strip $(HOST)),other)
  MODOUT=-J$(MODDIR)
endif

#--------------------------------------------------
# Rules
#--------------------------------------------------
all: dirs $(TARGET)
$(TARGET): $(OBJS_ALL)
	$(FC) $(FFLAGS) $(DFLAGS) -o $(TARGET).exe $^ $(LFLAGS)

$(OBJDIR)/%.o:$(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

$(OBJDIR)/%.o:$(SRCDIR_MAIN)/%.F90
	$(FC) $(FFLAGS) $(DFLAGS) $(MODOUT) -o $@ -c $<

dep:
	$(FDEP) $(SRCS_ALL) -b $(OBJDIR)/ > $(DEPDIR)/makefile.d

dirs:
	if test -d $(OBJDIR); then \
		: ; \
	else \
		mkdir $(OBJDIR); \
	fi
	if test -d $(MODDIR); then \
		: ; \
	else \
		mkdir $(MODDIR); \
	fi

clean:
	rm -f $(TARGET).exe
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJS_ALL)

#--------------------------------------------------
-include $(wildcard $(DEPDIR)/*.d)
