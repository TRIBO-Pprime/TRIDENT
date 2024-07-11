include ../../makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src
LIBDIR = lib

EXE = main

ifeq ($(FORT),x86_64-w64-mingw32-gfortran-win32)
	LIBDIR = libwin64
	EXE = main.exe
endif

ALL_LIBS = $(addprefix $(fftw3_dir)/$(LIBDIR)/, *.a)

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ       = $(SRC:.$(EXT)=.o)
ALL_OBJS  = $(addprefix $(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(algen_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(chole_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(digis_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(fftw3_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(files_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(gplot_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(intpl_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(least_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(minim_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(qsort_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(splin_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(tchev_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(utils_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(abbot_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(aniso_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(asfc2_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(deriv_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(filtr_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(morph_dir)/$(OBJDIR)/, mod_*.o)
ALL_OBJS += $(addprefix $(stats_dir)/$(OBJDIR)/, mod_*.o)

ALL_MODS  = -I$(MODDIR)
ALL_MODS += -I$(algen_dir)/$(MODDIR)
ALL_MODS += -I$(chole_dir)/$(MODDIR)
ALL_MODS += -I$(digis_dir)/$(MODDIR)
ALL_MODS += -I$(fftw3_dir)/$(MODDIR)
ALL_MODS += -I$(files_dir)/$(MODDIR)
ALL_MODS += -I$(gplot_dir)/$(MODDIR)
ALL_MODS += -I$(intpl_dir)/$(MODDIR)
ALL_MODS += -I$(least_dir)/$(MODDIR)
ALL_MODS += -I$(minim_dir)/$(MODDIR)
ALL_MODS += -I$(qsort_dir)/$(MODDIR)
ALL_MODS += -I$(splin_dir)/$(MODDIR)
ALL_MODS += -I$(tchev_dir)/$(MODDIR)
ALL_MODS += -I$(utils_dir)/$(MODDIR)
ALL_MODS += -I$(abbot_dir)/$(MODDIR)
ALL_MODS += -I$(aniso_dir)/$(MODDIR)
ALL_MODS += -I$(asfc2_dir)/$(MODDIR)
ALL_MODS += -I$(deriv_dir)/$(MODDIR)
ALL_MODS += -I$(filtr_dir)/$(MODDIR)
ALL_MODS += -I$(morph_dir)/$(MODDIR)
ALL_MODS += -I$(stats_dir)/$(MODDIR)

CFLAGS  = $(ALL_MODS) -fopenmp -fPIC -fdec
CFLAGS += -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -finit-real=nan -ffpe-summary=none

LFLAGS  = $(ALL_LIBS)

ifneq ($(FORT),x86_64-w64-mingw32-gfortran-win32)
	LFLAGS += -lpthread -lm
	LFLAGS += $(ALL_LIBS) /lib/gcc/x86_64-linux-gnu/12/libgomp.a /lib/gcc/x86_64-linux-gnu/12/libquadmath.a
	LFLAGS += -static-libgfortran -static-libgcc
endif

#LFLAGS += -lm -lgomp -static-libgfortran -static-libgcc -lpthread

ifeq ($(FORT),x86_64-w64-mingw32-gfortran-win32)
	LFLAGS += -lm -lgomp -static -lpthread -static
endif

ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
	CFLAGS += -O2
endif

ifneq ('$(GPROF)','')
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

%.o:	%.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	@find . -maxdepth 1 -name '*.mod*' -type f -print0 | xargs -0r mv -t ./$(MODDIR)

$(EXE):	$(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)
	@rm $(OBJDIR)/main.o

mod_script.o :
main.o : mod_script.o

#--------------------------------------------------
.PHONY: clean debug gprof all last

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f $(EXE)

debug:
	make "DEBUG=TRUE"

gprof:
	make "GPROF=TRUE"

all:
	(make clean ; make)

