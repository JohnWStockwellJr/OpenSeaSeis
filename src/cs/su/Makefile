MODULE_NAME = su

MODDIR    = $(SRCDIR)/cs/su
MODULE    = $(MODDIR)/mod_$(MODULE_NAME).cc
OBJS      = $(OBJDIR)/mod_$(MODULE_NAME).o $(OBJDIR)/csSUTraceManager.o $(OBJDIR)/csSUArguments.o $(OBJDIR)/csSUGetPars.o
LIB_v1.0  = libmod_$(MODULE_NAME).so.1.0
LIB_SU    = libsumodules.so
LIB_SU_PAR = libcseis_supar.so

INCS =  -I"$(SRCDIR)/cs/geolib"  -I"$(SRCDIR)/cs/system"  -I"$(SRCDIR)/cs/segy" -I"$(MODDIR)" -I"$(CWPROOT)/include" -I"$(CWPROOT)/src/Complex/include"

# -Wno-write-strings: Suppress warnings about deprecated string assignments
ALL_FLAGS = $(INCS) $(COMMON_FLAGS) -fPIC -g -Wno-write-strings

default: $(OBJS) $(LIBDIR)/$(LIB_v1.0)

clean:
	rm -f $(OBJS)

bleach: clean
	rm -f $(LIBDIR)/$(LIB_v1.0)

$(OBJDIR)/mod_$(MODULE_NAME).o: $(MODULE) $(MODDIR)/su_modules.h
	$(CPP) -c $(ALL_FLAGS) $(MODULE) -o $(OBJDIR)/mod_$(MODULE_NAME).o

$(OBJDIR)/csSUTraceManager.o: $(MODDIR)/csSUTraceManager.h $(MODDIR)/csSUTraceManager.cc
	$(CPP) -c $(ALL_FLAGS) $(MODDIR)/csSUTraceManager.cc -o $(OBJDIR)/csSUTraceManager.o

$(OBJDIR)/csSUArguments.o: $(MODDIR)/csSUArguments.h $(MODDIR)/csSUArguments.cc
	$(CPP) -c $(ALL_FLAGS) $(MODDIR)/csSUArguments.cc -o $(OBJDIR)/csSUArguments.o

$(OBJDIR)/csSUGetPars.o: $(MODDIR)/csSUGetPars.h $(MODDIR)/csSUGetPars.cc
	$(CPP) -c $(ALL_FLAGS) $(MODDIR)/csSUGetPars.cc -o $(OBJDIR)/csSUGetPars.o

$(LIBDIR)/$(LIB_v1.0): $(OBJS) $(OBJS_SU) $(LIBDIR)/$(LIB_SU) $(LIBDIR)/$(LIB_SU_PAR)
	$(CPP) -shared -Wl,-soname,$(LIB_v1.0) -o $(LIBDIR)/$(LIB_v1.0) $(OBJS) $(CWPROOT)/lib/libcwp.a $(CWPROOT)/lib/libpar.a -L$(LIBDIR) -L$(CWPROOT)/lib -lc -lgeolib -lcseis_system -lsegy -lsu -lm -lpthread -lsumodules
