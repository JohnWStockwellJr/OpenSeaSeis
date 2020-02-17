MODULE_NAME = trc_math_ens

MODULE    = $(SRCDIR)/cs/modules/$(MODULE_NAME)/mod_$(MODULE_NAME).cc
MODDIR    = $(SRCDIR)/cs/modules/$(MODULE_NAME)/
OBJS      = $(OBJDIR)/mod_$(MODULE_NAME).o
LIB_v1.0  = libmod_$(MODULE_NAME).so.1.0

INCS =  -I"$(SRCDIR)/cs/geolib"  -I"$(SRCDIR)/cs/system" -I"$(SRCDIR)/cs/segd"  -I"$(SRCDIR)/cs/io" -I"$(SRCDIR)/cs/segy" -I"$(MODDIR)"

ALL_FLAGS = $(INCS) $(COMMON_FLAGS) -fPIC

default: $(OBJS) $(LIBDIR)/$(LIB_v1.0)

clean:
	rm -f $(OBJS)

bleach: clean
	rm -f $(LIBDIR)/$(LIB_v1.0)

$(OBJDIR)/mod_$(MODULE_NAME).o: $(MODULE)
	$(CPP) -c $(ALL_FLAGS) $(MODULE) -o $(OBJDIR)/mod_$(MODULE_NAME).o

$(LIBDIR)/$(LIB_v1.0): $(OBJS)
	$(CPP) -shared -Wl,-$(SONAME),$(LIB_v1.0) -o $(LIBDIR)/$(LIB_v1.0) $(OBJS) -L$(LIBDIR) -lc -lgeolib -lcseis_system
