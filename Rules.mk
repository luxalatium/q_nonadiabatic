
#------------------------------------

# QTR source code

OBJECTS += Qtr.o Error.o Main.o InputFile.o Job.o Job_Nonadiabatic1d.o Job_Nonadiabatic2d.o Job_Nonadiabatic3d.o Job_Nonadiabatic4d.o Job_Nonadiabatic5d.o Log.o Parameters.o RandNum.o Nonadiabatic1d.o Nonadiabatic2d.o Nonadiabatic3d.o Nonadiabatic4d.o Nonadiabatic5d.o

TEMPOBJ := $(OBJECTS)
DEPOBJECTS := $(addprefix ../,$(TEMPOBJ))
DEPLIBS := $(addprefix ../,$(LIBS))

#------------------------------------

# Build rules

all: $(POTDIRS) $(FPOTDIRS) qtr
        #@echo
        #@echo "QTR Compilation"

qtr: $(OBJECTS) $(LIBS)
	$(CXX) -o $(TARGET_NAME) $^ $(LDFLAGS)

$(LIBS):
	$(MAKE) -C $@

$LIBS: $(POTDIRS) $(FPOTDIRS)

$(POTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(AR)" RANLIB="$(RANLIB)" CXXFLAGS="$(CXXFLAGS)"

$(FPOTDIRS):
	$(MAKE) -C $@ CC="$(CC)" CXX="$(CXX)" LD="$(LD)" AR="$(FAR)" FC="$(FC)" FFLAGS="$(FFLAGS)" RANLIB="$(RANLIB)" CXXFLAGS="$(CXXFLAGS)"

clean:
	rm -f $(OBJECTS) $(DEPENDS) qtr

clean-all: clean

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $<

DEPENDS= $(wildcard *.d)
-include $(DEPENDS)

.PHONY : all $(POTDIRS) $(FPOTDIRS) clean clean-all version.h
# DO NOT DELETE
