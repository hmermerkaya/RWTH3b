ARCH			:= $(shell $(ROOTSYS)/bin/root-config --arch)

ifeq ($(ARCH),linuxx8664gcc)
CXX				:= g++
LIBRARYDIR		:= $(CMSSW_BASE)/lib/$(SCRAM_ARCH)
ACLIBS			:= $(shell find $(LIBRARYDIR) -name 'libACDataFormats*' -prune -printf '%f ')
ACLIBS			:= $(subst lib,-l,$(ACLIBS))
ACLIBS			:= $(subst .so,,$(ACLIBS))
endif

ROOTCXXFLAGS	:= $(shell $(ROOTSYS)/bin/root-config --cflags) 
ROOTLIBS		:= $(shell $(ROOTSYS)/bin/root-config --libs) 

IDIR			:= $(CURDIR)/interface
SDIR			:= $(CURDIR)/source
ODIR			:= ../build
LDIR			:= ../lib

$(shell [ -d '$(ODIR)' ] || mkdir -p $(ODIR))
$(shell [ -d '$(LDIR)' ] || mkdir -p $(LDIR))

CXXFLAGS		:= -O -Wall -fPIC $(ROOTCXXFLAGS) -I$(IDIR)

TARGET			:= ACEvent
_OBJ			:= $(addsuffix .o, $(TARGET))
_DEPS			:= $(addsuffix .h, $(TARGET))

DEPS			= $(patsubst %,$(IDIR)/%,$(_DEPS))
OBJ				= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(LDIR)/lib%.so: $(ODIR)/%.o
	@echo " --> $(TARGET): Building shared library."
	$(CXX) $(CXXFLAGS) -shared $< -o $@
	ln -f -s $(CURDIR)/$@ $(LIBRARYDIR)/

all: $(LDIR)/lib$(TARGET).so

.PHONY: clean

clean:
	rm -Rf $(ODIR)
	rm -Rf $(LDIR)
	rm -f *.exe
	rm -f $(LIBRARYDIR)/lib$(TARGET).so

print-%:
	@echo $* = $($*)
