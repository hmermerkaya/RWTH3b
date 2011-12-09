ARCH              = $(shell $(ROOTSYS)/bin/root-config --arch)

ifeq ($(ARCH),linuxx8664gcc)
CXX               := g++
LIBRARYDIR        := $(CMSSW_BASE)/lib/$(SCRAM_ARCH)
ACLIBS            := $(shell find $(LIBRARYDIR) -name 'libACDataFormats*' -prune -printf '%f ')
ACLIBS            := $(subst lib,-l,$(ACLIBS))
ACLIBS            := $(subst .so,,$(ACLIBS))
BOOST_VERSION     := 1.44.0-cms3
BOOST_INCLUDE     := $(CMS_PATH)/$(SCRAM_ARCH)/external/boost/$(BOOST_VERSION)/include
BOOST_LIBRARY     := $(CMS_PATH)/$(SCRAM_ARCH)/external/boost/$(BOOST_VERSION)/lib
BOOST_LIBNAMES    := 
endif
ACTREELIBS        := -lACAnalyzer -lACEvent -lACTreeReader
ACTREELIBSDIR	  := ../lib

ROOTCXXFLAGS      := $(shell $(ROOTSYS)/bin/root-config --cflags) 
ROOTLIBS          := $(shell $(ROOTSYS)/bin/root-config --libs) 

IDIR              := $(CURDIR)/interface
SDIR              := $(CURDIR)/source
ODIR              := ../build

$(shell [ -d '$(ODIR)' ] || mkdir -p $(ODIR))

CXXFLAGS          := -O -Wall $(ROOTCXXFLAGS) -I$(IDIR)

TARGET			:= ACExampleAnalyzer
_OBJ			:= $(addsuffix .o, $(TARGET))
_DEPS			:= $(addsuffix .h, $(TARGET))

DEPS              = $(patsubst %,$(IDIR)/%,$(_DEPS))
OBJ               = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: $(SDIR)/%.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET).exe: $(OBJ)
	@echo " --> $(TARGET): Building executable."
	$(CXX) $(CXXFLAGS) $(ROOTLIBS) -L$(BOOST_LIBRARY) $(BOOST_LIBNAMES) -L$(LIBRARYDIR) $(ACLIBS) -L$(ACTREELIBSDIR) $(ACTREELIBS) $^ -o $@

all: $(TARGET).exe

.PHONY: clean

clean:
	rm -Rf $(ODIR)
	rm -f *.exe

print-%:
	@echo $* = $($*)
