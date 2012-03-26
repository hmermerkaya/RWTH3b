ARCH			:= $(shell $(ROOTSYS)/bin/root-config --arch)

ifeq ($(ARCH),linuxx8664gcc)
CXX				:= g++
LIBRARYDIR		:= $(CMSSW_BASE)/lib/$(SCRAM_ARCH)
endif

ROOTCXXFLAGS	:= $(shell $(ROOTSYS)/bin/root-config --cflags) 
ROOTLIBS		:= $(shell $(ROOTSYS)/bin/root-config --libs) 
BOOST_VERSION	:= 1.44.0-cms3
BOOST_INCLUDE	:= $(CMS_PATH)/$(SCRAM_ARCH)/external/boost/$(BOOST_VERSION)/include
BOOST_LIBRARY	:= $(CMS_PATH)/$(SCRAM_ARCH)/external/boost/$(BOOST_VERSION)/lib
BOOST_LIBNAMES	:= 

IDIR			:= $(CURDIR)/interface
SDIR			:= $(CURDIR)/source
ODIR			:= ../build
LDIR			:= ../lib

$(shell [ -d '$(ODIR)' ] || mkdir -p $(ODIR))
$(shell [ -d '$(LDIR)' ] || mkdir -p $(LDIR))

CXXFLAGS		:= -O -Wall -fPIC $(ROOTCXXFLAGS) -I$(BOOST_INCLUDE) -I$(IDIR)

TARGET			:= ACHarvester
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
