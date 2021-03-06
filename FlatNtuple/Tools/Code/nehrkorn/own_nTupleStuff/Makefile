#################################################
### Defining Compiling variables              ###
#################################################



CXX           = g++ -Wno-deprecated
#-m32  -Wno-deprecated
LD            = g++ -Wno-deprecated
#-m32  -Wno-deprecated 
INSTALL	      = /usr/bin/install
SHELL = /bin/sh

####################################
###                              ###
####################################

SHAREDLIBFLAGS =-L./CommonUtils/lib -lTauNtuple \
		-L./TauSpiner/lhapdf/lib -lLHAPDF \
		-L./TauSpiner/TAUOLA/lib -lTauolaCxxInterface -lTauolaFortran -lTauolaTauSpinner \
		-L./TauSpiner/hepmc/install/lib -lHepMC 

SHAREDCXXFLAGS =-I./CommonUtils \
		-I./TauSpiner/hepmc/install/include \
		-I./TauSpiner/TAUOLA/include/Tauola \
                -I./TauSpiner/TAUOLA/include/TauSpinner \
		-I./TauSpiner/lhapdf/include/ \
		-I./TauSpiner/pythia8135/include

SUBDIRS = CommonUtils/ TauSpiner/



CXXFLAGS      = -I$(ROOTSYS)/include -c -I/include $(SHAREDCXXFLAGS) -c -I./ -c
LDFLAGS       =  
ROOTLIBS      =  -L$(ROOTSYS)/lib -L/usr/lib/ -L/lib/i686/  -lCore -lCint -lHist -lGraf  -lGraf3d -lGpad -lTree -lRint -lReflexDict -lReflex -lPostscript -lMatrix -lPhysics -lGui -LObj -lThread -rdynamic -Wl,--rpath $(ROOTSYS)/lib

# -ldl -lpthread

DEFS  	      = -DSTANDALONE=1
LIBS          = $(ROOTLIBS) $(SHAREDLIBFLAGS)

TARGETS       =	Ntuple_Controller \
		NtupleReader \
		Parameters \
		Selection_Base \
		Selection_Factory \
		Selection \
		Tables \
		Plots \
		DoubleEventRemoval \
		PDG_Var \
		Tools \
		SkimConfig \
		HistoConfig \
		Example \
		TauSpinerInterface \
		ZtoEMu \
		ZtoEMu_Skim

HDRS          = $(TARGETS:=.h) 		
SRCS          = $(TARGETS:=.cxx)  Analysis.cxx
OBJS          = $(SRCS:.cxx=.o)

PROGRAM       = Analysis.exe

$(PROGRAM):     $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@echo $(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@$(LD) $(LDFLAGS) i386_linux/*.o $(LIBS) -o $(PROGRAM)
	@echo "done"

VPATH = utilities:i386_linux
vpath %.cxx inugent
vpath %.h inugent
vpath %.cxx cherepanov
vpath %.h cherepanov
vpath %.cxx nehrkorn
vpath %.h nehrkorn

#	@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o ../$(PROGRAM)

##########################################################
###                                                    ###
##########################################################

$(OBJS): %.o : %.cxx
	$(CXX) $(CXXFLAGS) $(DEFS) $< -o i386_linux/$@ 

.PHONY: clean cleanall all install sharedlib 

install: Analysis.exe


clean:
	@rm i386_linux/*.o
	@rm Analysis.exe

cleanall:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone clean && cd ../../; \
	done
	@rm i386_linux/*.o
	@rm Analysis.exe

all: sharedlib install

sharedlib:
	@for subdir in $(SUBDIRS); do \
		echo "Making all in $$subdir"; \
		cd $$subdir && $(MAKE) -f Makefile.standalone && cd ../ && pwd; \
	done
	$(install)

