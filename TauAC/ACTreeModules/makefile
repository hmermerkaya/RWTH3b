TARGETS		= ACAnalyzer ACEvent ACTreeReader ACExampleAnalyzer
OBJ			= $(addsuffix .dir, $(TARGETS))

%.dir:
	cd ./$(subst .dir,,$@); make

all: $(OBJ)
	@echo " --> Done!"

.SILENT: 
.PHONY: clean

clean:
	@echo " --> Cleaning shared objects, symbolic links, and executables."
	-for d in $(TARGETS); do (cd $$d; make clean); done
	@echo " --> Done!"
