# LAMMPS multiple-machine Makefile

SHELL = /bin/bash
#.IGNORE:

# Subtargets for makefile auto
SUBTARGET_PRE_DEBUG = $(@)
OPTIONSPOSTFIX_PRE_DEBUG =
ifneq ($(origin postfix), undefined)
    ifneq (${postfix}, )
        ifneq (,$(filter $(MAKECMDGOALS), auto clean-auto))
            ifeq (${postfix}, default)
                $(error A postfix value of "default" is not allowed)
            endif
            SUBTARGET_PRE_DEBUG = $@_${postfix}
            OPTIONSPOSTFIX_PRE_DEBUG = _${postfix}
        endif
    endif
endif
SUBTARGET = $(SUBTARGET_PRE_DEBUG)
OPTIONSPOSTFIX = $(OPTIONSPOSTFIX_PRE_DEBUG)
ifneq ($(origin debug), undefined)
    ifneq (${debug}, )
        ifneq (,$(filter $(MAKECMDGOALS), auto clean-auto))
            ifeq (,$(filter $(debug), ON FULL))
                $(error Only values for debug are "ON" and "FULL")
            endif
            ifeq (${debug},ON)
                SUBTARGET = $(SUBTARGET_PRE_DEBUG)-debug
                OPTIONSPOSTFIX = $(OPTIONSPOSTFIX_PRE_DEBUG)-debug
            endif
            ifeq (${debug},FULL)
                SUBTARGET = $(SUBTARGET_PRE_DEBUG)-fulldebug
                OPTIONSPOSTFIX = $(OPTIONSPOSTFIX_PRE_DEBUG)-fulldebug
            endif
        endif
    endif
endif

# Definitions

ROOT =	lmp
EXE =	$(ROOT)_${SUBTARGET}
SRC =	$(wildcard *.cpp)
INC =	$(wildcard *.h)
OBJ = 	$(SRC:.cpp=.o)

# Package variables

PACKAGE = asphere cohesionmodels molecule pascal poems rigid surfacemodels voronoi

PACKUSER = 

PACKLIB = poems voronoi 

PACKALL = $(PACKAGE) $(PACKUSER)

PACKAGEUC = $(shell echo $(PACKAGE) | tr a-z A-Z)
PACKUSERUC = $(shell echo $(PACKUSER) | tr a-z A-Z)

YESDIR = $(shell echo $(@:yes-%=%) | tr a-z A-Z)
NODIR  = $(shell echo $(@:no-%=%) | tr a-z A-Z)

# List of all targets

help:
	@echo ''
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'
	@echo 'make purge               purge obsolete copies of package sources'
	@echo 'make tar                 create lmp_src.tar.gz of src dir and packages'
	@echo 'make makelib             create Makefile.lib for static library build'
	@echo 'make makeshlib           create Makefile.shlib for shared library build'
	@echo 'make makelist            create Makefile.list used by old makes'
	@echo 'make -f Makefile.lib machine      build LAMMPS as static library for machine'
	@echo 'make -f Makefile.shlib machine    build LAMMPS as shared library for machine'
	@echo 'make -f Makefile.list machine     build LAMMPS from explicit list of files'
	@echo 'make stubs               build dummy MPI library in STUBS'
	@echo 'make install-python      install LAMMPS wrapper in Python'
	@echo ''
	@echo 'make package             list available packages'
	@echo 'make package-status      status of all packages'
	@echo 'make yes-package         install a single package in src dir'
	@echo 'make no-package          remove a single package from src dir'
	@echo 'make yes-all             install all packages in src dir'
	@echo 'make no-all              remove all packages from src dir'
	@echo 'make yes-standard        install all standard packages'
	@echo 'make no-standard         remove all standard packages'
	@echo 'make yes-user            install all user packages'
	@echo 'make no-user             remove all user packages'
	@echo 'make no-lib              remove all packages with external libs'
	@echo ''
	@echo 'make package-update      replace src files with package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo 'make package-diff        diff src files against package files'
	@echo ''
	@echo 'make machine             build LAMMPS where machine is one of:'
	@echo ''
	@files="`ls MAKE/ -I Makefile.user* -I Makefile.user_default -I CMakeLists_Catalyst.txt -I Makefile.auto.options* -I deprecated`"; \
	  for file in $$files; do head -1 MAKE/$$file; done
	@echo ''

# Build the code

.DEFAULT:
	@test -f MAKE/Makefile.$@
	@if [ ! -d Obj_${SUBTARGET} ]; then mkdir Obj_${SUBTARGET}; fi
	@$(SHELL) Make.sh style
	@$(SHELL) Make.sh models
	@cp MAKE/Makefile.$@ Obj_${SUBTARGET}/Makefile
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@cp Makefile.package Makefile.package.settings Obj_${SUBTARGET}
	@cd Obj_${SUBTARGET}; \
	$(MAKE) $(MFLAGS) "OBJ = $(OBJ)" "INC = $(INC)" "SHFLAGS =" \
	  "EXE = ../$(EXE)" ../$(EXE)

# Remove machine-specific object files

clean:
	@echo 'make clean-all           delete all object files'
	@echo 'make clean-machine       delete object files for one machine'

clean-all:
	rm -rf Obj_* *.a *.so MAKE/Makefile.auto.options* $(ROOT)_*

clean-%:
	rm -rf Obj_$(@:clean-%=%)$(OPTIONSPOSTFIX) MAKE/Makefile.$(@:clean-%=%).options$(OPTIONSPOSTFIX) $(ROOT)_$(@:clean-%=%)$(OPTIONSPOSTFIX) lib$(ROOT)_$(@:clean-%=%)$(OPTIONSPOSTFIX).a lib$(ROOT)_$(@:clean-%=%)$(OPTIONSPOSTFIX).so

purge: Purge.list
	@echo 'Purging obsolete and auto-generated source files'
	@for f in `grep -v '#' Purge.list` ;		\
	    do test -f $$f && rm $$f && echo $$f || : ;		\
	done

# Create a tarball of src dir and packages

tar:
	@cd STUBS; make clean
	@cd ..; tar cvzf src/$(ROOT)_src.tar.gz \
	  src/Make* src/Package.sh src/MAKE src/*.cpp src/*.h src/STUBS \
	  $(patsubst %,src/%,$(PACKAGEUC)) $(patsubst %,src/%,$(PACKUSERUC)) \
          --exclude=*/.svn
	@cd STUBS; make
	@echo "Created $(ROOT)_src.tar.gz"

# Make MPI STUBS library

stubs:
	@cd STUBS; make clean; make

# Create Makefile.lib, Makefile.shlib, and Makefile.list

makelib:
	@$(SHELL) Make.sh style

makeshlib:
	@$(SHELL) Make.sh style

makelist:
	@$(SHELL) Make.sh style

# install LAMMPS shared lib and Python wrapper for Python usage

install-python:
	@python ../python/install.py

# Package management

package:
	@echo 'Standard packages:' $(PACKAGE)
	@echo ''
	@echo 'User-contributed packages:' $(PACKUSER)
	@echo ''
	@echo 'make package             list available packages'
	@echo 'make package-status      status of all packages'
	@echo 'make yes-package         install a single package in src dir'
	@echo 'make no-package          remove a single package from src dir'
	@echo 'make yes-all             install all packages in src dir'
	@echo 'make no-all              remove all packages from src dir'
	@echo 'make yes-standard        install all standard packages'
	@echo 'make no-standard         remove all standard packages'
	@echo 'make yes-user            install all user packages'
	@echo 'make no-user             remove all user packages'
	@echo 'make no-lib              remove all packages with external libs'
	@echo ''
	@echo 'make package-update      replace src files with package files'
	@echo 'make package-overwrite   replace package files with src files'
	@echo 'make package-diff        diff src files against package file'

yes-all:
	@for p in $(PACKALL); do $(MAKE) yes-$$p; done

no-all:
	@for p in $(PACKALL); do $(MAKE) no-$$p; done

yes-standard:
	@for p in $(PACKAGE); do $(MAKE) yes-$$p; done

no-standard:
	@for p in $(PACKAGE); do $(MAKE) no-$$p; done

yes-user:
	@for p in $(PACKUSER); do $(MAKE) yes-$$p; done

no-user:
	@for p in $(PACKUSER); do $(MAKE) no-$$p; done

no-lib:
	@for p in $(PACKLIB); do $(MAKE) no-$$p; done

yes-%:
	@if [ ! -e Makefile.package ]; \
	  then cp Makefile.package.empty Makefile.package; fi
	@if [ ! -e Makefile.package.settings ]; \
	  then cp Makefile.package.settings.empty Makefile.package.settings; fi
	@if [ ! -e $(YESDIR) ]; then \
	  echo "Package $(@:yes-%=%) does not exist"; \
	elif [ -e $(YESDIR)/Install.sh ]; then \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	else \
	  echo "Installing package $(@:yes-%=%)"; \
	  cd $(YESDIR); $(SHELL) ../Install.sh 1; cd ..; \
		$(SHELL) Depend.sh $(YESDIR) 1; \
	fi;

no-%:
	@if [ ! -e $(NODIR) ]; then \
	  echo "Package $(@:no-%=%) does not exist"; \
	elif [ -e $(NODIR)/Install.sh ]; then \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
	else \
	  echo "Uninstalling package $(@:no-%=%)"; \
	  cd $(NODIR); $(SHELL) ../Install.sh 0; cd ..; \
		$(SHELL) Depend.sh $(NODIR) 0; \
        fi;

# status = list src files that differ from package files
# update = replace src files with newer package files
# overwrite = overwrite package files with newer src files
# diff = show differences between src and package files

package-status:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p status; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p status; done

package-update:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p update; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p update; done

package-overwrite:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p overwrite; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p overwrite; done

package-diff:
	@for p in $(PACKAGEUC); do $(SHELL) Package.sh $$p diff; done
	@echo ''
	@for p in $(PACKUSERUC); do $(SHELL) Package.sh $$p diff; done
