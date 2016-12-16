
OBJECTS_DIR= obj
REALMAKEFILE=../Makefile.in

all: solver

compile_source:
	@$(MAKE) -C ../ --no-print-directory solver

solver: FORCE compile_source
	@(cd $(OBJECTS_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

clean: FORCE
	@rm -rf obj bin *~
	
FORCE:
	@mkdir -p obj bin
