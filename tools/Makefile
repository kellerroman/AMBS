
OBJECTS_DIR= obj
REALMAKEFILE=../Makefile.in

all: solver

solver: FORCE
	@(cd $(OBJECTS_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

clean: FORCE
	@rm -rf obj bin *~
	
FORCE:
	@mkdir -p obj bin
