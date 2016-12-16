
OBJECTS_DIR= ../obj
BIN_DIR= ../bin
REALMAKEFILE=../Makefile.in

all: solver

solver: FORCE
	@(cd $(OBJECTS_DIR) && $(MAKE) -f $(REALMAKEFILE) --no-print-directory)

clean: FORCE
	@rm -rf $(OBJECTS_DIR) $(BIN_DIR) ../*~ *~
	
FORCE:
	@mkdir -p $(OBJECTS_DIR) $(BIN_DIR)  
