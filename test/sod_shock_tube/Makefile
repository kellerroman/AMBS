SUBDIRS := ${shell find ./*  -maxdepth 0 -type d -print }
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))
CALC_ERROR=calc_error
.PHONY: subdirs $(SUBDIRS) clean $(SUBCLEAN)
 
subdirs: $(SUBDIRS)

clean:  $(SUBCLEAN) 
	@rm -vf *.o *.mod *.ufo *.dat *.plt *.png *.vtk *.txt *.bin gridgen *.silo *.dx stubex $(CALC_ERROR)

list_targets:
	@echo $(SUBDIRS)

$(CALC_ERROR):$(CALC_ERROR).F90
	@gfortran $< -o $@
     
$(SUBDIRS):
	@$(MAKE) -C $@ #--no-print-directory

$(SUBCLEAN): %.clean:
	@$(MAKE) -C $* clean --no-print-directory
	
