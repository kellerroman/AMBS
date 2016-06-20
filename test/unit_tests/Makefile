#SUBDIRS = airfoil
SUBDIRS := ${shell find ./* -maxdepth 0 -type d -print }
SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))
.PHONY: subdirs $(SUBDIRS) clean $(SUBCLEAN)
subdirs: $(SUBDIRS)

clean: $(SUBCLEAN)

test:
	@echo $(SUBCLEAN)
$(SUBDIRS):
	@$(MAKE) -C $@ --no-print-directory
$(SUBCLEAN): %.clean:
	@$(MAKE) -C $* clean --no-print-directory
