# -*-Makefile-*-
# Only used for generating HTML listings (not pages)

all:
	chmod 700 page2html ch2html dir2html search

%.html: $(BASILISK)/cgi-bin/Makefile 		\
	$(BASILISK)/templates/page.static 	\
	$(BASILISK)/page2html
	$(BASILISK)/page2html $* > $@
