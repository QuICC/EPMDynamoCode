#
# This is only a very simple make file for documentation generation
# The real makefile has to be generated through a call to cmake
#

#SOURCE CODE DIRECTORY LAYOUT
#Set the layout of the source code directory
DOCDIR		= ./Doc/
CODEDOCDIR	= ./CodeDoc/

#DOCUMENTATION CREATION RULES
#Set make rules for documentation generation
#--------------------------------------------------------------------
.PHONY: CodeDoc
CodeDoc:	
	doxygen

.PHONY: doc
doc: $(DOCDIR)
#--------------------------------------------------------------------

#CLEANING RULES
#Set make rules for cleaning source directory
#--------------------------------------------------------------------
.PHONY: cleanCodeDoc
cleanCodeDoc:
	rm -r -f $(CODEDOCDIR)latex
	rm -r -f $(CODEDOCDIR)html
#--------------------------------------------------------------------
