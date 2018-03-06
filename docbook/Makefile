# Makefile to generate the documentation of raytrace
# $Id: Makefile,v 1.2 2008-11-21 20:02:32 thernis Exp $

# programs
XMLTO = xmlto
LATEX = latex

# main file
DOC = raytraceman.docbook_06
# XSL style sheet modification file
XSL = config.xsl


all: html

# Rules to generate html doc
html: $(DOC).xml $(XSL)
	$(XMLTO) html -m $(XSL) $(DOC).xml
	$(LATEX) tex-math-equations.tex
# 	dvipng -O 0cm,-0.1cm -Q 10 tex-math-equations.dvi 
# 	dvipng tex-math-equations.dvi
	dvi2bitmap --magnification=5 --scale=5 tex-math-equations.dvi

clean:
	rm -f tex-math-equations.*
	
xslt: $(DOC).xml $(XSL)
	xsltproc -o test.html /usr/share/xml/docbook/stylesheet/nwalsh/xhtml/docbook.xsl $(DOC).xml

