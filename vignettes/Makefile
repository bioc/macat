all:	pdf

macat.tex:    macat.Rnw
	"$(R_HOME)/bin/R" CMD Sweave $<

pdf:	macat.tex
	"$(R_HOME)/bin/R" CMD pdflatex $<
	make clean

clean:
	rm -f *.out *.bbl *.aux *.log *.blg *.toc  *~ *.eps
	rm -f *.nfs*
	rm -f *.dvi *.ps macat-* macat.tex
