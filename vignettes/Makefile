LATEXFILES = *.aux\
	*.bbl\
	*.blg\
	*.ilg\
	*.log\
	*.nlo\
	*.nls\
	*.toc\
	*.aux\
	*.dvi\
	*.map\
	*.out\
	*.figlist\
	*.dep\
	*.dpth\
	*.makefile\
	MetCirc.tex\

	
all: 
	make vignette

	
vignette: MetNet.Rnw
	pdflatex MetNet.tex
	bibtex MetNet
	pdflatex MetNet.tex
	pdflatex MetNet.tex

	
clean:	
	rm -f $(LATEXFILES)
	rm -f .Rhistory
	rm -rf figure/
	rm -rf cache/
