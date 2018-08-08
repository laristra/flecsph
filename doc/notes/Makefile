#Tex Makefile
FILE=notes

all: $(FILE).pdf

view: $(FILE).pdf display

.PHONY: clean

clean:
	rm -rf *~ *.idx *.lot *.toc *.lof *.aux *.blg *.out *.bbl *.log $(FILE).pdf

display:
	open $(FILE).pdf

$(FILE).pdf: $(FILE).tex
	pdflatex $(FILE)
	#bibtex $(FILE) # Uncommented once we have bib file
	pdflatex $(FILE)
