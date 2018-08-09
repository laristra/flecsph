#Tex Makefile
FILE=notes

all: $(FILE).pdf

$(FILE).pdf : notes.tex \
  FleCSI \
  Tree \
  General_algorithm \
  IO \
  Formulations \
  Kernels \
  Equations_of_state \
  Initial_data \
  Gravitation \
  Time_integration \
  Resolution_order \
  Applications
	

view: $(FILE).pdf display

.PHONY: clean

clean:
	rm -rf *~ *.idx *.lot *.toc *.lof *.aux *.blg *.out *.bbl *.log $(FILE).pdf

display:
	open $(FILE).pdf

$(FILE).pdf: $(FILE).tex
	pdflatex $(FILE)
	#bibtex $(FILE)   # Uncomment for once we have bib file
	pdflatex $(FILE)
	#pdflatex $(FILE) # Uncomment for once we have bib file
