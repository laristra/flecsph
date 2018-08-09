#Tex Makefile
FILE=notes

all: $(FILE).pdf

view: $(FILE).pdf display

.PHONY: clean

clean:
	rm -rf *~ *.idx *.lot *.toc *.lof *.aux *.blg *.out *.bbl *.log $(FILE).pdf

display:
	open $(FILE).pdf

$(FILE).pdf: $(FILE).tex \
              FleCSI.tex \
                Tree.tex \
   General_algorithm.tex \
                  IO.tex \
        Formulations.tex \
             Kernels.tex \
  Equations_of_state.tex \
        Initial_data.tex \
         Gravitation.tex \
    Time_integration.tex \
    Resolution_order.tex \
        Applications.tex
	pdflatex $(FILE)
	#bibtex $(FILE)   # Uncomment for once we have bib file
	pdflatex $(FILE)
	#pdflatex $(FILE) # Uncomment for once we have bib file
