all: pdf

pdf:
	pdflatex --shell-escape dgp_pres.tex
	pdflatex dgp.tex
	bibtex dgp.aux
	pdflatex dgp.tex
	pdflatex dgp.tex

clean:
	rm -rf *log *out *snm *toc *aux \#*\# *nav *blg *bbl

.PHONY: all clean pdf

