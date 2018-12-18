LATEXMKFLAGS = -quiet -pdf -pdflatex='pdflatex -synctex=1'
CLEANEXTS = fdb_latexmk,aux,bcf,fls,lof,log,lot,nlo,out,run.xml,synctex.gz,toc,bbl,blg,nlg,nls,vrb,nav,snm
.PHONY: all

all: dissertation-renfro-embedded.pdf presentation-renfro-slides.pdf presentation-renfro-notes.pdf

dissertation-renfro.pdf:
	latexmk $(LATEXMKFLAGS) dissertation-renfro

presentation-renfro-slides.pdf:
	latexmk $(LATEXMKFLAGS) presentation-renfro-slides

presentation-renfro-notes.pdf:
	latexmk $(LATEXMKFLAGS) presentation-renfro-notes

dissertation-renfro-embedded.pdf: dissertation-renfro.pdf
	gs -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress \
		-dEmbedAllFonts=true -sOutputFile=dissertation-renfro-embedded.pdf \
		-f dissertation-renfro.pdf

clean:
	rm -f dissertation-renfro.{$(CLEANEXTS)}
	rm -f presentation-renfro-slides.{$(CLEANEXTS)}
	rm -f presentation-renfro-notes.{$(CLEANEXTS)}

realclean: clean
	rm -f dissertation-renfro.pdf
	rm -f dissertation-renfro-embedded.pdf
	rm -f presentation-renfro-notes.pdf presentation-renfro-slides.pdf
