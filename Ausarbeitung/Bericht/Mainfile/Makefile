DIRNAME               =  $(notdir $(CURDIR))
MAINFILE              = hdllab2016_01

LATEX                 = pdflatex
LATEXPARAMS           = --synctex=1 --shell-escape --interaction=nonstopmode
PDFVIEWER             = okular

SUFFIX_NORMAL         = .pdf

WORDCOUNT_PAGE_OFFSET = 2



FILECOMMAND	= $(MAINFILE)


build:
	$(LATEX) $(LATEXPARAMS) $(FILECOMMAND)
	$(LATEX) $(LATEXPARAMS) $(FILECOMMAND)
	$(LATEX) $(LATEXPARAMS) $(FILECOMMAND)

all: build view
.PHONY: all

view:
	$(PDFVIEWER) $(MAINFILE)$(SUFFIX_NORMAL) &> /dev/null &

.PHONY: cleanall
cleanall:
	rm -f *.dvi *.log *.bak *.aux *.bbl *.blg *.idx *.ps *.eps *.toc *.out *~ *.synctex.gz *.ind *.ilg *.lof *.lot *-blx.bib *.run.xml *.xmpi *.acn *.acr *.alg *.glg *.glo *.gls *.slg *.slo *.sls *.ist $(MAINFILE)$(SUFFIX_NORMAL)

.PHONY: clean
clean:
	rm -f *.log *.bak *.aux *.bbl *.blg *.idx *.toc *.out *~ *.synctex.gz  *.ind *.ilg *.lof *.lot *-blx.bib *.run.xml *.xmpi *.acn *.acr *.alg *.glg *.glo *.gls *.slg *.slo *.sls *.ist
	
count: $(MAINFILE)$(SUFFIX_NORMAL)
	@$(eval WORD_COUNT=$(shell pdftops -f 2 $(MAINFILE)$(SUFFIX_NORMAL) - | ps2ascii | wc -w))
	@echo "Word count (approx.): $(WORD_COUNT)"
