name=readme
${name}.dvi: 
	latex ${name}.tex
	latex ${name}.tex
dvi: 
	latex ${name}.tex
	latex ${name}.tex
${name}.ps: dvi ${name}.tex
	dvips ${name}.dvi -o ${name}.ps
ps: ${name}.ps
${name}.pdf: ${name}.ps
	ps2pdf ${name}.ps
pdf: ${name}.pdf
	gv ${name}.pdf
bib: 
	latex ${name}.tex
	bibtex ${name}
	latex ${name}.tex
	latex ${name}.tex
clean: 
	rm -f ${name}.ps ${name}.pdf ${name}.dvi ${name}.aux
