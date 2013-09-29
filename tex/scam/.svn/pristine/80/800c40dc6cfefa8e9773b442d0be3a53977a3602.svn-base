all: scam_nips.pdf paper.pdf supplement.pdf scam.tgz 

scam_nips.pdf: scam_nips.tex local.bib additive_faithful.tex experiments.tex appendix.tex introduction.tex related_work.tex sparsistency.tex discussion.tex optimization.tex summary.tex
	pdflatex scam_nips
	bibtex scam_nips
#	/bin/cat scam_nips.bbl | /usr/bin/sed -e 's/\\begin\(.*\)/\\begin\1\\small\\\parskip=2pt\\itemsep=3pt/' > t.bbl; /bin/mv t.bbl scam_nips.bbl
	pdflatex scam_nips
	pdflatex scam_nips

paper.pdf: scam_nips.pdf paper.tex
	pdflatex paper

supplement.pdf: scam_nips.pdf supplement.tex
	pdflatex supplement

scam.tgz: paper.pdf supplement.pdf matlab-code-and-data.tar.gz
	tar -cvf - paper.pdf supplement.pdf matlab-code-and-data.tar.gz | gzip > scam.tgz