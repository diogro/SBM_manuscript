logfile := expressionVariance_$(shell date +%F)
all:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o SBM.pdf
	mv SBM.pdf out/
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o SBM.docx
	mv SBM.docx out/
	cp out/SBM.pdf out/archive/$(logfile).pdf
	cp out/SBM.docx out/archive/$(logfile).docx
pdf:
	pandoc main.md --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o SBM.pdf
	mv SBM.pdf out/
doc:
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o SBM.docx
	mv SBM.docx out/
log:
	cp out/SBM.pdf out/archive/$(logfile).pdf
	cp out/SBM.docx out/archive/$(logfile).docx
cleanall:
	ls outfile* | grep -v main.md | xargs rm