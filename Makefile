logfile := SBM_$(shell date +%F)
all:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex --metadata date="`date +%D`" -o SBM.pdf
	mv SBM.pdf out/
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o SBM.docx 
	mv SBM.docx out/
pdf:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o SBM.pdf
	mv SBM.pdf out/
tex:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o SBM.tex
xelatex: 
	xelatex SBM.tex
	xelatex SBM.tex
	mv SBM.pdf out/SBM_tex.pdf
	rm SBM.*
doc:
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o SBM.docx 
	mv SBM.docx out/
log:
	cp out/SBM.pdf out/archive/$(logfile).pdf	
	cp out/SBM.docx out/archive/$(logfile).docx
