main := SBM
logfile := $(main)_$(shell date +%F)
all:
	pandoc main.md --template=pandoc/default.latex --verbose --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o $(main).pdf
	mv $(main).pdf out/
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o $(main).docx
	mv $(main).docx out/
	pandoc main.md --template=pandoc/default.latex --verbose --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o $(main).tex
	mv $(main).tex out/
	cp out/$(main).pdf out/archive/$(logfile).pdf
	cp out/$(main).docx out/archive/$(logfile).docx
	cp out/$(main).tex out/archive/$(logfile).tex
pdf:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o $(main).pdf
	mv $(main).pdf out/
doc:
	pandoc main.md --reference-doc=pandoc/custom-reference.docx --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc -o $(main).docx
	mv $(main).docx out/
tex:
	pandoc main.md --template=pandoc/default.latex --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref --citeproc --pdf-engine=xelatex -o $(main).tex
	mv $(main).tex out/
log:
	cp out/$(main).pdf out/archive/$(logfile).pdf
	cp out/$(main).docx out/archive/$(logfile).docx
	cp out/$(main).tex out/archive/$(logfile).tex
cleanall:
	ls outfile* | grep -v main.md | xargs rm