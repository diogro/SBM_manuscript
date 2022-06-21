all:
	pandoc main.md  --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref  --citeproc  --pdf-engine=xelatex -o main.pdf
	mv main.pdf out/
doc:
	pandoc main.md  --lua-filter=pandoc/scholarly-metadata.lua --lua-filter=pandoc/author-info-blocks.lua  --filter pandoc-crossref  --citeproc -o main.docx 
	mv main.docx out/
clean:
	ls main* | grep -v main.md | xargs rm