overleaf_repo := x'https://git@git.overleaf.com/${OVERLEAF_ID:-}'
tex_file := "overleaf/index.tex"

[private]
clone-overleaf:
  [ -d overleaf ] || git clone {{overleaf_repo}} overleaf

update-overleaf: sync-overleaf
  cd overleaf; git add .; git commit -am "update"; git push

sync-overleaf: tex-render tex-keep-document tex-update-math-env tex-update-qmd tex-export-bib tex-clean

tex-render:
  #!/usr/bin/env bash
  touch files/bibexport.bib
  quarto render --profile man --to $OVERLEAF_FMT
  rsync _manuscript/_tex/ overleaf/ -r
  rsync article_files overleaf -a --size-only

[private]
tex-export-bib:
  #!/usr/bin/env bash
  command -v bibexport >/dev/null || tlmgr install bibexport
  bibexport -o files/bibexport.bib --nosave article.aux

[private]
tex-keep-document file = tex_file:
  #!/usr/bin/env bash
  temp_file=$(mktemp)
  awk '/\\begin\{document\}/ {flag=1; next} /\\end\{document\}/ {flag=0} flag' {{file}} | awk '!/\\maketitle/' > $temp_file
  mv $temp_file {{file}}

tex-update-math-env file = tex_file:
    #!/usr/bin/env bash
    sed -i '' 's/\\\[/\$\$/g' {{ file }}
    sed -i '' 's/\\\]/\$\$/g' {{ file }}
    sed -i '' 's/\\(/\$/g' {{ file }}
    sed -i '' 's/\\)/\$/g' {{ file }}
  
tex-update-qmd file = tex_file:
    #!/usr/bin/env bash
    sed -E -i '' 's/\\pandocbounded\{([^}]*)\}/\1/g' {{ file }}

[private]
tex-clean:
  -rm *.{aux,bbl,blg,log,tex}
  -rm trackchanges.sty *.cls