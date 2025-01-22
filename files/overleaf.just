mod tex 'tex.just'

overleaf_repo := x'https://git@git.overleaf.com/${OVERLEAF_ID:-}'
tex_file := "overleaf/index.tex"

[private]
clone-overleaf:
    [ -d overleaf ] || git clone {{overleaf_repo}} overleaf

update-overleaf: sync-overleaf
    cd overleaf; git add .; git commit -am "update"; git push

sync-overleaf: tex-render post-render
    -rsync _manuscript/_tex/ overleaf/ -r
    -rsync article_files overleaf -a --size-only

tex-render:
    #!/usr/bin/env bash
    touch files/bibexport.bib
    quarto render --profile man --to $OVERLEAF_FMT

post-render file=tex_file:
    just tex::keep-document {{file}}
    just tex::update-mathenv {{file}}
    just tex::update-qmd {{file}}
    just tex::clean
