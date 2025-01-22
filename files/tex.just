# remove "\pandocbounded{...}" keep "..."
[no-cd]
update-qmd file:
    #!/usr/bin/env bash
    sed -E -i '' 's/\\pandocbounded\{([^}]*)\}/\1/g' {{ file }}

# change from "\[ ... \]" to "$$ ... $$" and from "\(...\)" to "$...$"
[no-cd]
update-mathenv file:
    #!/usr/bin/env bash
    sed -i '' 's/\\\[/\$\$/g' {{ file }}
    sed -i '' 's/\\\]/\$\$/g' {{ file }}
    sed -i '' 's/\\(/\$/g' {{ file }}
    sed -i '' 's/\\)/\$/g' {{ file }}

[no-cd]
clean:
    -rm *.{aux,bbl,blg,log,tex}
    -rm trackchanges.sty *.cls

[no-cd]
keep-document file:
    #!/usr/bin/env bash
    temp_file=$(mktemp)
    awk '/\\begin\{document\}/ {flag=1; next} /\\end\{document\}/ {flag=0} flag' {{ file }} | awk '!/\\maketitle/' > $temp_file
    mv $temp_file {{ file }}

[private]
update-agu-format file:
    #!/usr/bin/env bash
    sed -i '' '/\\usepackage\[[^]]*\]{natbib}/d; /\\bibliographystyle{plainnat}/d' {{ file }}
    sed -E -i '' 's/\\citep\[([^]]*)\]\[([^]]*)\]\{([^}]*)\}/\\cite<\1>[\2]{\3}/g' {{ file }}
    sed -i '' 's/\\citep/\\cite/g' {{ file }}
    sed -i '' 's/\\citet/\\citeA/g' {{ file }}

_update-agu-format_linux file:
    sed -i '' '/\\usepackage\[\]\{natbib\}/d; /\\bibliographystyle\{plainnat\}/d' {{ file }}