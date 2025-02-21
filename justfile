import 'files/quarto.just'
import 'files/overleaf.just'

default:
    just --list

ensure-env: install-julia-deps clone-overleaf
    pixi install
    
install-julia-deps:
    #!/usr/bin/env julia --threads=auto --project=.
    using Pkg
    Pkg.develop([(;name="Speasy"), (;name="SpaceTools")])
    Pkg.resolve()
    Pkg.instantiate()

publish: quarto-publish

render:
    # quarto render presentations/index.qmd --to pptx
    # cp _site/presentations/index.pptx presentations/index.pptx

    quarto render presentations/SPARTHRB.qmd --to pptx
    cp _site/presentations/SPARTHRB.pptx presentations/_SPARTHRB.pptx

exec-scripts:
    julia --project scripts/scan.jl
    julia --project scripts/tm.jl
    # include("scripts/scan.jl")