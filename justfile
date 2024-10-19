import 'files/quarto.just'
import 'files/overleaf.just'

default:
    just --list

ensure-env: install-julia-deps clone-overleaf
    pixi install
    
install-julia-deps:
    julia --project -e 'using Pkg; Pkg.update();'

publish: quarto-publish

exec-scripts:
    julia --project scripts/scan.jl
    julia --project scripts/tm.jl
    # include("scripts/scan.jl")