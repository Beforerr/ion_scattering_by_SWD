import 'files/quarto.just'
import 'files/overleaf.just'

default:
    just --list

ensure-env: install-julia-deps clone-overleaf
    pixi install
    
install-julia-deps:
    julia --project -e 'using Pkg; Pkg.update();'

publish:
    nbqa isort notebooks/*.ipynb
    nbqa black notebooks/*.ipynb
    quarto publish gh-pages --no-prompt

exec-scripts:
    julia --project scripts/scan.jl
    # include("scripts/scan.jl")