import 'files/quarto.just'
import 'files/overleaf.just'

default:
    just --list

ensure-env: install-julia-deps clone-overleaf
    pixi install
    
install-julia-deps:
    #!/usr/bin/env -S julia --project
    using Pkg;
    Beforerr = PackageSpec(url="https://github.com/Beforerr/beforerr.jl");
    Pkg.develop([Beforerr]);
    Pkg.instantiate();

publish:
    nbqa isort notebooks/*.ipynb
    nbqa black notebooks/*.ipynb
    quarto publish gh-pages --no-prompt