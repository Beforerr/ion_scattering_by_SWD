import 'files/quarto.just'
import 'files/overleaf.just'

default:
    just --list

ensure-env: clone-overleaf
    pixi install

publish:
    nbqa isort notebooks/*.ipynb
    nbqa black notebooks/*.ipynb
    quarto publish gh-pages --no-prompt