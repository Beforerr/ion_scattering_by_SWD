env-install:
	mamba env create --file environment.yml

env-update:
	mamba env update --file environment.yml

publish:
	nbqa isort notebooks/*.ipynb
	nbqa black notebooks/*.ipynb
	quarto publish gh-pages --no-prompt