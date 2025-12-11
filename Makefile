lstm.md: lstm.qmd
	poetry run quarto render lstm.qmd --to gfm --output lstm.md

data:
	cd data && \
	poetry run quarto render covid_variants.qmd --to gfm --output covid_variants.md


.PHONY: data