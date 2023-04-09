README.md: README.Rmd
	Rscript -e "rmarkdown::render('README.Rmd', output_file = 'README.md')"
	rm README.html

check:
	Rscript -e "devtools::check()"


site: README.md
	Rscript -e "pkgdown::build_site()"


test:
	Rscript -e "devtools::test()"


doc:
	Rscript -e "devtools::document()"


build:
	Rscript -e "devtools::build()"