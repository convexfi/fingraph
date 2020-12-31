clean:
	rm -v src/*.so src/*.o
	rm -v R/RcppExports.R
	rm -v src/RcppExports.cpp

build:
	Rscript .roxygenize.R

install:
	R CMD INSTALL ../fingraph

test:
	Rscript -e "devtools::test()"

all:
	make build && make install && make test
