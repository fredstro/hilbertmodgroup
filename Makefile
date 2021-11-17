
build:
	sage setup.py build_ext --inplace $* 2>&1

install:
	sage setup.py install

test: install
	sage -t src/hilbert_modgroup

clean:
	rm -rf src/hilbert_modgroup/*.c
	rm -rf src/hilbert_modgroup/*.so
	rm -rf src/hilbert_modgroup/*.cpp
	rm -rf build
	rm -rf dist
PHONY: run clean