ifeq (, $(shell which sage))
 $(error "No sage in $(PATH), please install sage first or add the relevant path.")
endif
ifeq (docker-tox,$(firstword $(MAKECMDGOALS)))
  TOX_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
  $(eval $(TOX_ARGS):;@:)
endif

build:
	sage -python setup.py build_ext --inplace $* 2>&1

sdist:
	sage -python setup.py sdist

install:
	sage -python setup.py install

test:
	sage -t src/*

examples:
	./entrypoint.sh examples

tox:
	sage -tox src -e $(TOX_ARGS)

docker:
	docker build -t hilbertmodgroup .

docker-rebuild:
	docker build --no-cache -t hilbertmodgroup .

docker-test: docker
	docker run -it --init hilbertmodgroup test

docker-examples: docker
	docker run -p 8888:8888 -it --init hilbertmodgroup examples main

docker-tox: docker
	docker run -it --init hilbertmodgroup tox $(TOX_ARGS)

docker-shell: docker
	docker run -it --init hilbertmodgroup shell

docker-sage: docker
	docker run -it --init hilbertmodgroup run


clean:
	rm -rf src/hilbert_modgroup/*.c
	rm -rf src/hilbert_modgroup/*.so
	rm -rf src/hilbert_modgroup/*.cpp
	rm -rf hilbert_modular_group.egg-info
	rm -rf build
	rm -rf dist
	rm -rf Library
.PHONY: build sdist install docker test examples clean shell tox docker docker-examples docker-examples docker-rebuild \
		docker-shell docker-test docker-tox