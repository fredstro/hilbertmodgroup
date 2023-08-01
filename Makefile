# SageMath is needed in path unless we use docker
ifeq (docker,$(findstring docker, $(MAKECMDGOALS)))
  DOCKER = 1
endif
ifeq (,$(DOCKER)$(shell which sage))
  $(error "No sage in $(PATH). Please install sage first or add the relevant path.")
endif
# Command line arguments:
# Set to 0 to use local source instead of GIT
REMOTE_SRC=1
# Set to desired git branch
GIT_BRANCH=main
# Port for notebook
NBPORT=8888
# Arguments for tox
TOX_ARGS=doctest,coverage,pycodestyle,relint,codespell,flake8


ifeq (0,$(REMOTE_SRC))
  TAG=local
  GIT_BRANCH=
else
  TAG=$(GIT_BRANCH)
endif

EXAMPLES_ARGS:=$(GIT_BRANCH) $(NBPORT)

build:
	sage -python setup.py build_ext --inplace $* 2>&1

sdist:
	sage -python setup.py sdist

install:
	sage -python setup.py install

test:
	sage -t src/*

examples:
	sage -n jupyter --no-browser --ip='0.0.0.0' --port=$(NBPORT)\
                            --notebook-dir=examples\
                            --NotebookApp.custom_display_url=http://127.0.0.1:$(NBPORT)\
                            --NotebookApp.use_redirect_file=False\
                            --NotebookApp.browser=x-www-browser

tox:
	sage -pip install tox meson
	sage --python -m tox src -c tox.ini -e $(TOX_ARGS)

docker:
	docker build --build-arg GIT_BRANCH=$(GIT_BRANCH) --build-arg REMOTE_SRC=$(REMOTE_SRC) -t hilbertmodgroup-$(TAG) .

docker-rebuild:
	docker build --build-arg GIT_BRANCH=$(GIT_BRANCH) --build-arg REMOTE_SRC=$(REMOTE_SRC) --no-cache -t hilbertmodgroup-$(TAG) .

docker-test: docker
	docker run -it -e GIT_BRANCH=$(GIT_BRANCH) --init hilbertmodgroup-$(TAG) test

docker-examples: docker
	docker run -p $(NBPORT):$(NBPORT) -it -e GIT_BRANCH=$(GIT_BRANCH) -e NBPORT=$(NBPORT) --init hilbertmodgroup-$(TAG) examples $(EXAMPLES_ARGS)

docker-tox: docker
	docker run -it -e GIT_BRANCH=$(GIT_BRANCH) -e TOX_ARGS=$(TOX_ARGS) --init hilbertmodgroup-$(TAG) tox

docker-shell: docker
	docker run -it -e GIT_BRANCH=$(GIT_BRANCH) --init hilbertmodgroup-$(TAG) shell

docker-sage: docker
	docker run -it -e GIT_BRANCH=$(GIT_BRANCH) --init hilbertmodgroup-$(TAG) run


clean:
	rm -rf src/hilbert_modgroup/*.c
	rm -rf src/hilbert_modgroup/*.so
	rm -rf src/hilbert_modgroup/*.cpp
	rm -rf hilbert_modular_group.egg-info
	rm -rf build
	rm -rf dist
	rm -rf Library
.PHONY: build sdist install test examples clean shell tox docker docker-examples docker-examples docker-rebuild \
		docker-shell docker-test docker-tox