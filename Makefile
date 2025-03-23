# Command line arguments:
# Set to 1 to use GIT instead of local source for docker install
REMOTE_SRC=0
# Set to desired git branch if using GIT source
GIT_BRANCH=main
# Port for notebook
NBPORT=8888
# Arguments for tox
TOX_ARGS=doctest,coverage,pycodestyle,relint,codespell,flake8,passagemath

ifeq (0,$(REMOTE_SRC))
  TAG=local
  GIT_BRANCH=
else
  TAG=$(GIT_BRANCH)
endif

EXAMPLES_ARGS:=$(GIT_BRANCH) $(NBPORT)
SAGE_EXEC:=$(shell which sage)

sage-build: sage-check
	sage -pip install build
	sage -python -m build .

sage-sdist: sage-check
	sage -python -m build --sdist .

sage-install: sage-check
	sage -pip install .

build: venv-check
	pip install build
	python -m build

sdist: venv-check
	python -m build --sdist

install: venv-check
	pip install .

test: sage-check
	sage -t src/*

sage-examples: sage-check
	sage --notebook=jupyterlab --no-browser --ip='0.0.0.0' --port=$(NBPORT)\
                            --notebook-dir=examples\
                            --ServerApp.custom_display_url=http://127.0.0.1:$(NBPORT)\
                            --ServerApp.use_redirect_file=False\
                            --ServerApp.browser=x-www-browser

examples: venv-check
	pip install jupyterlab
	jupyter lab --no-browser --ip='0.0.0.0' --port=$(NBPORT)\
                            --notebook-dir=examples\
                            --ServerApp.custom_display_url=http://127.0.0.1:$(NBPORT)\
                            --ServerApp.use_redirect_file=False\
                            --ServerApp.browser=x-www-browser

tox: sage-check
	sage -pip install tox meson
	sage --python -m tox -c tox.ini -e $(TOX_ARGS)

docker-build:
	docker build --build-arg GIT_BRANCH=$(GIT_BRANCH) --build-arg REMOTE_SRC=$(REMOTE_SRC) -t hilbertmodgroup-$(TAG) .

docker-rebuild:
	docker build --build-arg GIT_BRANCH=$(GIT_BRANCH) --build-arg REMOTE_SRC=$(REMOTE_SRC) --no-cache -t hilbertmodgroup-$(TAG) .

docker-test: docker-build
	docker run --platform linux/amd64 -it -e GIT_BRANCH=$(GIT_BRANCH) --init hilbertmodgroup-$(TAG) test

docker-examples: docker-build
	docker run --platform linux/amd64 -p $(NBPORT):$(NBPORT) -it -e GIT_BRANCH=$(GIT_BRANCH) -e NBPORT=$(NBPORT) --init hilbertmodgroup-$(TAG) examples $(EXAMPLES_ARGS)

docker-tox: docker-build
	docker run --platform linux/amd64 -it -e GIT_BRANCH=$(GIT_BRANCH) -e TOX_ARGS=$(TOX_ARGS) --init hilbertmodgroup-$(TAG) tox

docker-shell: docker-build
	docker run --platform linux/amd64 -it -e GIT_BRANCH=$(GIT_BRANCH) --init hilbertmodgroup-$(TAG) shell

docker-sage: docker-build
	docker run --platform linux/amd64 -it -e GIT_BRANCH=$(GIT_BRANCH) --init hilbertmodgroup-$(TAG) run

venv-check:
	@if [ "${VIRTUAL_ENV}" = "" ]; then\
		echo "It seems you are trying to build this outside a virtual environment. If you know what you are doing you can comment out this check, otherwise read the README.md file" >&2;\
		exit 1;\
	fi

sage-check:
	@if [ "${SAGE_EXEC}" = "" ]; then\
		echo "No sage in PATH:$(PATH). Please install sage first, add the relevant path or use python build.";\
		exit 1;\
	fi

clean:
	rm -rf src/hilbert_modgroup/*.c
	rm -rf src/hilbert_modgroup/*.so
	rm -rf src/hilbert_modgroup/*.cpp
	rm -rf hilbert_modular_group.egg-info
	rm -rf build
	rm -rf dist
	rm -rf Library
	rm -rf __pycache__
	rm -rf src/*/__pycache__

.PHONY: sage-check venv-check clean shell\
		build sdist install test examples tox\
		sage-build sage-sdist sage-install sage-examples \
 		docker docker-examples docker-rebuild docker-shell docker-test docker-tox\
