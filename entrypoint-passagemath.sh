#!/bin/bash

set -e

cd /home/sage/hilbertmodgroup

case $1 in
    test)
        echo "Running passagemath doctests..."
        python -m sage.doctest -p --force-lib --environment=_doctest_environment src/
        ;;
    lint)
        echo "Running linting checks..."
        tox -e ruff,pycodestyle,relint,codespell
        ;;
    tox)
        shift
        echo "Running tox with args: $@"
        tox "$@"
        ;;
    shell)
        /bin/bash
        ;;
    *)
        echo "Unknown command: $1"
        echo "Usage: test | lint | tox [args] | shell"
        exit 1
        ;;
esac