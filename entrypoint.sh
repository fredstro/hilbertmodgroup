#!/bin/sh

branch="${2:=develop}"
git pull -ff origin $branch
echo "Pulling branch: $branch"
case $1 in
    test)
      echo "Docker container running Sage doctests for the hilbertmodgroup package."
      sage -t src
      ;;
    tox-all)
      echo "Docker container running tox with 'doctest', 'coverage', 'pycodestyle-minimal', 'relint', 'codespell'"
      sage -python -m tox src
      ;;
    examples)
      echo "Docker container with Jupyter Notebook interface to run example notebooks."
      echo "NOTE: The Jupyter Notebook server is only accessible using the URL from outside the container."
      sage -n jupyter --no-browser --ip='0.0.0.0' --port=8888\
                            --notebook-dir=/home/sage/hilbertmodgroup/examples\
                            --NotebookApp.custom_display_url=http://127.0.0.1:8888\
                            --NotebookApp.use_redirect_file=False\
                            --NotebookApp.browser=x-www-browser
      ;;
    run)
      sage
      ;;
    shell)
      /bin/bash
      ;;
    *)
      echo "Unknown command $1. Exiting..."
      ;;
esac
