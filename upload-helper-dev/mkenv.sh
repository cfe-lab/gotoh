#!/bin/bash

# this script sets up a virtualenv called 'env' in the local directory.
# This virtual env can be activated OUTSIDE of the docker container
# so that files can be edited on the host using 'almost the same'
# python environment as is used inside the container.
# this means, for example, that emacs will use the correct versions
# of python for syntax and type checking when editing files on the host.
# 'almost', because the container could well have a different
# python version as the host os.


# virtualenv -p /usr/local/bin/python3 env
virtualenv -p /usr/bin/python3 env

# set MYPYPATH env variable upon activation as well...
# echo "setenv MYPYPATH :/home/wscott/oracle-stuff/scoracle-dev:" >> env/bin/activate.csh
# echo "setenv PYTHONPATH :/home/wscott/oracle-stuff/scoracle-dev:" >> env/bin/activate.csh

source env/bin/activate
pip3 install --upgrade pip
pip3 install -r ulhelper/requirements.txt
