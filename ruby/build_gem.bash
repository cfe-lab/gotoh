#! /usr/bin/env bash

cp ../alignment/gotoh/gotoh.cpp ext/cfe_gotoh/cfe_gotoh.cpp
gem build cfe_gotoh.gemspec
rm ext/cfe_gotoh/cfe_gotoh.cpp
