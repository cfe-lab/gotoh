#! /usr/bin/env bash

cp ../alignment/gotoh/gotoh.cpp ext/gotoh
gem build gotoh.gemspec
rm ext/gotoh/gotoh.cpp
