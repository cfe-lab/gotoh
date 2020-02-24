# Gotoh

Gotoh is a simple sequence aligner used in several of the Centre's bioinformatics tool's.

By building extension modules, Gotoh can be called from Python or from Ruby.


# Python Bindings

To build Gotoh's Python bindings, you'll need to have the following
installed:

- A supported Python interpreter (e.g. Python 3.6)
- The Python header files (`apt install python-dev` or `yum install
  platform-python-devel` on Linux, usualy installed by default on Windows)
- A C++ compiler

With these installed, you should be able to build the extension module with

    python setupy.py build


# Ruby Bindings

To build Gotoh's Ruby bindings (which are called "alignment" instead of "gotoh"
for arbitrary historical reasons), you'll need to have the following installed:

- Ruby
- The Ruby header files (installed by default with
  [DevKit](https://rubyinstaller.org/add-ons/devkit.html) on Windows or
  [RVM](https://rvm.io/) on Linux)
- A C++ compiler

With these installed, you should be able to build the extenion module with

    ruby extconf.rb
    make
