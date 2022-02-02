# Gotoh

Gotoh is a simple sequence aligner used in several of the BC-CfE's bioinformatics tools.

The functionality is launched as a Docker container, and provided via a FastAPI interface.

By building extension modules, Gotoh can also be called from Python directly.

For details on the alignment algorithm used, see the following:

```
Gotoh, Osamu. "An improved algorithm for matching biological sequences."
    Journal of Molecular Biology, Volume 162, Issue 3 (1982): 705-708.

Gotoh, Osamu. "Optimal alignment between groups of sequences and its application
     to multiple sequence alignment." Computer applications in the biosciences: 
     CABIOS 9.3 (1993): 361-370.
```


# Python Bindings

To build Gotoh's Python bindings, you'll need to have the following
installed:

- A supported Python interpreter (e.g. Python 3.9)
- The Python header files (`apt install python-dev` or `yum install
  platform-python-devel` on Linux, usualy installed by default on Windows)
- A C++ compiler

With these installed, you should be able to build the extension module with
```
    python setupy.py build
```

Alternately, you may install the `gotoh` module using `pip` by adding this repo
to your `requirements.txt` file or similar method.
```
    git+https://github.com/cfe-lab/gotoh.git@v0.3a0#egg=gotoh&subdirectory=alignment/gotoh
```

# Ruby Bindings

To build Gotoh's Ruby bindings (which are called "alignment" instead of "gotoh"
for arbitrary historical reasons), you'll need to have the following installed:

- Ruby
- The Ruby header files (installed by default with
  [DevKit](https://rubyinstaller.org/add-ons/devkit.html) on Windows or
  [RVM](https://rvm.io/) on Linux)
- A C++ compiler

With these installed, you should be able to build the extension module with
```
    ruby extconf.rb
    make
```
