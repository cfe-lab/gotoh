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

# Launching the webserver Docker container

To build the Docker container, use the included Makefile in the root directory.
This is currently tested with `make build-base`. This builds the Docker image.
You can launch the Docker image with `make run-base-local`.
(Root permissions are needed to run these commands.)

Once the Docker image is running, you can launch the FastAPI using
```
cd /webserver/
./runserver.sh
```

There are two API endpoints `/align_it` (nucleotide) and `align_it_aa` (amino acid).
See the FastAPI docs `<local-address>/docs` for details.

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
    git+https://github.com/cfe-lab/gotoh.git@v0.3.0#egg=gotoh&subdirectory=alignment/gotoh
```

# Ruby Bindings

The `/ruby` directory contains the directory structure required to build a
gem containing Gotoh's Ruby bindings.  This gem, named `gotoh`, contains
a module called `Gotoh` that holds the `align_it` and `align_it_aa` functions.

To build Gotoh's Ruby bindings (which are called "alignment" instead of "gotoh"
for arbitrary historical reasons), you'll need to have the following installed:

- Ruby
- The Ruby header files (installed by default with
  [DevKit](https://rubyinstaller.org/add-ons/devkit.html) on Windows or
  [RVM](https://rvm.io/) on Linux)
- A C++ compiler

The "canonical" environment for testing and building this package is the 
dev container defined in the `.devcontainer` directory.  This system is based 
on the `cfe_ubuntu` Ruby image based on Ubuntu 20.04 (this old version is 
required to run Ruby 2.2.2), and stripped down to only the parts necessary 
for Ruby.  This environment is also used by the CI/CD pipeline for testing
and building.

In this environment (or another suitable environment), you should be able to 
build the extension module using the `build_gem.bash` script.  By default this 
will make a gem with the version number `0.1.0.pre`; set the environment variable 
`CFE_GOTOH_VERSION` before building to assign a proper version number.

## Manually building the Ruby bindings

If you want to manually build the bindings, e.g. for testing/debugging/development,
copy `/ruby/ext/gotoh/extconf.rb` to the `/alignment/gotoh/` directory and 
change the `create_makefile('gotoh/gotoh')` line to `create_makefile('gotoh')`;
then run
```
    ruby extconf.rb  # generates a Makefile and other supporting files
    make
```
which will build `gotoh.so`, which can be `require`d from Ruby.
