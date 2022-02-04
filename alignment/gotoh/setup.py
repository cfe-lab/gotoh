from setuptools import setup, Extension

about = {}
with open("about.py") as f:
    exec(f.read(), about)

gotoh = Extension(about['__title__'],
                  sources = ['gotoh.cpp'],
                  define_macros=[('__PYTHON__', None)])

setup( name = about['__title__'],
       version = about['__version__'],
       description=about['__description__'],
       url=about['__url__'],
       ext_modules = [gotoh],
       zip_safe = False) #avoid headache with permissions on ~/.python-eggs
