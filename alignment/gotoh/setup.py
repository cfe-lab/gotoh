from setuptools import setup, Extension

gotoh = Extension('gotoh',
                  sources = ['gotoh.cpp'],
                  define_macros=[('__PYTHON__', None)])

setup (name = 'gotoh',
       version = '0.3.0',
       description = "Wrapper for Gotoh alignment code",
       ext_modules = [gotoh],
       zip_safe = False) #avoid headache with permissions on ~/.python-eggs
