from setuptools import setup, Extension

gotoh = Extension(
    "gotoh", sources=["alignment.cpp"], define_macros=[("__PYTHON__", None)],
)

setup(
    name="gotoh",
    version="0.2",
    description="Wrapper for Conan's alignment.cpp code",
    ext_modules=[gotoh],
    zip_safe=False,
    extras_require={'dev': ['pytest']})
