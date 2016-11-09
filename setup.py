from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'phASER Read Variant Mapper',
  ext_modules = cythonize("read_variant_map.py"),
)
