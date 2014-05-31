from distutils.core import setup, Extension
import numpy.distutils.misc_util

setup(
        name = "Bloch Simulator C Extension",
        ext_modules=[Extension("bloch_simulator", ["bloch_simulator.c"])],
        include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
        )

