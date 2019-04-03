#!/usr/bin/env python
#
# Copyright (C) 2009-2011 University of Edinburgh
#
# This file is part of IMUSim.
#
# IMUSim is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# IMUSim is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with IMUSim.  If not, see <http://www.gnu.org/licenses/>.
import os

depsOK = True

try:
    import numpy
except ImportError:
    depsOK = False
    print("NumPy should be installed first from suitable binaries.")
    print("See http://numpy.scipy.org/")

try:
    import scipy
except ImportError:
    depsOK = False
    print("SciPy should be installed first from suitable binaries.")
    print("See http://www.scipy.org/")

try:
    import matplotlib
except ImportError:
    depsOK = False
    print("Matplotlib should be installed first from suitable binaries.")
    print("See http://matplotlib.sf.net/")

try:
    import cython

    if not all(map(os.path.exists, ['imusim/maths/quaternions.c', 'imusim/maths/natural_neighbour.c',
                                    'imusim/maths/quat_splines.c', 'imusim/maths/vectors.c'])):
        depsOK = False
        print("You need to generate c modules by running 'cython -a imusim/maths/*.pyx'")
except ImportError:
    depsOK = False
    print("Cython need to be installed, we need it to generate quaternions lib")
try:
    import mayavi
except ImportError:
    try:
        import enthought.mayavi
    except ImportError:
        depsOK = False
        print("Mayavi should be installed first from suitable binaries.")
        print("See http://code.enthought.com/projects/mayavi/")

try:
    from setuptools import setup, find_packages
    from setuptools.extension import Extension
    if depsOK:
        setup(
            name = "imusim",
            version = "0.2",
            author = "Alex Young and Martin Ling",
            license = "GPLv3",
            url = "http://www.imusim.org/",
            install_requires = ["simpy>=2.3,<3", "pyparsing"],
            packages = find_packages(),
            include_dirs = [numpy.get_include()],
            ext_modules = [
                Extension("imusim.maths.quaternions",
                    ['imusim/maths/quaternions.c']),
                Extension("imusim.maths.quat_splines",
                    ['imusim/maths/quat_splines.c']),
                Extension("imusim.maths.vectors",['imusim/maths/vectors.c']),
                Extension("imusim.maths.natural_neighbour",[
                    'imusim/maths/natural_neighbour/utils.c',
                    'imusim/maths/natural_neighbour/delaunay.c',
                    'imusim/maths/natural_neighbour/natural.c',
                    'imusim/maths/natural_neighbour.c'])]
        )
except ImportError:
    print("Setuptools must be installed - see http://pypi.python.org/pypi/setuptools")
