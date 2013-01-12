Installing IMUSim
=================

These instructions describe how to install IMUSim on your system. IMUSim has
been tested on Linux, Mac OS X and Windows - these instructions cover all three
operating systems.

Python prerequisites
--------------------

IMUSim has several prerequisites which must be installed first: `Python
<http://www.python.org/>`_ 2.6 or 2.7, plus the following major libraries:

    - `NumPy <http://numpy.scipy.org>`_
    - `SciPy <http://www.scipy.org>`_
    - `Matplotlib <http://matplotlib.sf.net/>`_
    - `MayaVi <http://code.enthought.com/projects/mayavi/>`_

Whilst all of these libraries can be compiled and installed from source, by far
the easiest and best way to install them is from a pre-packaged distribution.
There are three main options:

    - Python(x,y) - a free scientific Python distribution with a
      `Windows version <http://www.pythonxy.com/>`_ and a
      `Linux version <http://code.google.com/p/pythonxy-linux/>`_.
    - The `Enthought Python Distribution (EPD)
      <http://www.enthought.com/products/epd.php>`_ - a comercial scientific
      Python distribution for Windows, Mac OS X and Linux, which is `free for
      academic use <http://www.enthought.com/products/edudownload.php>`_.
    - For Linux systems, using your distribution's own Python packages.

The standard EPD install includes all of the components required for IMUSim.

The Python(x,y) installer prompts you to select the components you would like
to install. You can install as much as you like, but the minimum components
required for IMUSim are:

    - NumPy
    - SciPy
    - Matplotlib
    - ETS (the Enthought Tool Suite, including Mayavi)
    - SetupTools
    - IPython

If you are using a Linux system, you can install the libraries from your
package manager. The following command will install all IMUSim's Python
dependencies on Debian or Ubuntu systems::

    $ sudo apt-get install python-numpy python-scipy python-matplotlib \
        mayavi2 ipython python-setuptools python-simpy python-pyparsing

A similar combination of packages should be required on other distributions.

C compiler
----------

IMUSim is partly written in C, so the setup script needs a compiler.

On Linux, install gcc from your package manager if it is not already installed,
and the Python headers (in Debian/Ubuntu, these are in the `python-dev`
package).

On Mac OS X, make sure the `Apple Developer Tools
<http://developer.apple.com/technologies/tools/>`_ are installed.

On Windows, `install the free MinGW compiler
<http://www.mingw.org/wiki/Getting_Started>`_. You may be able to use another
compiler, but only MinGW has been tested.

Installing IMUSim
-----------------

Make sure you have downloaded the latest release from the `IMUSim homepage
<http://www.imusim.org/>`_. Unpack the contents, which should consist of a
single `imusim` directory. Then open a command prompt and change to that
directory.

If you are using the MinGW compiler on Windows, you will need to run the
following commands first to build the C components of IMUSim correctly::

    > PATH=C:\MinGW\bin;%PATH%
    > python setup.py build_ext --compiler=mingw32 -liberty

Other systems should not need this step.

Then, you can build and install IMUSim with the following command::

    > python setup.py install

If all the necessary prerequisites were set up correctly, this should complete
the installation of IMUSim.

If you have problems, you can get help on the `IMUSim mailing list
<http://groups.google.com/group/imusim-users>`_.

Once you have installed IMUSim, see the :doc:`IMUSim tutorial <tutorial>` for
how to get started using it.
