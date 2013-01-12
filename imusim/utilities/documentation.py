"""
Utilities for documenting code.
"""
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

class prepend_method_doc(object):
    """
    Decorator to prepend documentation of a matching named method.

    Usage::

        class A:
            def foo():
                \"\"\" Documentation of A.foo \"\"\"

        class B:
            @prepend_method_doc(A)
            def foo():
                \"\"\" Documentation of B.foo \"\"\"

    The documentation of B.foo will have the documentation of A.foo prepended.

    """
    def __init__(self, cls):
        self.cls = cls

    def __call__(self, f):
        methodname = f.__name__
        srcmethod = getattr(self.cls, methodname)
        srcdoc = getattr(srcmethod, '__doc__', '')
        f.__doc__ = srcdoc + f.__doc__ if f.__doc__ is not None else ''
        return f
