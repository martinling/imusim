"""
Test utilities for vectors.
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

import numpy as np

def assert_vectors_correlated(actual, desired, targetCorrelation=0.95):
    """
    Assert that the correlation between two column vector arrays is greater
    than a specified target.

    Correlations calculated as the diagonal of the upper right quadrant of
    the full correlation matrix. Each of the diagonal elements of this
    sub-matrix represents the correlation between an actual vector component
    and its corresponding desired value.

    @param actual: Actual column vector array.
    @param desired: Desired column vector array.
    @param targetCorrelation: Minimum correlation to accept.
    """
    correlationMatrix = np.corrcoef(actual, desired)
    s = correlationMatrix.shape[0] // 2
    correlationMatrix = correlationMatrix[s:, :s]
    correlations = np.diag(correlationMatrix)

    assert np.all(correlations > targetCorrelation), "Vectors are not " \
    "sufficiently correlated. \nCorrelation values: %s\nTarget %.2f"\
        %(correlations, targetCorrelation)
