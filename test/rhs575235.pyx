# This file is generated automatically by QuTiP.
# (C) Paul D. Nation & J. R. Johansson

from numpy import *
cimport libc.math as cmath

import numpy as np
cimport numpy as np
cimport cython
from qutip.cy.spmatfuncs import spmv_csr, spmvpy

ctypedef np.complex128_t CTYPE_t
ctypedef np.float64_t DTYPE_t



@cython.boundscheck(False)
@cython.wraparound(False)

def cy_td_ode_rhs(double t, np.ndarray[CTYPE_t, ndim=1] vec, np.int_t hb, np.int32_t w0_A0T0, np.int_t Ohm_L0A0T0, np.int_t fL_L0, np.int_t wL_L0):
    
    cdef Py_ssize_t row
    cdef int num_rows = len(vec)
    cdef np.ndarray[CTYPE_t, ndim=1] out = np.zeros((num_rows),dtype=np.complex)
     
    return out
