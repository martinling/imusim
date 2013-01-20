"""
Quaternion maths utility module cython header

(c) University of Edinburgh 2010
"""

ctypedef struct quaternion_t:
    double w
    double x
    double y
    double z

cdef inline void mult_quat_quat(quaternion_t *q, quaternion_t *p,
        quaternion_t *dest)
cdef inline void mult_quat_scalar(quaternion_t *q, double scalar, 
        quaternion_t *dest)
cdef inline void quaternion_add(quaternion_t *q, quaternion_t *p,
        quaternion_t *dest)
cdef inline void quaternion_sub(quaternion_t *q, quaternion_t *p,
        quaternion_t *dest)
cdef inline void quaternion_log(quaternion_t *q, quaternion_t *dest)
cdef inline void quaternion_exp(quaternion_t *q, quaternion_t *dest)

cdef class Quaternion:
    cdef quaternion_t _components
