# cython: experimental_cpp_class_def=True
# tag: cpp

import cython
from libcpp cimport bool

# import both numpy and the Cython declarations for numpy
import numpy as np
cimport numpy as np

# Declaring a C++ class interface
cdef extern from "dispersion.h":
    cdef cppclass Disperser:

        # constructor
        Disperser(double, double, double) except +

        # methods
        double lower_y, upper_y
        double lower_y, upper_y
        double get_f()
        double get_u(double)
        double get_v(double)
        double get_l(double)
        void update(int, double*, double*, double, int)
        
        # members
        double u_f, v_m, displacement
        bool bouncing_particles, direction_likelyhood, bursting_process, heavy
        bool gaussian_u, gaussian_v
        double fall_vel

# Create Cython wrapper class:
# makes this accessible from external Python code
cdef class PyDisperser:

    # hold a C++ instance which we're wrapping
    cdef Disperser *thisptr

    # constructor
    def __cinit__(self, Re, lower_y, upper_y): 
        self.thisptr = new Disperser(Re, lower_y, upper_y)

    #  destructor
    def __dealloc__(self):      
        del self.thisptr

    # methods

    def get_f(self):
        return self.thisptr.get_f()

    def get_u(self, y):
        return self.thisptr.get_u(y)

    def get_v(self, y):
        return self.thisptr.get_v(y)

    def get_l(self, y):
        return self.thisptr.get_l(y)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def update(self, N, 
                    np.ndarray[double, ndim=1, mode="c"] X not None, 
                    np.ndarray[double, ndim=1, mode="c"] Y not None,
                    tF, max_step):
        self.thisptr.update(N, &X[0], &Y[0], tF, max_step)
        return None

    # Properties
    property u_f:
        def __get__(self): return self.thisptr.u_f
        def __set__(self, u_f): self.thisptr.u_f = u_f

    property v_m:
        def __get__(self): return self.thisptr.v_m
        def __set__(self, v_m): self.thisptr.v_m = v_m

    property bouncing_particles:
        def __get__(self): return self.thisptr.bouncing_particles
        def __set__(self, bouncing_particles): 
            self.thisptr.bouncing_particles = bouncing_particles

    property direction_likelyhood:
        def __get__(self): return self.thisptr.direction_likelyhood
        def __set__(self, direction_likelyhood): 
            self.thisptr.direction_likelyhood = direction_likelyhood

    property bursting_process:
        def __get__(self): return self.thisptr.bursting_process
        def __set__(self, bursting_process): 
            self.thisptr.bursting_process = bursting_process

    property gaussian_u:
        def __get__(self): return self.thisptr.gaussian_u
        def __set__(self, gaussian_u): 
            self.thisptr.gaussian_u = gaussian_u

    property gaussian_v:
        def __get__(self): return self.thisptr.gaussian_v
        def __set__(self, gaussian_v): 
            self.thisptr.gaussian_v = gaussian_v

    property displacement:
        def __get__(self): return self.thisptr.displacement
        def __set__(self, displacement):
            self.thisptr.displacement = displacement

    property fall_vel:
        def __get__(self): return self.thisptr.fall_vel
        def __set__(self, fall_vel):
            self.thisptr.fall_vel = fall_vel
