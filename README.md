* Introduction
Monte Carlo method for studying the behavior of a dispersant inserted
into a steady boundary layer channel flow. The approach is the
one-particle analysis, where the random path of a single particle
through the flow is simulated in a given time interval. The particle
is subjected to flow mechanisms in a one-way coupled manner including
vertical turbulent diffusion, ejections and boundary interactions. The
one-particle analysis is repeated a number of times and the
statistical properties, ensemble mean and variance, are calculated
based on the complete set of outcomes.
The model is implemented in C++ and Python because of poor runtime
performance for this kind of models if the code is pure Python or
Matlab. The implementation is quite efficient and bears no comparison
with pure Python or Matlab.
* The mixing of C++ and Python through Cython
The C++ class "Disperser" defined in "dispersion.h"/"dispersion.cpp"
is made available for easy interaction in Python using wrapper
functions which are defined in the Cython file "mymod.pyx".

To process it use the command:
```
cython --cplus mymod.pyx
```

Compile the C++ code: `code`
```
python setup.py build_ext --inplace
```
* Requirements
libboost-python