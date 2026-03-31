# num_meth/methods/linsys/all.py

"""Provides all methods regarding linear systems.

The module imports the following modules:

- `factorization` - Provides methods for factorizing a matrix.
- `gauss` - Provides a method for echelonizing and solving using Gauss.
- `remonte` - Provides methods for solving by Remonte's methods.
- `iterative` - Provides iterative methods of Jacobi, Gauss-Seidel and Relaxation.
- `other` - Provides some other methods.
"""

from .factorization import *
from .gauss import *
from .remonte import *
from .iterative import *
from .other import *