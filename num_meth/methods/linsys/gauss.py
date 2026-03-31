# num_meth/methods/linsys/gauss.py

"""Provides the gauss method.

The module contains the following functions:

- `gaussian(A, b, method, tol)` - Applies gaussian method to linear system Au = b.
"""

from sage.all import matrix, vector, QQ
from typing import Union
from .remonte import upperBackSubs
from ...verify import nonNull

def gaussian(A: matrix, b: Union[vector, None] = None, method: Union[str, None] = None, tol: Union[float, None] = None) -> tuple[matrix, vector]:
    
    """ Gaussian method for making A be upper triangular. If b is given, linear system is solve by this way.

    Examples:
        >>> A = matrix(QQ, [[1, 1, 0, 0], [2, -1, 5, 0], [0, 3, -4, 2], [0, 0, 2, 6]])
        >>> b = vector(QQ, [5, -9, 19, 2])

        >>> try:
        >>>    show("No pivot:", gaussian(A, b, method=None, tol=None))
        >>>    show("Partial pivot:", gaussian(A, b, method="partial", tol=None))
        >>>    show("Partial scaled pivot:", gaussian(A, None, method="scaled"))
        >>>    show("Full pivot:", gaussian(A, b, method="full", tol=None))
        >>> except ZeroDivisionError as e:
        >>>     print(e)
        
    Args:
        A (matrix): An invertible square n-matrix. 
        b (Optional[vector]): A n-vector.
        method (Optional[string]): How will pivots be taken:
            - None (or other stuff): Without pivot. Default.
            - "partial": Partial pivot.
            - "scaled": Partial scaled pivot.
            - "full": Full pivot..
        tol (Optional[float]): Checks if given solution doesn't exceed this tolerance. Default: 10^{-10}

    Returns:
        tuple: A tuple with two components,
            - 0 (matrix): The upper triangular matrix after applying gaussian method to A.
            - 1 (vector): Linear system's solution. If b is None, then this is None.

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If some non-evitable 0 appears as a pivot.
        ResidualToleranceError: If real solution and calculated solution's difference exceeds tolerance.
    """
    nonNull(A)

    n = A.nrows()

    if b is None:
        bUse = vector(QQ, [0]*n)
    else:
        bUse = b

    M = A.augment(bUse, subdivide=True).change_ring(QQ)

    method = method.lower() if method is not None else None
    if method == "scaled":
        escalas = [max(map(abs, row)) for row in A]
    
    if method == "full":
        p_cols = list(range(n))

    for i in range(n):
        if method == "partial":
            piv = i
            numPivote = abs(M[i][i])

            for k in range(i+1, n): # Find pivote
                if(abs(M[k][i]) > numPivote):
                    piv = k
                    numPivote = abs(M[k][i])

            if(piv != i): # Permute
                M.swap_rows(i, piv)
                
        elif method == "scaled":
            piv = i
            numPivote = abs(M[i][i]/escalas[i])
            for k in range(i+1, n):
                posNumPiv = abs(M[k][i]/escalas[k])
                if(posNumPiv > numPivote):
                    piv = k
                    numPivote = posNumPiv
            if(piv != i): # Permuta
                M.swap_rows(i, piv)
                escalas[i], escalas[piv] = escalas[piv], escalas[i]
        
        elif method == "full":
            piv_i, piv_j = i, i
            numPivote = abs(M[i, i])
            for i2 in range(i, n):
                for j2 in range(i, n):
                    if abs(M[i2, j2]) > numPivote:
                        numPivote = abs(M[i2, j2])
                        piv_i, piv_j = i2, j2
            if piv_i != i:
                M.swap_rows(i, piv_i)
            if piv_j != i:
                M.swap_columns(i, piv_j)
                p_cols[i], p_cols[piv_j] = p_cols[piv_j], p_cols[i]
                

        if M[i, i] == 0:
            if b is None:
                return (M.submatrix(0, 0, n, n), None)
            raise ZeroDivisionError(f"Zero detected on (M[{i},{i}])! Solution might not be unique!")
        
        for isig in range(i+1, n): # Elimination.

            mult = - M[isig][i] / M[i][i]
            
            M[isig, i] = 0
            for j in range(i+1, n+1): 
                M[isig, j] = M[i, j]*mult + M[isig, j]
    
    if b is None:
        return (M.submatrix(0, 0, n, n), None)

    u = upperBackSubs(M.submatrix(0, 0, n, n), M.column(n), tol)
    
    if method == "full": # Undo permutations.
        u_final = [0]*n
        for idx_actual, idx_original in enumerate(p_cols):
            u_final[idx_original] = u[idx_actual]
        return (M, vector(QQ, u_final))
    
    return (M, u)