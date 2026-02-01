from sage.all import *
from typing import Union
from .remonte import remonteArriba
from ...verify import nonNull

def gaussian(A: matrix, b: Union[vector, None], method: Union[str, None] = None, tol: Union[float, None] = None) -> tuple[matrix, vector]:
    """Aplica el método de eliminación gaussiana sin pivotar para resolver el sistema Au = b

        Si el argumento "pause" se omite, por defecto será Verdad.

        Parameters:
        ----------
        A : matrix
            Una matriz invertible cuyas diagonales son no nulas.
        b : vector
            El término independiente del sistema.
        method : str
            none -: No hay pivotaje.
            partial -: Pivotaje parcial.
            scaled -: Pivotaje parcial escalado.
            full -: Pivotaje total
        steps : bool, optional
            Si se debería o no mostrar cada paso por pantalla.

        Raises:
        ------
        ZeroDivisionError
            Si la matriz es singular en alguna columna.
        ValueError
            Si faltan los argumentos "A" o "b".
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

    u = remonteArriba(M.submatrix(0, 0, n, n), M.column(n), tol)
    
    if method == "full": # Undo permutations.
        u_final = [0]*n
        for idx_actual, idx_original in enumerate(p_cols):
            u_final[idx_original] = u[idx_actual]
        return (M, vector(QQ, u_final))
    
    return (M, u)