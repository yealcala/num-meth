from sage.all import *
from typing import Optional
from ...utils import pnorm
from ...exceptions import ConvergenceError
from .remonte import remonteAbajo
from ...verify import nonNull, matrixMustBeSquare, matrixMustBeInvertible, matrixMustHaveNonNullDiagonal, vectorMustBeDimension

def __DEF(A: matrix) -> tuple:
    n = A.nrows()
    D = matrix.zero(QQ, nrows=n, ncols=n)
    E = matrix.zero(QQ, nrows=n, ncols=n)
    F = matrix.zero(QQ, nrows=n, ncols=n)
    for i in range(0, n):
        for j in range(0, n):
            if i == j: # Diagonal
                D[i,j] = A[i,j]
            elif i > j: # Triangular inferior
                E[i,j] = -A[i,j]
            else:
                F[i,j] = -A[i,j]
    return (D,E,F)


def __powerIteration(M: matrix, N: matrix, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo):    
    z = z0.change_ring(M.base_ring())/pnorm(z0, norm)

    lOld = 0
    k = 0
    error = tol + 1
    while error > tol and k < itMax:
        y = remonteAbajo(M, N*z)
        zOld = z

        normaY = pnorm(y, norm)
        if normaY == 0: return (0, z, k)
        z = y / normaY

        j = 0
        maxVal = -1
        for i in range(len(z)):
            if abs(z[i]) > maxVal:
                maxVal = abs(z[i])
                j = i

        lCurr = y[j] / zOld[j]

        error = abs(lCurr - lOld)
        lOld = lCurr
        k += 1


    return (lOld, z, k)

def jacobiIteration(A: matrix, b: vector, x0: Optional[vector] = None, tol: float = 1e-10, maxIt: int = 100) -> tuple[vector, float, int]:
    nonNull(A, b, -1, tol, maxIt)
    matrixMustBeSquare(A)
    matrixMustBeInvertible(A)
    matrixMustHaveNonNullDiagonal(A)
    n = A.ncols()
    vectorMustBeDimension(b, n)
    if x0 is not None: vectorMustBeDimension(x0, n)


    (D, E, F) = __DEF(A)
    M = D
    N = E+F

    rho, _, _ = __powerIteration(M, N, vector(A.base_ring(), [1]*A.ncols()), itMax=50)
    if abs(rho) >= 1:
        raise ConvergenceError("Jacobi method has no convergence for input matrix.")
    
    
    xCurr = x0 if x0 is not None else vector(A.base_ring(), [0]*n)
    k = 0
    error = tol + 1
    while error > tol and k < maxIt:
        xNext = remonteAbajo(M, N*xCurr + b)
        error = float(pnorm(xNext - xCurr, oo))
        xCurr = xNext
        k += 1
    
    return (xCurr, error, k)




def __getRho(w: float, D: matrix, E: matrix, F: matrix, n: int) -> float:
    """ Calcula el radio espectral de la matriz de iteración para un w dado. """
    try:
        M = (D/w - E)
        N = F + (1-w)/w * D
        

        v_init = vector(D.base_ring(), [1]*n)
        rho, _, _ = __powerIteration(M, N, v_init, itMax=100)
        
        return float(abs(rho))
    except:
        return 1.0

def __findBestw(D: matrix, E: matrix, F: matrix, n: int) -> float:
    """ 
        Busca el w óptimo mediante una búsqueda ternaria.
    """
    # El límite superior se baja ligeramente de 1.99 a 1.85 para asegurar.
    #   Si queda muy cerca de 2, la convergencia suele ser muchísimo peor
    #   debido a la inestabilidad númerica de la matriz de Sobrerelajación.
    izq = 0.01
    der = 1.85
    
    for _ in range(50):
        m1 = izq + (der - izq) / 3
        m2 = der - (der - izq) / 3
        
        rho1 = abs(__getRho(m1, D, E, F, n))
        rho2 = abs(__getRho(m2, D, E, F, n))
        
        if rho1 < rho2:
            der = m2
        else:
            izq = m1
    w_opt = (izq + der) / 2
   
    if __getRho(w_opt, D, E, F, n) >= 1.0:
        show("Wopts es muy grande")
        return 1.0 # Si wopt es "muy grande", mejor Gauss-Seidel...
        
    return w_opt


_defaultOptions = {"tolerance": 1e-10, "max_iterations": 100, "w": None, "norm": oo}

def gaussSeidelIteration(A: matrix, b: vector, x0: Optional[vector], opts: Optional[dict[str, any]] = _defaultOptions) -> tuple[vector, float, int, float]:
    if opts is None: opts = {"w": 1}
    return relaxIteration(A, b, x0, opts)


def relaxIteration(A: matrix, b: vector, x0: Optional[vector], opts: Optional[dict[str, any]] = _defaultOptions) -> tuple[vector, float, int, float]:
    """
        Calcula iterativamente la solución del sistema lineal Ax = b mediante el Método de Sobrerelajación
            calculando para ello, si fuera necesario, el valor óptimo de w.
        
        Precondiciones: El invocador debe asegurarse de que el método iterativo de sobrerelajación asociado a ese
            sistema lineal sea convergente. Además, si hay alguna opción declarada en opts es incorrecta, el mal funcionamiento 
            del programa es completamente responsabilidad del invocador.
        
        Parámetros:
        - A (matrix): La matriz del sistema.
        - b (vector): El vector del sistema.
        - x0 (vector): El vector sobre el cual empezar a iterar. Default: Vector nulo.
        - opts (dict): Opciones que se pueden añadir al método
        -> tolerance (number): El error que estamos dispuestos a asumir. Default: 10e-10
        -> max_iterations (number): El número máximo de iteraciones. Default: 200
        -> w (number): El w a usar. Default: None (el óptimo es calculado)
        -> norm (int | oo): La norma a usar. Default: oo (norma infinito).

        Return (vector): La solución aproximada al sistema usando los parámetros introducidos.

        Raises:
        - ValueError: Si algun parámetro obligatorio es incorrecto.
    """

    # <----- Comprobaciones iniciales ----->
    nonNull(A, b)
    matrixMustBeSquare(A)
    matrixMustBeInvertible(A)
    matrixMustHaveNonNullDiagonal(A)

    n = A.ncols()
    vectorMustBeDimension(b, n)
    if x0 is not None: vectorMustBeDimension(x0, n)


    # <----- Gestión de opciones ----->
    if opts is None: opts = {}
    current_opts = _defaultOptions.copy()
    current_opts.update(opts)
        
    (D, E, F) = __DEF(A)

    # <----- Establecer el w óptimo ----->
    wopt = current_opts["w"]
    if wopt is None:
        wopt = __findBestw(D, E, F, n)
    
    if wopt <= 0 or wopt >= 2:
        raise ConvergenceError(f"Method won't converge because w is outside (0, 2) interval")
    elif wopt <= 0 or wopt >= 2:
        wopt = 1

    # <---- Método de Sobrerelajación ----->
    M = (D/wopt - E)
    N = (F + (1-wopt)/wopt * D)

    if wopt == 1: # Gauss-Seidel
        rho, _, _ = __powerIteration(M, N, vector(D.base_ring(), [1]*A.ncols()), itMax=50)
        if abs(rho) >= 1:
            raise ConvergenceError("Gauss-Seidel or Relaxation method (w = 1) doesn't converge!")
    
        

    if x0 is None: x0 = vector(A.base_ring(), [0]*n)

    xCurr = x0
    k = 0
    error = current_opts["tolerance"] + 1
    while error > current_opts["tolerance"] and k < current_opts["max_iterations"]:
        xNext = remonteAbajo(M, N*xCurr + b)
        error = pnorm(xNext - xCurr, current_opts["norm"])
        xCurr = xNext
        k += 1

    return (xCurr, error, k, wopt)