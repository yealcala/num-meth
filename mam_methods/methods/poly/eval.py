from sage.all import *
from ...utils import extract_coefficients
from ...verify import nonNull

def horner(p, x: float) -> tuple[float, callable]:
    nonNull(p, x)
    coeffs = extract_coefficients(p)
    n = len(coeffs) 
    x = RR(x)
    

    Q = [RR(coeffs[0])]

    pVal = RR(coeffs[0])
    
    for i in range(1, n):
            
        pVal = RR(coeffs[i]) + x * pVal
        if i < n - 1:
            Q.append(pVal)


    if n == 1:
        return pVal, []
        
    return pVal, Q