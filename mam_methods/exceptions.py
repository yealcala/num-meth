

class ResidualToleranceError(Exception):
    """Excepción lanzada cuando la solución dada contiene error significativo"""
    def __init__(self, res, mensaje="La solución no es lo suficientemente precisa"):
        self.res = res
        super().__init__(f"{mensaje}. Norma del residuo: {res}")

class ConvergenceError(Exception):
    """Excepción lanzada cuando un método no converge"""
    def __init__(self, mensaje="El método no converge"):
        super().__init__(mensaje)