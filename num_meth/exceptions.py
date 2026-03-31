

class ResidualToleranceError(Exception):
    """This error is thrown where solution is not accurate or precise enough"""
    def __init__(self, res, mensaje="Solution is not accurate enough"):
        self.res = res
        super().__init__(f"{mensaje}. Norma del residuo: {res}")

class ConvergenceError(Exception):
    """This error is thrown if a method doesn't converge."""
    def __init__(self, mensaje="Method doesn't converge"):
        super().__init__(mensaje)