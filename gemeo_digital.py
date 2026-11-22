import numpy as np
import modelo_gemeo


class GemeoDigital:

    def __init__(self, x0):

        self.x = np.array(x0, dtype=float)

    def ler(self):

        return self.x.copy()

    def atualizar(self, u):

        self.x = modelo_gemeo.passo_rk4(self.x, u)