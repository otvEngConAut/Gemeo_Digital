import numpy as np
import modelo


class PlantaFalsa:
    """Finge ser a bancada (ESP32 + carrinho + pêndulo)."""

    def __init__(self, x0):
        self.x = np.array(x0, dtype=float)
        self.u = 0.0

    def ler(self):
        """Estado medido: [posição, ângulo, vel. do carrinho, vel. angular]."""
        return self.x.copy()

    def escrever(self, u):
        """Recebe a tensão [V] e avança a bancada um passo de TS."""
        self.u = u
        self.x = modelo.passo_rk4(self.x, u)