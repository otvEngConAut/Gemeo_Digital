import numpy as np
import parametros as p
import controle

from planta import PlantaFalsa
from gemeo_digital import GemeoDigital


# =========================================================
# CONFIGURAÇÃO DA SIMULAÇÃO
# =========================================================

TEMPO_SIMULACAO = 30.0

X0 = np.array([0.0, 0.05, 0.0, 0.0])


# =========================================================
# INICIALIZAÇÃO
# =========================================================

def criar_simulacao():

    K = controle.ganho_lqr()

    N = round(TEMPO_SIMULACAO / p.TS)

    planta = PlantaFalsa(X0.copy())

    gemeo = GemeoDigital(X0.copy())

    return K, N, planta, gemeo


# =========================================================
# EXECUÇÃO DE UM PASSO DA SIMULAÇÃO
# =========================================================

def executar_passo(planta, gemeo, K, controlador_ativo=True):

    x = planta.ler()

    if controlador_ativo:
        u = controle.controle_hibrido(x, K)
    else:
        u = 0.0          # sem controle: tensão zero

    planta.escrever(u)
    gemeo.atualizar(u)

    x_planta = planta.ler()
    x_gemeo = gemeo.ler()

    return x_planta, x_gemeo, u

def executar_passo_centralizar(planta, gemeo):
    """Um passo levando o carro ao centro do trilho (sem controle do pêndulo)."""
    x = planta.ler()
    u = controle.centralizar(x)

    planta.escrever(u)
    gemeo.atualizar(u)     # o gêmeo recebe a mesma tensão

    return planta.ler(), gemeo.ler(), u