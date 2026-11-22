import numpy as np
import parametros as p
import scipy.linalg as sl
from scipy.signal import cont2discrete


E_DES = 2*p.m_pend*p.g*p.l_pend     # Energia do carrinho quando ele está no topo
KE = 5    # GANHO -- ajustar se ficar muito doido --

# valores de ganho: 0.22 s

# ===========================================================================================================================================
# código calcula a energia de um pêndulo a partir do estado atual, compara com uma energia desejada e gera uma ação para empurrar o carrinho 
# e bombear energia. Depois converte essa força em volts para o motor.
# ===========================================================================================================================================


def energia(x):

    theta = x[1]    # Posição angular
    theta_dot = x[3]    # Velocidade Angular

    return p.m_pend*p.g*p.l_pend*(1 - np.cos(theta)) + 1/2*(p.I_pend+p.m_pend*p.l_pend**2)*theta_dot**2



def aceleracao_desejada(x):
    
    theta = x[1]        # Posição angular
    theta_dot = x[3]    # Velocidade Angular
    diff_E = energia(x) - E_DES    # Erro da energia

    sinal = np.sign(theta_dot*np.cos(theta))        # Diz para qual lado empurar

    return KE*p.g*diff_E*sinal - 8*x[0]         # -> 8*x[0] puxa o carrinho de volta para o centro do trilho


def force_swingup(x, a_des):
    
    theta = x[1]    # posição angular
    x_dot = x[2]    # posição do carro
    theta_dot = x[3]    # velocidade angular do pendulo

    m, l, I, b = p.m_pend, p.l_pend, p.I_pend, p.b_pend
    M, c, g = p.M_car, p.C_car, p.g

    theta_2dot = (-b*theta_dot - m*l*np.cos(theta)*a_des - m*g*l*np.sin(theta)) / (I + m*l**2)  # aceleração angular

    return (M+m)*a_des + m*l*np.cos(theta)*theta_2dot - m*l*np.sin(theta)*theta_dot**2 + c*x_dot


#=========================================================
# A função recebe a força desejada e a velocidade do carrinho
#=========================================================
def force2volt(F, x_dot):
    """Converte força no carrinho [N] em tensão no motor [V]."""
    return (F * p.Rm_car * p.R_polia**2 + p.Kt_car * p.Kb_car * x_dot) / (p.Kt_car * p.R_polia)


# ==================================================================
# O swing up decide o empurrão, calcula a força para dar esse empurrão
# e traduz em volts para o motor
# ==================================================================

def swingup(x):

    a_des = aceleracao_desejada(x)  # È a aceleração que que o carrinho deve ter para bombear energia
    F = force_swingup(x, a_des)     # Aqui a aceleração vira força

    return force2volt(F, x[2])      # A força vira tensão, usando a velocidade do carrinho, que é o índice 2







# ==================================================================
# LQR discreto (params_ss_pendulum + script híbrido do autor)
# ==================================================================
Q_LQR = np.diag([500, 100, 1, 1])
R_LQR = np.array([[0.001]])


def matrizes_continuas():
    g = p.g
    m, l, I, b = p.m_pend, p.l_pend, p.I_pend, p.b_pend
    M, c = p.M_car, p.C_car
    Rm, Kb, Kt, r = p.Rm_car, p.Kb_car, p.Kt_car, p.R_polia

    alpha = (I + m*l**2)*(M + m) - (m*l)**2
    d = c + Kt*Kb/(Rm*r**2)

    Ac = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, (m**2*l**2*g)/alpha, -(I + m*l**2)*d/alpha, -b*m*l/alpha],
        [0, m*g*l*(M + m)/alpha, -m*l*d/alpha, -b*(M + m)/alpha]])

    Bc = np.array([[0],
                   [0],
                   [(I + m*l**2)*Kt/(alpha*Rm*r)],
                   [m*l*Kt/(alpha*Rm*r)]])
    return Ac, Bc


def matrizes_discretas():
    Ac, Bc = matrizes_continuas()
    # c2d do MATLAB usa segurador de ordem zero (zoh) por padrão
    A, B, _, _, _ = cont2discrete((Ac, Bc, np.eye(4), np.zeros((4, 1))), p.TS, method='zoh')
    return A, B


def dlqr(A, B, Q, R):
    P = sl.solve_discrete_are(A, B, Q, R)
    return np.linalg.solve(R + B.T @ P @ B, B.T @ P @ A)


def ganho_lqr():
    A, B = matrizes_discretas()
    return dlqr(A, B, Q_LQR, R_LQR)


# ==================================================================
# Chaveamento entre swing up e LQR (script híbrido do autor)
# ==================================================================
THETA_SWITCH = 15 * np.pi/180        # limite angular [rad]
THETA_DOT_SWITCH = 100 * np.pi/180   # limite de velocidade angular [rad/s]
X_DES = np.array([0.0, np.pi, 0.0, 0.0])


def sat(u, u_max, u_min):
    return min(max(u, u_min), u_max)


def controle_hibrido(x, K):
    if abs(x[0]) >= 0.23:
        u = -15*x[0] - 80*x[2]                          # proteção do trilho
    elif abs(x[1] - np.pi) < THETA_SWITCH and abs(x[3]) < THETA_DOT_SWITCH:
        u = -(K @ (x - X_DES))[0]                       # LQR
    else:
        u = swingup(x)                                  # swing up por energia
    return sat(u, 12, -12)


# ==================================================================
# Retorno ao centro do trilho (usado pelo botão Reiniciar)
# Lei proporcional simples; não vem do repositório do autor.
# ==================================================================
KP_CENTRO = 100    # [V/m]


def centralizar(x):
    return sat(-KP_CENTRO*x[0], 12, -12)

if __name__ == "__main__":
    A, B = matrizes_discretas()
    K = ganho_lqr()
    print("K =", K)
    print("|autovalores de A|:", np.abs(np.linalg.eigvals(A)))
    print("|autovalores de A-BK|:", np.abs(np.linalg.eigvals(A - B @ K)))
