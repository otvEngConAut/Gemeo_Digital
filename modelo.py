import numpy as np
import parametros as p

# =================================================================================================
#   Função das derivadas:
#       Dado como o sistema está agora e a força que estou aplicando, para onde ele vai se mexer?
# =================================================================================================

def dinamica(x,f):

    theta = x[1]
    x_dot = x[2]
    theta_dot = x[3]

    m, l, I, b = p.m_pend, p.l_pend, p.I_pend, p.b_pend
    M, c, g = p.M_car, p.C_car, p.g

    # Denominador comum

    alpha = m**2*l**2*np.sin(theta)**2 + M*m*l**2 + (M+m)*I

    # ========================================
    # Aceleração do carrinho
    # ========================================
    x_2dot  = ( b*m*l*theta_dot*np.cos(theta)
                + m**2*l**2*g*np.sin(theta)*np.cos(theta)
                + (I+m*l**2)*( f-c*x_dot + m*l*np.sin(theta)*theta_dot**2))/alpha
    
    # ========================================
    # Aceleração angular
    # ========================================

    theta_2dot = -(f*m*l*np.cos(theta)
                   -c*m*l*x_dot*np.cos(theta)
                   +m**2*l**2*theta_dot**2*np.sin(theta)*np.cos(theta)
                   +(M+m)*(b*theta_dot+m*g*l*np.sin(theta)))/alpha


    # ========================================================================================================================================
    # Retorno das quatro derivadas na mesma ordem dos estados: Velocidade do carrinho, velocidade angular do pendulo, aceleração do carrinho,
    # aceleração angular
    # ========================================================================================================================================
    return np.array([x_dot, theta_dot, x_2dot, theta_2dot])





# ==================================================================
# A tensão gera corrente, mas o motor girando gera uma tensão contrária 
# proporcional à velociade, oque reduz a força. Por isso, quanto mais 
# rápido o carrinho, menos força o motor entrega.
# ==================================================================
def volt2force(v,x_dot):    # Recebe a tensão e a velocidade do carrinho
        
        return (p.Kt_car*v*p.R_polia - p.Kt_car*p.Kb_car*x_dot)/(p.Rm_car*p.R_polia**2)




# ==================================================================
# Ela recebe o estado do pendulo e a tensão do motor. Primeiro converte a
# tensão em força, usando a velocidade do carrinho naquele instante.
# Depois entrega essa força para a função dinâmica, que responde: 
# "Com esse estado e essa força, quão rápido cada coisa esta mudando."
# Então deriv responde a pergunta para onde o sistema vai se mexer agora?
# ==================================================================
def deriv(x, v):
    """Derivadas do sistema dado o estado x e a tensão v [V]."""
    f = volt2force(v, x[2])   # x[2] é a velocidade do carrinho
    return dinamica(x, f)




# ==================================================================
# A função RK4 usa essa resposta para avançar no tempo, 10 milisseguntos.
# O Euler, olharia a derivada só no começo do passo e andaria em linha reta.
# O probllema é que a derivada muda durante o passo, e o erro se acumula.

# O RK4 olha 4 vezes. O k1 é a derivada no começo do passo. O k2 é a
# derivada no meio, usando k1 para estimar como o estado estará lá.
# O k3 é outra derivada no meio agora corrigida com k2. E o k4 é a derivvada no
# fim do passo, usando k3.
# ==================================================================


def passo_rk4(x, v):
    """Avança a simulação um passo de p.TS segundos (Runge-Kutta de 4ª ordem)."""
    dt = p.TS
    k1 = deriv(x, v)
    k2 = deriv(x + dt/2 * k1, v)
    k3 = deriv(x + dt/2 * k2, v)
    k4 = deriv(x + dt * k3, v)

    return x + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

