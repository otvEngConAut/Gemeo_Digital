# Declaração das variáveis do pêndulo invertido


# ==================================================
# DADOS DO PÊNDULO
# ==================================================

g = 9.81    # [m/s²] aceleração da gravidade

m_pend = 76.4/1000   # massa [kg]
l_pend = 0.18  # centro_de_gravidade [m] 
I_pend = 0.0012    # inercia  [kg*m²]
b_pend = 0   # coeficiente_amortecimento_viscoso [N*m*s] 
guia = 0.60 # [m] tamanho da guia linear


# ==================================================
# DADOS DO CARRO
# ==================================================
M_car = 370.2/1000  # [kg] Massa do Carro
C_car = 8.0060       # [N*s/m] Coeficinete de amortecimento viscoso entre carro e guia


# ==================================================
# DADOS DO MOTOR
# ==================================================

Rm_car= 10.6919    # [homs] resistencia do enrolamento
Kb_car = 0.2842     # [V*s / rad] Constante de força contra-eletromotriz
Kt_car = 0.1336     # [N*m/A] Constante de torque (orque/Corrente)
R_polia = 0.0198      # [m] raio da polia em 

# ==================================================
# ATRITO SECO
# ==================================================

Fc = 0.035      # [N]

# ==================================================
# SIMULAÇÃO
# ==================================================

TS = 0.01    # [s] de exportação
TF = 30     # [s] tempo total da simulação
