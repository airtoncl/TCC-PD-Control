from FULL_LFP import full
from matplotlib import pyplot as plt
from scipy.signal import welch
import numpy as np

A = 1.0

def fitnessD( x, target ):
    d = min( A, ( abs( A * (x - target) ) / target ) )
    return d

def computeFitness( ref_dict, target_dict ):
    n = len( ref_dict )
    fit = 0.0
    for key in ref_dict.keys():
        fit += fitnessD( ref_dict[key] , target_dict[key] )
    print( 'fitness:' , n - fit )
    return n - fit

def extract_beta_power(signal, fs=1000):
    f, Pxx = welch(signal, fs=fs, nperseg=1024)
    beta_band = (f >= 13) & (f <= 30)
    beta_power = np.trapz(Pxx[beta_band], f[beta_band])
    return beta_power

# 1. Obter o valor de referência (Healthy)
_, lfp_healthy_graph = full(pd=0, t_sim=4000)
beta_target = extract_beta_power(lfp_healthy_graph[7])  # STN é o eletrodo 7
print(beta_target)

# 2. Configuração do PID
Kp, Ki, Kd = 1e4, 5e3, 2e4  # Ajuste fino necessário!
integral = 0
prev_error = 0
stim_amp = 2e-3  # valor inicial da corrente de estímulo
stim_history = []
beta_history = []

# 3. Loop de controle
for step in range(100):
    print(f"\n--- Iteração {step+1} ---")
    
    # Simulação com estímulo atual no STN
    _, lfp_pd_graph = full(pd=1, t_sim=4000, Input_STN_amp=stim_amp)
    beta_power = extract_beta_power(lfp_pd_graph[7])
    error = beta_power - beta_target
    integral += error
    derivative = error - prev_error

    # PID
    output = Kp*error + Ki*integral + Kd*derivative
    prev_error = error

    # Atualização da corrente de estímulo
    stim_amp = max(0.0, min(3e-3, stim_amp - output))  # saturação entre 0 e 3nA

    print(f"Erro: {error:.2e} | Potência Beta (PD): {beta_power:.2e} | Corrente STN: {stim_amp:.2e} A")

    stim_history.append(stim_amp)
    beta_history.append(beta_power)

# 4. Visualização dos resultados
plt.figure(figsize=(10,4))
plt.plot(beta_history, label='Potência Beta - PD')
plt.axhline(beta_target, color='g', linestyle='--', label='Target (Healthy)')
plt.xlabel("Iterações")
plt.ylabel("Potência Beta")
plt.title("Controle PID sobre a Potência Beta no STN")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("PID_STN_beta_control.png")
plt.show()