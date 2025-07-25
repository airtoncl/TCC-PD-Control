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
_, lfp_healthy_graph = full(pd=0, t_sim=10000)
beta_target = extract_beta_power(lfp_healthy_graph[7])  # STN é o eletrodo 7
print(beta_target)

# 2. Início do controle
#stim_amp = 2e-3
stim_freq = 50  # Hz inicial
freq_min, freq_max = 0, 200  # limites realistas
integral = 0
prev_error = 0

# Inicialização dos ganhos base
Kp_base, Ki_base, Kd_base = 5e3, 1e3, 1e3
Kp, Ki, Kd = Kp_base, Ki_base, Kd_base

# Parâmetros de adaptação
max_amp = 3e-3
min_amp = 0.0
gain_factor = 1e-3
damping_factor = 0.5

freq_history = []
beta_history = []
error_history = []

# 3. Loop PID adaptativo
for step in range(500):
    print(f"\n--- Iteração {step+1} ---")
    
    _, lfp_pd_graph = full(pd=1, t_sim=8000, freq=stim_freq, stim_type='STN')
    beta_power = extract_beta_power(lfp_pd_graph[7])

    error = beta_power - beta_target
    delta_error = error - prev_error
    integral += error

    # Adaptação dos ganhos
    Kp = Kp_base * (1 + gain_factor * abs(error))
    Ki = Ki_base * (1 + gain_factor * abs(error))
    Kd = Kd_base * (1 + gain_factor * abs(delta_error))

    # Damping para evitar overshoot quando erro muda de sinal
    if np.sign(error) != np.sign(prev_error):
        Kp *= damping_factor
        Kd *= damping_factor

    # Controle PID
    output = Kp * error + Ki * integral + Kd * delta_error
    prev_error = error

    # Atualiza o estímulo (atuador)
    #stim_amp = np.clip(stim_amp - output, min_amp, max_amp)
    # Atualizar frequência (com saturação)
    stim_freq = np.clip(stim_freq - output, freq_min, freq_max)

    print(f"Erro: {error:.2e} | Potência Beta: {beta_power:.2e} | Freq STN: {stim_freq:.2e} A")
    print(f"Gains: Kp={Kp:.1f}, Ki={Ki:.1f}, Kd={Kd:.1f}")

    freq_history.append(stim_freq)
    beta_history.append(beta_power)
    error_history.append(error)

# 4. Plotar resultados
plt.figure(figsize=(12, 4))
plt.subplot(1, 2, 1)
plt.plot(beta_history, label='Beta Power - PD')
plt.axhline(beta_target, color='g', linestyle='--', label='Target (Healthy)')
plt.xlabel("Iterações")
plt.ylabel("Potência Beta")
plt.title("Controle Adaptativo - Potência Beta (STN)")
plt.grid(True)
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(freq_history, label='Frequência DBS (Hz)')
plt.xlabel("Iterações")
plt.ylabel("Frequência (Hz)")
plt.title("Frequência de Estímulo Adaptada")
plt.grid(True)
plt.legend()

plt.tight_layout()
plt.savefig("Controle_Frequencia_DBS.png")
plt.show()