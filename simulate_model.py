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

# pd=0, it_num=1, t_sim=100.000
# freq=0, cortstim=0, str_stim=False
# std=0, stim_type = 'None', faixa = 'high'
for i in range(1,2): #range(1,21)
    lfp_pd, lfp_pd_graph = full(pd=1, it_num=i, t_sim=2000) # parkinsonian without DBS

    lfp_healthy, lfp_healthy_graph = full(pd=0, it_num=i, t_sim=2000) # healthy without DBS

    fitness =  computeFitness( ref_dict=lfp_pd, target_dict=lfp_healthy )

    idx = 7  # índice do STN

    print('PD:: \n', lfp_pd)

    print('Healthy:: \n', lfp_healthy)

     # Plot do sinal bruto
    plt.figure(figsize=(10, 3))
    plt.plot(lfp_pd_graph[idx], color='r', label='pd')
    plt.plot(lfp_healthy_graph[idx], color='g', label='healthy')
    plt.title("LFP - Eletrodo 0")
    plt.xlabel("Tempo (passos)")
    plt.ylabel("Potencial (uV)")
    plt.grid(True)
    plt.tight_layout()
    plt.legend()
    plt.savefig('Comparacao.png')

    

    fs = 1000  # taxa de amostragem (Hz)
    f_pd, Pxx_pd = welch(lfp_pd_graph[idx], fs=fs, nperseg=1024)
    f_healthy, Pxx_healthy = welch(lfp_healthy_graph[idx], fs=fs, nperseg=1024)

    plt.figure(figsize=(10, 4))
    plt.semilogy(f_pd, Pxx_pd, color='r', label='PD')
    plt.semilogy(f_healthy, Pxx_healthy, color='g', label='Healthy')
    plt.title("Espectro de Potência - LFP (Eletrodo 0)")
    plt.xlabel("Frequência (Hz)")
    plt.ylabel("Potência (uV²/Hz)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.xlim(xmin=0, xmax=50)
    plt.savefig("Espectro_Potencia_Comparacao.png")

    beta_power_pd = np.trapz(Pxx_pd[(f_pd >= 13) & (f_pd <= 30)], f_pd[(f_pd >= 13) & (f_pd <= 30)])
    beta_power_healthy = np.trapz(Pxx_healthy[(f_healthy >= 13) & (f_healthy <= 30)], f_healthy[(f_healthy >= 13) & (f_healthy <= 30)])

    print(f"Potência Beta no STN - PD: {beta_power_pd:.2e}")
    print(f"Potência Beta no STN - Healthy: {beta_power_healthy:.2e}")

