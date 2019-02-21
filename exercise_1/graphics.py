import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def energy(string_name):
    """
    Genaración de los gráficos de energía en función del tiempo.
    """
    data = np.loadtxt(string_name, dtype=float)
    time = data[:,0]
    energy_K = 0.01*data[:,1]
    energy_U = 0.01*data[:,2]
    total_energy = 0.01*data[:,3]

    plt.plot(time, energy_K, label="kinetic energy")
    plt.plot(time, energy_U, label="potential energy")
    plt.plot(time, total_energy, label="total energy")
    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 0.9))
    plt.grid(True, linestyle="--", color="0.5")
    plt.xlim(xmin=time[0], xmax=time[-1])
    #plt.title("Energía del Sistema", fontsize=20)
    plt.ylabel("$E$[$x10^2$]", fontsize=20)
    plt.xlabel("$t$", fontsize=20)
    plt.show()

def histogram(string_name):
    """
    Genaración del histograma de velocidades.
    """
    archive = open(string_name, "r")
    data = np.array([float(value.replace("\n", "")) for value in archive.readlines()])
    
    media = data.mean()
    sigma = data.std()

    num_bins = 19

    n, bins, patches = plt.hist(data, num_bins, normed=1, facecolor="blue", alpha=1.0, edgecolor="black", 
		label=" $\sigma$ = {0:.3f}".format(sigma))
    
    y = mlab.normpdf(bins, media, sigma)
    
    plt.plot(bins, y, "k--")
    plt.xlabel('$|v|$', fontsize=20)
    plt.ylabel('$P(|v|)/N$', fontsize=20)
    #plt.title("Histograma", fontsize=20)
    plt.grid(True, linestyle="--", color="0.5")
    plt.legend()
    plt.subplots_adjust(left=0.15)
    plt.show()

def main():
    energy("Energy")
    histogram("histogr_1")

if __name__ == '__main__':
    main()
