import numpy as np
import matplotlib.pyplot as plt


def energy(string_name):
    """
    Funcion encargada de generar los graficos de energia del sistema de particulas.
    Keyword Arguments:
    string_name: str (nombre del archivo de datos).
    """
    data = np.loadtxt(string_name, dtype=float)
    time = data[:,0]
    energy_K = 0.001 * data[:,1]
    energy_U = 0.001 * data[:,2]
    total_energy = 0.001 * data[:,3]

    plt.plot(time, energy_K, label="kinetic energy")
    plt.plot(time, energy_U, label="potential energy")
    plt.plot(time, total_energy, label="total energy")
    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 0.9))
    plt.grid(True, linestyle="--", color="0.5")
    plt.xlim(xmin=time[0], xmax=time[-1])
    plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 
        ["$0$", "$1$", "$2$", "$3$", "$4$", "$5$", "$6$", "$7$", "$8$", "$9$"])
    #plt.title("Energia del Sistema", fontsize=20)
    plt.ylabel("$E$[$x10^3$]", fontsize=20)
    plt.xlabel("$t$", fontsize=20)
    plt.show()


def pressure(string_name):
    """
    Funcion encargada de generar el grafico de presion del sistema de particulas.
    Keyword Arguments:
    string_name: str (nombre del archivo de datos). 
    """
    data = np.loadtxt(string_name, dtype=float)
    time = data[:,0]
    press = data[:,4]

    plt.plot(time, press, label="presion")
    plt.legend(loc='upper right', bbox_to_anchor=(1.0, 0.9))
    plt.grid(True, linestyle="--", color="0.5")
    plt.xlim(xmin=time[0], xmax=time[-1])
    plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 
        ["$0$", "$1$", "$2$", "$3$", "$4$", "$5$", "$6$", "$7$", "$8$", "$9$"])
    plt.xlim(xmin=time[0], xmax=time[-1])
    #plt.title("Presion del Sistema", fontsize=20)
    plt.ylabel("P", fontsize=20)
    plt.xlabel("t", fontsize=20)
    plt.show()



def mean_value(vector, time):
    """
    Funcio que calcula el valor medio de vector, respecto a time.
    Keyword Arguments:
    vector: numpy.array (funcion a calcular valor medio).
    time: numpy.array (funcion con respecto a quien se toma el valor medio).
    """
    account = 0
    
    for i in range(len(time)):
        account += vector[i] / (len(time))
       
    return account



def dispertion(vector, time):
    """
    Funcion que calcula la desviacon standar de vector, respecto a time.
    Keyword Arguments:
    vector: numpy.array (funcion a calcular valor de dispersion).
    time: numpy.array (funcion con respecto a quien se toma el valor de dispersion).
    """
    auxl = [value ** 2 for value in vector]
    return np.sqrt(mean_value(auxl, time) - mean_value(vector, time) ** 2)


def structure(name_1, name_2):
	"""
    Funcion que grafica el parametro de estructura.
    Keyword Arguments:
    name_1: str (archivo 1).
    name_2: str (archivo 2).
    """
	data_1 = np.loadtxt(name_1, dtype=float)
	data_2 = np.loadtxt(name_2, dtype=float)
    
	time_1 = data_1[:,0]
	stru_1 = data_1[:,1]

	time_2 = data_2[:,0]
	stru_2 = data_2[:,1]

	plt.plot(time_1, stru_1, label="liquido")
	plt.plot(time_1, stru_1, label="solido")
	plt.grid(True, linestyle="--", color="0.5")
	plt.xscale("log")
	plt.yscale("log")
	plt.xlabel("$t$")
	plt.ylabel("$S(k)$")
	plt.show()


def graficador(date_1, date_2, defaull):
	"""
    Funcion que grafica funcion radial y coef de difucion.
    Keyword Arguments:
    name_1: numpy.array (archivo 1).
    name_2: numpy.array (archivo 2).
    """
    plt.plot(date_1, date_2, "k.")
    plt.grid(True, linestyle="--", color="0.5")
    plt.xlabel("$t$",fontsize=20)
    plt.ylabel(defaull,fontsize=20)
    plt.show()
