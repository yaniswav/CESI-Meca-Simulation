## Importation des bibliothèques nécessaires au calcul
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib import cm
from pylab import *

class Simulation:

    def __init__(self):
        
        self.m = 0.03  #Masse de la voiture en kg
        self.hp = 0.93  #Hauteur maximale de la pente en mètres
        self.alpha = 40   #Angle de la pente en degrés  
        self.alpharad = (self.alpha*np.pi)/180   #Angle de la pente en radians  
        self.hl = 0.23  #Hauteur du looping en mètres 
        self.r = self.hl/2   #Rayon du looping en mètres
        self.hr = 0.1   #Hauteur du ravin en mètres 
        self.lr = 0.7   #Largeur du ravin en mètres 
        self.g = 9.81   # constante gravitationnelle en m.s^-2
        self.mu = 0.002  # coefficient de roulement plastique / plastique 
        self.Cx = 0.04   # coefficient de traînée 
        self.S1 = 3*10**-4  # Surface de référence de la voiture en m^2
        self.pair = 1.225   # masse volumique de l'air en kg.m^-3
        self.S1Cx = 0.001
        self.S2Cz = 0.01

    #Le calcul de la vitesse en bas de la pente sans puis avec frottements du sol
    def vitesse_pente(self):
        self.baspente_sf = np.sqrt(2*self.g*self.hp)
        print(f"La vitesse en bas de la pente sans frottements sera de {round(self.baspente_sf,4)} m/s.")

        self.baspente_af = np.sqrt(2*self.hp*self.g*(1-(self.mu/np.tan(self.alpharad))))
        print(f"La vitesse en bas de la pente avec frottements sera de {round(self.baspente_af, 4)} m/s.")

    #Le calcul et le tracé de la vitesse de sortie en fonction de la hauteur de la pente sans et avec frottements
    def graphique_pente(self, hauteur):

        vitesse = [np.sqrt(i*self.g*2) for i in hauteur]
        plt.title("Vitesse en fonction de la hauteur sans frottements")
        plt.plot(hauteur, vitesse, color="red")
        plt.xlabel('Hauteur h (m)')
        plt.ylabel('Vitesse de la pente v (m/s)')
        plt.plot(93,4.271, marker="o", color="red")
        plt.show()

        vitesse = [np.sqrt(2*i*self.g*(1-(self.mu/np.tan(self.alpharad)))) for i in hauteur]
        plt.title("Vitesse en fonction de la hauteur avec frottements")
        plt.plot(hauteur, vitesse, color="green")
        plt.xlabel('Hauteur h (cm)')
        plt.ylabel('Vitesse de la pente v (m/s)')
        plt.plot(93,4.271, marker="o", color="red")
        plt.show()


    def vitesse_minimale_looping(self):
        #sans frottement => racine de 5gr

        self.v_min_loop_f = np.sqrt(5*self.g*self.r)
        print(f"La vitesse minimale en bas de la pente sans frottements sera de {round(self.v_min_loop_f, 4)} m/s.")

        self.v_min_loop_af = np.sqrt(self.r*self.g*(5+2*self.mu*self.alpharad))
        print(f"La vitesse minimale en bas de la pente avec frottements sera de {round(self.v_min_loop_af, 4)} m/s.")

    def equation_looping(self,y,t):
        theta, thetap, thetapp = y[0],y[1],y[2]
        dthetadt = [thetapp,thetap,-((self.pair*self.Cx*self.S1*(thetap*self.r)**2)/(2*self.m))-(self.mu*(self.r*thetap**(2)+self.g*np.cos(theta)))-(self.g*np.sin(theta))-self.r*thetapp,]
        return dthetadt

    def vitesse_minimal_ravin(self):
            self.vitesse_minimal_ravin_sf = np.sqrt((self.g*self.lr**2)/(2*self.hr))
            print(f"La vitesse minimale pour passer le ravin sans frottement est de {self.vitesse_minimal_ravin_sf}")

    def trace_ravin(self):

        # ravin
        V = np.sqrt((self.g*self.lr**2)/(2*self.hr))
        tr = np.linspace(0,0.5,100)
        # sans frottements

        xravin = V*tr
        yravin = -0.5*self.g*tr**2

        # avec frottements
        k1 = 0.5*self.pair*self.S1Cx
        k2 = 0.5*self.pair*self.S2Cz

        def afravin(vfravin, tr):
            
            df = [vfravin[2],vfravin[3], -((k1)/(self.m))*np.sqrt(vfravin[2]**2+vfravin[3]**2)*vfravin[2], -((k2)/(self.m))*np.sqrt(vfravin[2]**2+vfravin[3]**2)*vfravin[2]-self.g]
            return df

        vfravin0 = [0,0,V,0]
        vfravin = odeint(afravin, vfravin0, tr)

        # affichage du graphique
        plt.plot(xravin, yravin, color = "blue")
        plt.plot(vfravin[:,0], vfravin[:,1], color = "red")
        plt.title("Trajectoire dans le ravin")
        plt.xlabel("x(m)")
        plt.ylabel("y(m)")
        plt.show()




simulation = Simulation()

simulation.trace_ravin()
simulation.vitesse_minimal_ravin()
t = np.linspace(0,100,100)
y0 = [0,(37*np.pi)/180,(159*np.pi)/180]
sol = odeint(simulation.equation_looping,y0,t)

plt.title('Position angulaire de theta')
plt.plot(t, sol[:, 0], 'b', label='thetapp(t)')
plt.plot(t, sol[:, 1], 'r', label='thetap(t)')
plt.plot(t, sol[:, 2], 'g', label='theta(t)')
plt.legend(loc='best')
plt.xlabel('Temps')
plt.ylabel('Position en degré')
plt.grid()
plt.show()

print("Projet créé par Yanis")
# simulation.vitesse_pente()
# simulation.graphique_pente(hauteur = [i*(10**-2) for i in range(0, 96, 5)])
# simulation.vitesse_minimale_looping()
# simulation.trace_ravin()
