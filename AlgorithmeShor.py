# -*- coding: utf-8 -*-
"""
Created on Sun May 10 22:25:39 2020

@author: Hugo
"""

from qiskit import QuantumCircuit, execute, Aer, IBMQ
from qiskit.compiler import transpile, assemble
from qiskit.tools.jupyter import *
from qiskit.visualization import *
import timeit
import math
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
import numpy as np
from qiskit.tools.visualization import plot_histogram
# Loading your IBMQ account(s)
provider = IBMQ.load_account()

#----------- Importation des différentes bibliothèques utiles dans le code -----------

def pgcd(x,y):
    while(y!=0):
        t = y
        y = x % y
        x = t
    return x

# Cette fonction sert à calculer le pgcd entre x et y via l'algorithme d'Euclide

def shorAlgorithm(a):
    #----------- On cherche dans cette fonction a créer un circuit adapté à "a" afin de pouvoir calculer
    # les résultats de la fonction f -----------
        
    n=3
    
    p = QuantumRegister(4,"p")
    w = QuantumRegister(n,"w")
    c = ClassicalRegister(n, "c")
    qc = QuantumCircuit(w,p, c)
    
    # On applique la fonction : "a**x mod 15"
    if (a == 2):
        qc.x(p[0])
        qc.h(w[0])
        qc.h(w[0])
        qc.measure(w[0],c[0])
        
        qc.h(w[1])
        qc.cx(w[1],p[0])
        qc.cx(w[1],p[2])
        
        if (c[0]):
            qc.u1(math.pi/2,w[0])
        qc.h(w[1])
        qc.measure(w[1],c[1])
        
        qc.h(w[2])
        
        qc.cx(p[0],p[3])
        qc.cx(p[3],p[0])
        qc.cx(p[0],p[3])
        
        qc.cx(p[1],p[0])
        qc.cx(p[0],p[1])
        qc.cx(p[1],p[0])
        
        qc.cx(p[1],p[2])
        qc.cx(p[2],p[1])
        qc.cx(p[1],p[2])
        
        if (c[1]):
            qc.u1(math.pi/2,w[2])
        if (c[0]):
            qc.u1(math.pi/4,w[2])
        
        qc.h(w[2])
        qc.measure(w[2],c[2])
        
    
    elif (a == 7):
        qc.x(p[0])
        qc.h(w[0])
        qc.h(w[1])
        qc.cx(w[0], p[1])
        qc.cx(w[0], p[2])

        qc.x(p[2])
        qc.ccx(w[1], p[2], p[0])
        qc.cx(p[0], p[2])

        qc.ccx(w[1], p[1], p[3])

    elif (a == 11):
        qc.h(w[0])
        qc.h(w[1])
        qc.cx(w[0], p[1])
        qc.cx(w[0], p[3])
    if (a in [7,11]):
        #----------- On effectue la transformée de Fourier quantique -----------
        for k in range(n):
            j = n - k
            if (j -  1) !=2:
                qc.h(w[j-1])

            for i in reversed(range(j-1)):
                qc.cu1(2*np.pi/2**(j-i),w[i], w[j-1])

        #----------- On mesure les valeurs -----------
        qc.measure(w[0],c[0])
        qc.measure(w[1],c[1])
        qc.measure(w[2],c[2])
    return qc
    
for a in [7,11]: # On teste ici toutes les valeurs de a possible pour cet implémentation
    
    start = timeit.default_timer() # Début du calcul du temps
    N = 15 # Nombre maximal à factoriser
    print("Le nombre choisi est :",a,"pour N =",N)

    #----------- Début de l'algorithme de Shor -----------
    
    gcd = pgcd(N,a) # Calcul du pgcd afin de déterminer facilement si N et a sont multiples

    if (gcd != 1): # Si N et a sont multiples alors on a directement les diviseurs
        diviseur1 = gcd
        diviseur2 = N/gcd
        print(diviseur1,diviseur2)
    else:
        qc = shorAlgorithm(a) # On exécute la fonction shor afin de récupérer le circuit 
                                   # pour mesurer le rayon de shor
            
        #----------- Exécution sur le seul processeur quantique possible car il possède 15 qubits ---------
        my_provider = IBMQ.get_provider()
        backend  = my_provider.get_backend('ibmq_16_melbourne')
        shots = 2**3 # On choisit un nombre relativement faible par rapport aux capacités du processeur 
                     # afin de ne pas avoir de temps d'exécution trop long
        job = execute(qc, backend,shots = shots)
        counts = job.result().get_counts()
        
        #----------- On récupère les probabilités des valeurs de la fonction f(x)
        print(counts)

        r = "111" # Initialisation du rayon sur la valeur la plus grande possible
        rnum = 0  # Rayon en int

        #----------- On cherche la clé la plus petite qui ne correspond pas à du bruit -----------
        for key,val in counts.items():
            if val > shots/(2**3): # Seuillage par rapport shots/8 (car on a 8 valeurs possibles pour r)
                if key < r and key != "000": # Recherche du plus petit pic dans counts
                    r = key

        #----------- On passe le rayon du cycle en entier -----------
        for i in range(3):
            rnum += (int(r[i]))*(2**(len(r)-i-1))

        #----------- On affiche les résultats si la valeur de r nous donne un résultat convenable ---------
        if (rnum % 2 != 1) or ((a ** rnum) % N == -1):
            print("Les diviseurs sont:",pgcd((a**rnum) +1,N),pgcd((a**rnum) -1,N))
        else:
            print("Erreur sur le calcul de r")
        stop = timeit.default_timer()
        print(start-stop)