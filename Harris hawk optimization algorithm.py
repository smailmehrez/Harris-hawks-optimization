from Fonction import*
import numpy as np
import random
import matplotlib.pyplot as plt
import math
import warnings
warnings.filterwarnings('ignore')

def HHO(func,Solution,nb_solu,taille_solu,itération):
    T = itération
    taille_solution = taille_solu
    nb_solution = nb_solu
    Matrice_Solution =np.copy(Solution)
    #********************HHO prinssipale*******************************
    idxMin = np.argmin(Matrice_Solution[:,taille_solution])
    solution_best_in_itra=np.copy(Matrice_Solution[idxMin,:])
    vecteur_Fonc_objectiv_meillure=np.full(shape=(T+1), fill_value=Matrice_Solution[idxMin,taille_solution],dtype=float)
    _,UB,_,LB=func(([1,2]))
    #******Inialiser les parametre utilise dans algorithme HHO******
    t = 1
    Smin = LB
    Smax = UB
    while t <= T:
        E1 = 1 - (t/T)
        for i in range (nb_solution):
            E0 = random.uniform(-1,1)
            E= 2 * E0 * E1 #****calcule E avec  l'Equation 3
            # --------  /* Phase d'exploration */
        #******Mettre à jour la position 𝑆𝑖𝑡 en utilisant l'Equation 1
            if abs(E) >= 1 :
                q = random.uniform(0,1)
                if q >= 0.5 :
                    r1 = random.uniform(0,1)
                    r2 = random.uniform(0,1)
                    Srand = random.randint(0,nb_solution-1)
                    Matrice_Solution[i,:-1] = Matrice_Solution[Srand,:-1] - r1 * abs(Matrice_Solution[Srand,:-1] - 2 * r2 * Matrice_Solution[i,:-1])
                else:
                    r3 = random.uniform(0,1)
                    r4 = random.uniform(0,1)
                    Matrice_Solution[i,:-1] = (solution_best_in_itra[:-1] - Matrice_Solution[i,:-1].mean()) -r3 * (Smin + r4 *(Smax - Smin))
            else:
                r = random.uniform(0,1)
                if abs(E) >= 0.5 and r >= 0.5:
                    #******Mettre à jour la position 𝑆𝑖𝑡 en utilisant l'Equation 4
                    r5 = random.uniform(0,1)
                    Matrice_Solution[i,:-1] = (solution_best_in_itra[:-1] - Matrice_Solution[i,:-1])   - E * abs(2 * (1 - r5) * solution_best_in_itra[:-1] - Matrice_Solution[i,:-1])    
                else :
                    if abs(E) < 0.5 and r >= 0.5:
                        #******Mettre à jour la position 𝑆𝑖𝑡 en utilisant l'Equation 6
                        Matrice_Solution[i,:-1] = solution_best_in_itra[:-1] - E * abs(solution_best_in_itra[:-1] - Matrice_Solution[i,:-1])
                    else :
                        if abs(E) >= 0.5 and r < 0.5:
                            #******Mettre à jour la position 𝑆𝑖𝑡 en utilisant l'Equation 10
                            r6 = random.uniform(0,1)
                            Y = solution_best_in_itra[:-1] - E * abs(2 * (1 - r6) * solution_best_in_itra[:-1] - Matrice_Solution[i,:-1])
                            beta=2
                            sigma=(math.gamma(1+beta)*math.sin(math.pi*beta/2)/(math.gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(1/beta)
                            u = random.uniform(0,1)
                            v = random.uniform(0,1) 
                            LF = 0.01 * (u * sigma) / (abs(v)**(1/beta))
                            S = random.uniform(0,1)
                            Z = Y + S * LF
                            FY,_,_,_ = func(Y)
                            FZ,_,_,_ = func(Z)
                            if FY < Matrice_Solution[i,taille_solution]:
                                Matrice_Solution[i,:-1] = Y

                            else:
                                if FZ < Matrice_Solution[i,taille_solution]:
                                    Matrice_Solution[i,:-1] = Z
                        else :
                            #******Mettre à jour la position 𝑆𝑖𝑡 en utilisant l'Equation 11
                            if abs(E) < 0.5 and r < 0.5:
                                r7 = random.uniform(0,1)
                                Y = solution_best_in_itra[:-1] - E * abs(2 * (1 - r7) * solution_best_in_itra[:-1] - Matrice_Solution[i,:-1].mean())
                                beta=1.5
                                sigma=(math.gamma(1+beta)*math.sin(math.pi*beta/2)/(math.gamma((1+beta)/2)*beta*2**((beta-1)/2)))**(1/beta)
                                u = random.uniform(0,1)
                                v = random.uniform(0,1) 
                                LF = 0.01 * (u * sigma) / (abs(v)**(1/beta))
                                S = random.uniform(0,1)
                                Z = Y + S * LF
                                FY,_,_ ,_= func(Y)
                                FZ,_,_,_ = func(Z)
                                if FY < Matrice_Solution[i,taille_solution]:
                                    Matrice_Solution[i,:-1] = Y
                                else:
                                    if FZ < Matrice_Solution[i,taille_solution]:
                                        Matrice_Solution[i,:-1] = Z
        # *********************traiter les borne*****************************
            for j in range(taille_solution):
                if Matrice_Solution[i,j] < LB:
                    Matrice_Solution[i,j] = "{0:.3f}".format(random.uniform(LB,0))
                if Matrice_Solution[i,j] > UB:
                    Matrice_Solution[i,j] = "{0:.3f}".format(random.uniform(0,UB))
        # ***********************calcule la fonction objective***************     
            Matrice_Solution[i,taille_solution],_,_,_ = func(Matrice_Solution[i,:-1])
        #********************************************************************
        idxMin = np.argmin(Matrice_Solution[:,taille_solution])
        if (solution_best_in_itra[taille_solution]> Matrice_Solution[idxMin,taille_solution]):
            solution_best_in_itra=np.copy(Matrice_Solution[idxMin,:])
        vecteur_Fonc_objectiv_meillure[t] = solution_best_in_itra[taille_solution]
        t = t + 1
 
    return(vecteur_Fonc_objectiv_meillure)



T=100
n=10
fonction=Quartic
_,UB,taille_solution,LB=fonction(([1,2]))
nb_solution=n
# crer matrice du Solution
M1 = np.full(shape=(nb_solution, taille_solution+1), fill_value=0,dtype=np.float64)
# inialiser le Matrice du solution random
for j in range(nb_solution):
    for i in range(taille_solution):
        M1[j,i] = "{0:.3f}".format(random.uniform(LB, UB)) 
# calcule la fonction objective 
for j in range(nb_solution):
    M1[j,taille_solution],_,_,_=fonction(M1[j,:-1])
vecteur_best=HHO(fonction,M1,nb_solution,taille_solution,T)
ax1 =plt.axes()
ax1.plot(vecteur_best,label="HHO",c="#4B2991")
ax1.set_ylabel('Coût F')
ax1.set_xlabel("Nombre d'itération")
if fonction== Ackley or fonction == Michalewicz or fonction== Schwefel1 or fonction== Schwefel2_26 or fonction== Himmelblau:
    ax1.legend()
    ax1.set_title("Le Coût f en fonction de Nombre d'itération")
    plt.show()
else:
    ax1.set_yscale('log')
    ax1.legend()
    ax1.set_title("Le Coût f en fonction de Nombre d'itération")
    plt.show()