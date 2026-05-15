import numpy as np 
import matplotlib.pyplot as plt 

def p(mu):
    return ((3/2)*mu - (1/2)*(mu**3))
 
def f_k(mu,k):
    resultat = mu 
    for i in range(k):
        resultat=p(resultat)
    return resultat

def s_k(mu,k=3):
    return (1/2)*(1-f_k(mu,k))

#création des valeurs de mu
mu_values = np.linspace(-1,1,500)

plt.figure(figsize=(8, 6))

#on teste pour différentes vaelurs de k ce que donne le polynôme s_k
for k in range(1,6):
    s_values = s_k(mu_values,k)
    plt.plot(mu_values,s_values,label=f'k={k}')
#on affiche la figure     
plt.title('Profils de coupure de Becke $s_k(\mu)$', fontsize=14)
plt.xlabel('$\mu$ (Coordonnée elliptique)', fontsize=12)
plt.ylabel('$s_k(\mu)$ (Poids spatial)', fontsize=12)

plt.xlim(-1, 1)
plt.ylim(0.0, 1.0) 

plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.show()