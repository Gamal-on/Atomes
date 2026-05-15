# Atomes

Fichier de référence pour la strucutre électronique de la molécule d'urée et la position des centres des atomes : urea.wfn

Étapes : 

## Grille

Paramètres choisis : $p = 20$, lebedev à 53.

## Calcul des poids de Hirshfeld

Weights_Urea_Grid : Dossiers contenant les densités électroniques des atomes isolés (sphéricalisées), sur notre grille, avec les  centres définis comme sur urea.wfn

$\omega_j(r) = \frac{\rho_j(r)}{\rho^M(r)}$

### Calcul des densités sphériques pour atomes isolés

On calcule le terme $\rho_j(r)$, qui concerne la densité des atomes de la molécule pris isolément.

### Calcul des poids : division par la densité électronique totale de la molécule

On calcule le terme $\rho^M(r)$, qui conerne la densité électronique de la molécule entière. 

## Calcul des facteurs de diffusion


### Faire les transformées de Fourier

## Calculs des facteurs de structure
