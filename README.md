# Atomes

On génère par gaussian un fichier texte contenant les informations nécessaires au calcul de la densité électronique de la molécule d'urée, ainsi que la position des centres des atomes 

Fichier de référence pour la strucutre électronique de la molécule d'urée et la position des centres des atomes : urea.wfn

## Data

Fonctions permettant l'extraction des paramètres de la molécule d'urée depuis le fichier wfn

## Grille

Paramètres choisis : $p = 20$, lebedev à 53.

urea_grid_multiwfn.txt : Fichier texte contenant les valeurs $x, y, z$ pour la grille

## Calcul des poids de Hirshfeld

Les scripts concernant cette partie se trouvent dans le dossier /Partitionnement.

Les poids de Hirshfeld d'un atome dans une molécule, en fonction de la position, est donné par : 

$\omega_j(r) = \frac{\rho_j(r)}{\rho^M(r)}$

Avec : 

- $\rho_j(r)$ la densité électronique sphérique d"un atome isolé
- $\rho^M(r)$  la densité électronique de la promolécule

### Calcul des densités sphériques des atomes isolés et de la promolécule d'urée

Pour cela, on utilise Multiwfn, qui nous renvoie un fichier texte contenant, pour chaque point de la grille, la densité électronique de l'atome, ou de la promolécule.

Ces densités sont présentes dans le dossier Partitionnement/Weights_Urea_Grid

Les densités présentes dans le dossier Partitionnement/Weights_Urea_BI sont calculées avec la grille de Multiwfn et ne sont donc pas utilisés

### Calcul des poids : division par la densité électronique totale de la molécule

On calcule le terme $\rho^M(r)$, qui conerne la densité électronique de la molécule entière. 

- weights_arrays : script contenant les np-arrays de densités électroniques et des poids de Hirshfeld
- export_weights : script qui exporte les poids de H pour chaque atome dans un fichier texte
- calcul_weights : script de fonction calculant les poids de H
- plot_weights : script permettant  de visualiser les poids en 3D, pour chaque atome

## Calcul des facteurs de diffusion


### Faire les transformées de Fourier

## Calculs des facteurs de structure
