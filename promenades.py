"""
Promenades aléatoires de mouches entre n pièces disposées en cercle.

Chaque pièce communique avec ses deux voisines par des portes de largeurs
différentes ; la probabilité de passer par une porte est proportionnelle à
sa largeur. On construit la matrice de transition de cette chaîne de Markov,
on cherche sa distribution stationnaire par la méthode de la puissance, et
on compare deux départs (mouches réparties partout, ou toutes dans la
pièce 1) pour n = 5 (impair) et n = 6 (pair).

Projet L3 Physique, CY Cergy Paris Université (octobre 2024).

Utilisation :
    python promenades.py            # résultats + figure
    python promenades.py --sauver   # enregistre la figure dans figures/
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt

SAUVER = "--sauver" in sys.argv

# Affichage à 4 chiffres après la virgule
np.set_printoptions(precision=4, suppress=True)

# Largeurs des portes (en cm) : L[i] = porte entre la pièce i et la pièce i+1
parametres = {
    5: np.array([125, 100, 75, 75, 75]),
    6: np.array([125, 100, 75, 75, 75, 75]),
}

N_TRANSITIONS = 1000


def matrice_transition(n, L):
    """P[i, j] = probabilité d'aller de la pièce i à la pièce j."""
    P = np.zeros((n, n))
    for i in range(n):
        P[i, (i + 1) % n] = L[i]                # porte vers la pièce suivante
        P[i, (i - 1 + n) % n] = L[(i - 1 + n) % n]  # porte vers la pièce précédente
    return P / P.sum(axis=1, keepdims=True)     # chaque ligne somme à 1


def methode_puissance(P, pi0, tol=1e-6, max_iter=10000):
    """
    Itère pi <- pi . P jusqu'à convergence.
    Renvoie la distribution obtenue, le nombre d'itérations et un booléen
    indiquant si l'itération a convergé.
    """
    pi = pi0.copy()
    for iteration in range(1, max_iter + 1):
        pi_nouveau = pi @ P
        if np.linalg.norm(pi_nouveau - pi, 1) < tol:
            return pi_nouveau, iteration, True
        pi = pi_nouveau
    return pi, max_iter, False


# ---------------------------------------------------------------------------
# Calculs
# ---------------------------------------------------------------------------

resultats = {}
for n, L in parametres.items():
    P = matrice_transition(n, L)

    # Deux états initiaux : mouches réparties uniformément, ou toutes en pièce 1
    init_eq = np.ones(n) / n
    init_un = np.zeros(n)
    init_un[0] = 1

    pi_eq, it_eq, cv_eq = methode_puissance(P, init_eq)
    pi_un, it_un, cv_un = methode_puissance(P, init_un)

    # Distribution après N_TRANSITIONS pas : p(t) = p(0) . P^t
    Pt = np.linalg.matrix_power(P, N_TRANSITIONS)
    resultats[n] = {
        "P": P,
        "pi_eq": (pi_eq, it_eq, cv_eq),
        "pi_un": (pi_un, it_un, cv_un),
        "dist_eq": init_eq @ Pt,
        "dist_un": init_un @ Pt,
        "dist_un_plus1": init_un @ Pt @ P,
    }


def etat(convergence, iterations):
    return f"convergé en {iterations} itérations" if convergence else "NON convergé : oscillation"


for n, res in resultats.items():
    print(f"\n===== {n} pièces =====")
    print("Matrice de transition P :")
    print(res["P"])
    pi, it, cv = res["pi_eq"]
    print(f"\nMéthode de la puissance, départ équiprobable ({etat(cv, it)}) :")
    print(pi)
    pi, it, cv = res["pi_un"]
    print(f"Méthode de la puissance, départ pièce 1 ({etat(cv, it)}) :")
    print(pi)
    print(f"\nAprès {N_TRANSITIONS} transitions, départ équiprobable : {res['dist_eq']}")
    print(f"Après {N_TRANSITIONS} transitions, départ pièce 1      : {res['dist_un']}")
    print(f"Après {N_TRANSITIONS + 1} transitions, départ pièce 1      : {res['dist_un_plus1']}")

# ---------------------------------------------------------------------------
# Figure : évolution de la probabilité d'être dans chaque pièce, départ pièce 1
# ---------------------------------------------------------------------------

T = 40
fig, axs = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for ax, (n, L) in zip(axs, parametres.items()):
    P = resultats[n]["P"]
    p = np.zeros(n)
    p[0] = 1
    historique = [p]
    for _ in range(T):
        p = p @ P
        historique.append(p)
    historique = np.array(historique)
    for piece in range(n):
        ax.plot(historique[:, piece], marker=".", label=f"Pièce {piece + 1}")
    parite = "impair : convergence" if n % 2 else "pair : oscillation sans fin"
    ax.set_title(f"n = {n} ({parite})")
    ax.set_xlabel("Nombre de transitions t")
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=8, loc="upper right")
axs[0].set_ylabel("Probabilité (départ : toutes les mouches en pièce 1)")
plt.tight_layout()

if SAUVER:
    os.makedirs("figures", exist_ok=True)
    plt.savefig("figures/promenades.png", dpi=90)
    print("\nFigure enregistrée : figures/promenades.png")
else:
    plt.show()
