# cellyod — génération de tissus cellulaires

![cellyod](celloyd.gif)

cellyod pave un domaine par un diagramme de Voronoï, régularisé par des itérations de Lloyd. Une fraction des
germes peut être des **bundles** (faisceaux vasculaires) qui attirent les autres cellules au sens d'une distance
de Mahalanobis : on obtient alors des tissus végétaux, du pavage hexagonal régulier aux tissus très hétérogènes
organisés en colonnes (modèle de Trivino *et al.*, *A continuous description of plant cellular tissues from
variational Voronoi morphogenesis*).

C++17, sans dépendance pour le cœur et la ligne de commande. L'interface graphique dépend en plus de raylib.

```bash
make                         # -> cellyod
cd examples
../cellyod tissu.txt         # -> tissu_cells.txt + tissu.svg
../cellyod tissu.txt alpha=0.1 "anisotropy=1000 90" svg=none
../cellyod -load tissu_cells.txt lloydIterations=100   # reprend les itérations là où elles s'étaient arrêtées
../cellyod -h                # liste des paramètres
```

Toutes les grandeurs sont **sans dimension** : l'unité de longueur est le diamètre équivalent visé des cellules.
Une cellule a donc une aire visée de π/4, et un domaine `Lx × Ly` contient environ `Lx·Ly/(π/4)` cellules.

## Algorithme

1. **Germes** (`lattice`) :
   - `hex` : réseau hexagonal de même densité, germes déplacés d'au plus `disorder` × le pas ;
   - `square` : grille carrée, même perturbation (avec `disorder 0`, grille exacte) ;
   - `random` (ou `poisson`) : tirage avec une distance minimale de 0.7 × l'espacement moyen ;
   - `uniform` : tirage uniforme sans distance minimale.
2. **Bundles** : une fraction `bundleFraction` des germes est tirée au hasard, plus le germe le plus proche de
   chaque point `bundle x y`.
3. **Itérations de Lloyd modifiées** (`lloydIterations`) : chaque germe ordinaire va en
   p' = (c + α b)/(1 + α), avec c le centroïde de sa cellule et b le bundle le plus proche au sens de
   Mahalanobis. C'est une descente de E = E_CVL − α E_att. Avec `alpha 0`, on retrouve Lloyd.
   Les bundles suivent `bundleMotion` : `lloyd` (vers le centroïde de leur cellule), `fixed` (immobiles) ou
   `attract` (attirés comme les autres cellules).
4. **Voronoï** découpé au domaine, par découpe de demi-plans (cellules convexes), puis fusion des sommets
   communs : le pavage est décrit par des sommets partagés.

Contrairement à nodeGen, il n'y a aucune reconformation (arêtes courtes, nombre d'arêtes) : le pavage est le
Voronoï exact des germes.

**Paramètres de forme.** `alpha` règle surtout l'hétérogénéité des aires (les cellules proches d'un bundle
rétrécissent) ; l'anisotropie `K` de Σ règle l'organisation spatiale (`K = 1` : amas autour de chaque bundle ;
`K ≫ 1` : colonnes orientées à `anisotropyAngle` de l'axe x). Les itérations ne convergent pas : à α et K fixés,
les cellules continuent de se resserrer, et `lloydIterations` est un troisième paramètre de forme.

| Tissu | Réglage | Exemple |
|---|---|---|
| hexagonal pur | `alpha 0` | |
| référence (amas autour des bundles) | `lattice random`, `bundleFraction 0.05`, `alpha 0.05`, `anisotropy 1` | `examples/tissu.txt` |
| colonnes | idem avec `anisotropy 1000` | `examples/colonnes.txt` |
| coupe tomographique | `lattice square`, `disorder 0`, `lloydIterations 5`, `alpha 0.02`, `anisotropy 10000` | `examples/tomographie.txt` |

## Paramètres

Fichier « mot-clé valeur(s) », `#` pour les commentaires. Sur la ligne de commande, `mot-clé=valeur` remplace la
valeur du fichier.

| Mot-clé | Défaut | Description |
|---|---|---|
| `Lx`, `Ly` | 40, 40 | dimensions du domaine [0, Lx] × [0, Ly] |
| `lattice` | hex | état initial : `hex`, `square`, `random` (ou `poisson`), `uniform` |
| `disorder` | 0.3 | perturbation des germes de `hex` et `square`, en fraction du pas |
| `seed` | 1 | graine du générateur aléatoire |
| `lloydIterations` | 10 | nombre d'itérations |
| `bundleFraction` | 0 | fraction des germes tirés au hasard comme bundles |
| `bundle x y` | — | bundle imposé : germe le plus proche du point (répétable) |
| `alpha` | 0 | intensité de l'attraction vers le bundle le plus proche |
| `anisotropy K [angle]` | 1 | anisotropie K de Σ (et orientation des colonnes, en degrés) |
| `anisotropyAngle` | 90 | orientation des colonnes, en degrés depuis l'axe x (90 : Σ = diag(1, K)) |
| `bundleMotion` | lloyd | `lloyd`, `fixed` ou `attract` |
| `output` | cells.txt | fichier de cellules produit |
| `svg` | cells.svg | vue SVG (`none` : pas de SVG) |
| `morphology` | — | ajoute une ligne de paramètres et de métriques au fichier (balayages) |

Σ = K a aᵀ + b bᵀ, avec a la direction des colonnes et b sa perpendiculaire. Toute matrice définie positive
s'écrit ainsi à un facteur près, et seule la forme de Σ compte.

## Fichier de cellules

```
cellyod 1
param <mot-clé> <valeur>        tous les paramètres du tissu
domain xmin ymin xmax ymax
iteration <n>                   nombre d'itérations effectuées
vertices <N>                    puis N lignes : x y
cells <M>                       puis M lignes : type sx sy nv v0 v1 ... v(nv-1)
```

Les sommets sont partagés entre cellules. Ceux d'une cellule sont donnés dans le sens trigonométrique. L'identifiant
d'une cellule est son rang, à partir de 0. `type` vaut 0 pour une cellule ordinaire et 1 pour un bundle.
`(sx, sy)` est le germe de la cellule : avec les paramètres, il suffit à reconstruire le pavage, ce qui permet
de reprendre les itérations (`-load`).

## Métriques

Affichées en fin de calcul et écrites dans le fichier `morphology`. Elles sont calculées sur toutes les cellules,
bords compris :

- ⟨SI⟩ = ⟨4πA/P²⟩ : indice de forme moyen (compacité ; 0.907 pour l'hexagone régulier) ;
- σ_A/⟨A⟩ : coefficient de variation des aires (hétérogénéité) ;
- ⟨d_NN⟩ : distance moyenne au plus proche voisin, entre centroïdes, rapportée à l'espacement visé √(π/4) ;
- ⟨d_M⟩ : distance de Mahalanobis moyenne d'une cellule ordinaire au bundle le plus proche, rapportée à
  l'espacement visé ;
- aire et SI moyens des cellules ordinaires et des bundles ;
- nombre moyen de côtés et histogramme (cellules intérieures) ;
- `lastMove` : déplacement moyen des germes à la dernière itération, rapporté à l'espacement visé.

## Organisation du code

| Fichier | Rôle |
|---|---|
| `src/geom.hpp` | vecteurs, polygones (aire, centroïde) |
| `src/params.hpp/.cpp` | paramètres et leur registre (`paramSpecs()`) |
| `src/tissue.hpp/.cpp` | domaine, germes, Voronoï, itérations, sommets partagés, métriques |
| `src/io.hpp/.cpp` | fichier de cellules, SVG, fichier `morphology` |
| `src/cli.cpp` | ligne de commande |

Le registre `paramSpecs()` décrit chaque paramètre une seule fois : nom, groupe, aide, bornes conseillées, et s'il
faut régénérer les germes après une modification. La lecture des paramètres, le fichier de cellules, l'aide
`-h` et le panneau de l'interface en sont déduits : **ajouter un paramètre = un champ dans `Params` + une ligne
dans `paramSpecs()`**. De même, `Metrics::items()` liste les métriques affichées.

La forme du domaine est isolée dans la structure `Domain` (`src/tissue.hpp`), pour l'instant un rectangle.
