# Étude : désordre d'un pavage de Voronoï au fil des itérations de Lloyd

Lloyd pur (sans bundles), trois états initiaux (`uniform`, `random`, `hex` avec `disorder 0.5`), domaine 40 × 40,
5 graines, de 0 à 1000 itérations. Les indicateurs suivis portent sur :

- les tailles et les formes : dispersion des aires et des longueurs de paroi, énergie de Lloyd normalisée G ;
- la topologie : nombre de côtés, fraction d'hexagones, ordre orientationnel local ψ6 ;
- l'orientation des parois : ordre orientationnel global S6, anisotropie S2, dispersion des angles aux sommets ;
- la connectivité : loi d'Aboav-Weaire, appariement des pentagones et heptagones, taux de transitions T1.

Un complément mesure la taille d'échantillon représentative (tissus de 80 à 130 000 cellules, indicateurs par
fenêtres, taille des domaines d'orientation) et en déduit des ordres de grandeur pour une éprouvette de
fissuration pré-entaillée.

```bash
./run.sh          # environ 1 minute : simulations, mesures, figures, rapport.pdf
./run_taille.sh   # complément (environ 6 minutes) ; lancé par run.sh si ses données manquent
```

| Fichier | Rôle |
|---|---|
| `run.sh` | réglages de l'étude (en tête), simulations (`cellyod`, reprises `-load`), agrégation, figures, rapport |
| `run_taille.sh` | complément : effet de taille et analyse par fenêtres |
| `mesure.cpp` | indicateurs de désordre d'un fichier de cellules (utilise le cœur de cellyod) ; `-hist` : histogrammes des angles ; `-windows w` : indicateurs par fenêtre |
| `rapport.tex` | rapport (LuaLaTeX, pgfplots) ; figures et tableau lus dans `data/` et `figures/` |

Seules les sources sont versionnées ; `work/`, `data/`, `figures/`, `mesure` et `rapport.pdf` sont produits par
`run.sh`.
