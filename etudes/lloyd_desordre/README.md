# Étude : désordre d'un pavage de Voronoï au fil des itérations de Lloyd

Lloyd pur (sans bundles), trois états initiaux (`uniform`, `random`, `hex` avec `disorder 0.5`), domaine 40 × 40,
5 graines, de 0 à 1000 itérations. Neuf indicateurs sont suivis : dispersion des aires et des longueurs de
paroi, nombre de côtés, fraction d'hexagones, ordre orientationnel ψ6, énergie de Lloyd normalisée G, etc.

```bash
./run.sh      # environ 1 minute : simulations, mesures, figures, rapport.pdf
```

| Fichier | Rôle |
|---|---|
| `run.sh` | réglages de l'étude (en tête), simulations (`cellyod`, reprises `-load`), agrégation, figures, rapport |
| `mesure.cpp` | indicateurs de désordre d'un fichier de cellules (utilise le cœur de cellyod) |
| `rapport.tex` | rapport (LuaLaTeX, pgfplots) ; figures et tableau lus dans `data/` et `figures/` |

Seules les sources sont versionnées ; `work/`, `data/`, `figures/`, `mesure` et `rapport.pdf` sont produits par
`run.sh`.
