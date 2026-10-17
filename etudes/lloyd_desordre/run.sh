#!/bin/bash
# Étude : désordre d'un pavage de Voronoï au fil des itérations de Lloyd (sans bundles).
#
# Pour chaque état initial et chaque graine, cellyod itère de point de contrôle en point de contrôle (reprise
# avec -load : identique à un calcul d'une traite) ; mesure calcule les indicateurs de désordre de chaque fichier ;
# les moyennes et écarts types sur les graines sont écrits dans data/, les figures dans figures/, puis le rapport
# est compilé (latexmk, LuaLaTeX).
#
# Usage : ./run.sh            (environ 2 minutes)
set -e
export LC_ALL=C
cd "$(dirname "$0")"
ROOT=../..
C=$ROOT/cellyod
SRC=$ROOT/src

# --- Réglages de l'étude
L=40                                   # domaine L x L (environ 2000 cellules)
STATES="uniform random hex"            # états initiaux (hex : avec disorder 0.5)
SEEDS="1 2 3 4 5"
CHECKPOINTS="0 1 2 3 5 7 10 15 20 30 50 70 100 150 200 300 500 700 1000"

make -s -C $ROOT cellyod
c++ -std=c++17 -O2 -Wall -Wextra -I$SRC -o mesure mesure.cpp $SRC/params.cpp $SRC/expr.cpp $SRC/tissue.cpp $SRC/io.cpp
mkdir -p work data figures

# --- 1. Simulations
for st in $STATES; do
  extra=""
  [ "$st" = hex ] && extra="disorder=0.5"
  for s in $SEEDS; do
    prev=""
    for k in $CHECKPOINTS; do
      out=$(printf "work/%s_s%d_%04d.txt" $st $s $k)
      if [ -z "$prev" ]; then
        $C Lx=$L Ly=$L lattice=$st $extra seed=$s lloydIterations=$k output=$out svg=none > /dev/null
      else
        $C -load $prev lloydIterations=$k output=$out svg=none > /dev/null
      fi
      prev=$out
    done
  done
  echo "simulations $st : terminées"
done

# --- 2. Mesures, puis moyenne et écart type sur les graines (une ligne par point de contrôle)
for st in $STATES; do
  : > work/mesures_$st.txt
  for s in $SEEDS; do
    ./mesure work/${st}_s${s}_*.txt | awk -v s=$s 'NR == 1 { if (s == 1) print "seed", $0; next } { print s, $0 }' \
      >> work/mesures_$st.txt
  done
  awk 'NR == 1 { for (i = 3; i <= NF; i++) name[i] = $i; nf = NF; next }
       { it = $2; if (!(it in cnt)) order[++no] = it; cnt[it]++
         for (i = 3; i <= nf; i++) { s[it, i] += $i; s2[it, i] += $i * $i } }
       END { printf "iteration"; for (i = 3; i <= nf; i++) printf " %s %s_sd", name[i], name[i]; print ""
             for (o = 1; o <= no; o++) { it = order[o]; printf "%d", it
               for (i = 3; i <= nf; i++) { m = s[it, i] / cnt[it]; v = s2[it, i] / cnt[it] - m * m
                 printf " %.6g %.6g", m, (v > 0 ? sqrt(v) : 0) }
               print "" } }' work/mesures_$st.txt > data/$st.dat
done
# tableau du rapport : valeurs moyennes à quelques points de contrôle
{
  echo '\begin{tabular}{@{}lrrrrrrr@{}}'
  echo '\toprule'
  echo 'État & $k$ & $\mathrm{cv}_A$ & $\mathrm{cv}_L$ & $\sigma_n$ & $p_6$ & $\psi_6$ & $G$ \\'
  for st in $STATES; do
    echo '\midrule'
    awk -v st=$st 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
      $1 == 0 || $1 == 10 || $1 == 100 || $1 == 1000 {
        printf "%s & %d & %.3f & %.3f & %.3f & %.3f & %.3f & %.5f \\\\\n", ($1 == 0 ? "\\code{" st "}" : ""), $1,
          $(c["cvA"]), $(c["cvL"]), $(c["sigmaN"]), $(c["p6"]), $(c["psi6"]), $(c["G"]) }' data/$st.dat
  done
  echo '\midrule'
  echo 'hexagonal parfait & & 0 & 0 & 0 & 1 & 1 & 0.08019 \\'
  echo '\bottomrule'
  echo '\end{tabular}'
} > data/tableau.tex
# second tableau : orientation des parois et connectivité
{
  echo '\begin{tabular}{@{}lrrrrrrr@{}}'
  echo '\toprule'
  echo 'État & $k$ & $S_6$ & $S_2$ & $\sigma_V$ (°) & $a$ & $f_{57}$ & $\tau_{T1}$ \\'
  for st in $STATES; do
    echo '\midrule'
    awk -v st=$st 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
      $1 == 0 || $1 == 10 || $1 == 100 || $1 == 1000 {
        printf "%s & %d & %.3f & %.4f & %.1f & %.2f & %.3f & %s \\\\\n", ($1 == 0 ? "\\code{" st "}" : ""), $1,
          $(c["S6"]), $(c["S2"]), $(c["sigmaV"]), $(c["aboav"]), $(c["f57"]),
          ($1 == 0 ? "--" : sprintf("%.1e", $(c["t1"]))) }' data/$st.dat
  done
  echo '\midrule'
  echo 'hexagonal parfait & & 1 & 0 & 0 & -- & 0 & -- \\'
  echo '\bottomrule'
  echo '\end{tabular}'
} > data/tableau2.tex
# histogrammes des angles aux sommets et des orientations de paroi (graine 1)
for f in uniform_s1_0000 uniform_s1_0010 uniform_s1_1000 hex_s1_1000; do
  ./mesure -hist work/$f.txt > data/hist_$f.dat
done
# référence : pavage hexagonal parfait
$C Lx=$L Ly=$L lattice=hex disorder=0 lloydIterations=0 output=work/parfait.txt svg=none > /dev/null
./mesure work/parfait.txt > data/parfait.dat

# --- 3. Vues : petit domaine, cellules colorées par nombre de côtés (à partir de l'export polygones)
side_svg() { # fichier de polygones -> SVG (palette de cellyod-gui : 4 rouge, 5 orange, 6 vert, 7 bleu clair...)
  awk -v s=40 -v L=16 'BEGIN { split("8e44ad d64541 f0963c cde4be 5faae1 2d6eb4 193c6e", pal, " ")
      printf "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n", L * s, L * s }
    $1 ~ /^[0-9]+$/ && NF > 3 { n = $3; c = n < 3 ? 3 : (n > 9 ? 9 : n)
      printf "<polygon fill=\"#%s\" stroke=\"#1a1a1a\" stroke-width=\"1\" points=\"", pal[c - 2]
      for (k = 0; k < n; k++) printf "%.1f,%.1f ", $(4 + 2 * k) * s, (L - $(5 + 2 * k)) * s
      print "\"/>" }
    END { printf "<rect width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"#000\" stroke-width=\"3\"/>\n</svg>\n", L * s, L * s }' "$1"
}
prev=""
for k in 0 3 20 200; do
  out=$(printf "work/vue_%04d.txt" $k)
  if [ -z "$prev" ]; then
    $C Lx=16 Ly=16 lattice=uniform seed=2 lloydIterations=$k output=$out svg=none polygons=work/vue_poly_$k.txt > /dev/null
  else
    $C -load $prev lloydIterations=$k output=$out svg=none polygons=work/vue_poly_$k.txt > /dev/null
  fi
  prev=$out
  side_svg work/vue_poly_$k.txt > work/vue_$k.svg
  inkscape work/vue_$k.svg --export-type=png --export-filename=figures/vue_$k.png -w 640 > /dev/null 2>&1
done

# --- 4. Rapport (la section « taille représentative » utilise les données de run_taille.sh)
[ -f data/taille_tableau.tex ] || ./run_taille.sh
latexmk -lualatex -interaction=nonstopmode -halt-on-error rapport.tex > work/latex.log 2>&1 \
  || { tail -30 work/latex.log; exit 1; }
echo "rapport : rapport.pdf"
