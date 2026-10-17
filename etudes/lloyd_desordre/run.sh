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
  inkscape work/vue_$k.svg --export-type=png --export-filename=figures/vue_$k.png -w 1280 > /dev/null 2>&1
done

# --- 3 bis. Cartes d'orientation : hexagones colorés selon l'orientation locale du réseau (arg psi6 / 6, teinte
#            cyclique sur 60°), défauts (n != 6) en gris foncé ; tissus 40 x 40 de la graine 1
orient_svg() { # fichier de valeurs par cellule (mesure -cells), fichier de polygones, côté L
  awk -v s=20 -v L=$3 '
    function hex(v) { return sprintf("%02x", int(255 * v + 0.5)) }
    function hsv(h, sa, va,   i, f, p, q, t, r, g, b) {
      h = (h % 360) / 60; i = int(h); f = h - i; p = va * (1 - sa); q = va * (1 - sa * f); t = va * (1 - sa * (1 - f))
      if (i == 0) { r = va; g = t; b = p } else if (i == 1) { r = q; g = va; b = p } else if (i == 2) { r = p; g = va; b = t }
      else if (i == 3) { r = p; g = q; b = va } else if (i == 4) { r = t; g = p; b = va } else { r = va; g = p; b = q }
      return "#" hex(r) hex(g) hex(b) }
    FNR == NR { if (FNR > 1) { n[$1] = $2; ph[$1] = $4 }; next }
    FNR == 1 { printf "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%d\" height=\"%d\">\n", L * s, L * s }
    $1 ~ /^[0-9]+$/ && NF > 3 { id = $1; col = (n[id] == 6 ? hsv(6 * ph[id], 0.7, 0.95) : "#2b2b2b")
      printf "<polygon fill=\"%s\" stroke=\"#ffffff\" stroke-width=\"0.6\" points=\"", col
      for (k = 0; k < $3; k++) printf "%.1f,%.1f ", $(4 + 2 * k) * s, (L - $(5 + 2 * k)) * s
      print "\"/>" }
    END { printf "<rect width=\"%d\" height=\"%d\" fill=\"none\" stroke=\"#000\" stroke-width=\"3\"/>\n</svg>\n", L * s, L * s }' "$1" "$2"
}
for f in uniform_s1_0010 uniform_s1_0100 uniform_s1_1000 hex_s1_1000; do
  $C -load work/$f.txt output=/dev/null svg=none polygons=work/poly_$f.txt > /dev/null
  ./mesure -cells work/$f.txt > work/cells_$f.txt
  orient_svg work/cells_$f.txt work/poly_$f.txt $L > work/orient_$f.svg
  inkscape work/orient_$f.svg --export-type=png --export-filename=figures/orient_$f.png -w 1600 > /dev/null 2>&1
done
# palette de la teinte (légende) : bande de 0 à 60°
awk 'BEGIN { print "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"360\" height=\"20\">"
  for (h = 0; h < 360; h++) { i = int(h / 60); f = h / 60 - i; p = 0.95 * 0.3; q = 0.95 * (1 - 0.7 * f); t = 0.95 * (1 - 0.7 * (1 - f)); v = 0.95
    if (i == 0) { r = v; g = t; b = p } else if (i == 1) { r = q; g = v; b = p } else if (i == 2) { r = p; g = v; b = t }
    else if (i == 3) { r = p; g = q; b = v } else if (i == 4) { r = t; g = p; b = v } else { r = v; g = p; b = q }
    printf "<rect x=\"%d\" width=\"1.2\" height=\"20\" fill=\"rgb(%d,%d,%d)\"/>\n", h, 255 * r, 255 * g, 255 * b }
  print "</svg>" }' > work/palette.svg
inkscape work/palette.svg --export-type=png --export-filename=figures/palette.png -w 360 > /dev/null 2>&1

# --- 4. Rapport (la section « taille représentative » utilise les données de run_taille.sh)
[ -f data/taille_tableau.tex ] || ./run_taille.sh
latexmk -lualatex -interaction=nonstopmode -halt-on-error rapport.tex > work/latex.log 2>&1 \
  || { tail -30 work/latex.log; exit 1; }
echo "rapport : rapport.pdf"
