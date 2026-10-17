#!/bin/bash
# Complément de l'étude : taille d'échantillon représentative (tirage uniforme + Lloyd).
#
#  A. Effet de taille : tissus de L = 10 à 160 (80 à 32 600 cellules), k = 100 itérations, nombreuses graines ;
#     moyenne et dispersion entre graines des indicateurs en fonction du nombre de cellules.
#  B. Fenêtres : dans de grands tissus (L = 160, 3 graines, k = 0, 10, 100, 1000 ; L = 320, 1 graine, k = 100),
#     indicateurs calculés dans des fenêtres carrées de côté w ; leur dispersion d'une fenêtre à l'autre donne la
#     taille à partir de laquelle une portion de tissu est représentative, et S6(w) la taille des domaines
#     d'orientation.
#
# Usage : ./run_taille.sh     (environ 6 minutes), puis ./run.sh ou latexmk pour le rapport
set -e
export LC_ALL=C
cd "$(dirname "$0")"
ROOT=../..
C=$ROOT/cellyod
SRC=$ROOT/src
make -s -C $ROOT cellyod
c++ -std=c++17 -O2 -Wall -Wextra -I$SRC -o mesure mesure.cpp $SRC/params.cpp $SRC/expr.cpp $SRC/tissue.cpp $SRC/io.cpp
mkdir -p work/taille data

# moyenne et écart type, par valeur de la première colonne, de toutes les colonnes suivantes (en-tête conservé)
aggregate() {
  awk 'NR == 1 { for (i = 2; i <= NF; i++) name[i] = $i; nf = NF; next }
       { g = $1; if (!(g in cnt)) order[++no] = g; cnt[g]++
         for (i = 2; i <= nf; i++) { s[g, i] += $i; s2[g, i] += $i * $i } }
       END { printf "%s nb", key; for (i = 2; i <= nf; i++) printf " %s %s_sd", name[i], name[i]; print ""
             for (o = 1; o <= no; o++) { g = order[o]; printf "%s %d", g, cnt[g]
               for (i = 2; i <= nf; i++) { m = s[g, i] / cnt[g]; v = s2[g, i] / cnt[g] - m * m
                 printf " %.6g %.6g", m, (v > 0 ? sqrt(v) : 0) }
               print "" } }' key="$1"
}

# --- A. Effet de taille (k = 100)
K=100
: > work/taille/A.txt
for spec in 10:40 20:20 40:10 80:6 160:3; do
  L=${spec%%:*}; NS=${spec##*:}
  for s in $(seq 1 $NS); do
    out=work/taille/A_L${L}_s$s.txt
    [ -f $out ] || $C Lx=$L Ly=$L lattice=uniform seed=$s lloydIterations=$K output=$out svg=none > /dev/null
    ./mesure $out | awk -v L=$L 'NR > 1 { print L, $0 }' >> work/taille/A.txt
  done
  echo "A : L = $L ($NS graines)"
done
{ echo "L $(./mesure work/taille/A_L10_s1.txt | head -1)"; cat work/taille/A.txt; } | aggregate L > data/taille_A.dat

# --- B. Fenêtres
WINDOWS="2 3 4 6 8 12 16 24 32 48 64"
for s in 1 2 3; do
  prev=""
  for k in 0 10 100 1000; do
    out=work/taille/B_L160_s${s}_$k.txt
    if [ ! -f $out ]; then
      if [ -z "$prev" ]; then
        $C Lx=160 Ly=160 lattice=uniform seed=$s lloydIterations=$k output=$out svg=none > /dev/null
      else
        $C -load $prev lloydIterations=$k output=$out svg=none > /dev/null
      fi
    fi
    prev=$out
  done
  echo "B : L = 160, graine $s"
done
prev=""
for k in 10 100; do
  out=work/taille/B_L320_s1_$k.txt
  if [ ! -f $out ]; then
    if [ -z "$prev" ]; then
      $C Lx=320 Ly=320 lattice=uniform seed=1 lloydIterations=$k output=$out svg=none > /dev/null
    else
      $C -load $prev lloydIterations=$k output=$out svg=none > /dev/null
    fi
  fi
  prev=$out
done
echo "B : L = 320"
# indicateurs globaux des grands tissus (k = 10, 100, 1000)
for k in 10 100 1000; do
  ./mesure work/taille/B_L160_s*_$k.txt | awk 'NR > 1' > work/taille/B_global_$k.txt
done
for k in 0 10 100 1000; do
  : > work/taille/B_win_$k.txt
  files="work/taille/B_L160_s1_$k.txt work/taille/B_L160_s2_$k.txt work/taille/B_L160_s3_$k.txt"
  [ $k = 100 ] && files="$files work/taille/B_L320_s1_$k.txt"
  for w in $WINDOWS 96 128; do
    for f in $files; do
      ./mesure -windows $w $f 2> /dev/null | awk 'NR > 1' >> work/taille/B_win_$k.txt || true
    done
  done
  { echo "w n cvA sigmaN p6 Amoy S6"; cat work/taille/B_win_$k.txt; } | aggregate w > data/taille_B_$k.dat
done
# tableau : nombre de cellules d'une fenêtre pour que la dispersion relative d'un indicateur d'une fenêtre à
# l'autre (écart type / moyenne) descende sous 5 % et 2 % (interpolation log-log entre fenêtres)
{
  echo '\begin{tabular}{@{}lrrrrrr@{}}'
  echo '\toprule'
  echo ' & \multicolumn{2}{c}{$k=10$} & \multicolumn{2}{c}{$k=100$} & \multicolumn{2}{c}{$k=1000$} \\'
  echo 'Indicateur & 5\,\% & 2\,\% & 5\,\% & 2\,\% & 5\,\% & 2\,\% \\'
  echo '\midrule'
  for spec in "Amoy:aire moyenne (densité)" "p6:\$p_6\$" "sigmaN:\$\\sigma_n\$" "cvA:\$\\mathrm{cv}_A\$"; do
    col=${spec%%:*}; label=${spec#*:}
    printf '%s' "$label"
    for k in 10 100 1000; do
      for thr in 0.05 0.02; do
        awk -v col=$col -v thr=$thr 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
          { n = $(c["n"]); cv = $(c[col "_sd"]) / $(c[col]); if ($(c["nb"]) < 6) next
            if (cv <= thr && !done) { if (pn > 0) { t = (log(thr) - log(pcv)) / (log(cv) - log(pcv)); r = exp(log(pn) + t * (log(n) - log(pn))) } else r = n
              printf " & %d", r; done = 1 }
            pn = n; pcv = cv }
          END { if (!done) printf " & $>$%d", pn }' data/taille_B_$k.dat
      done
    done
    echo ' \\'
  done
  echo '\bottomrule'
  echo '\end{tabular}'
} > data/taille_tableau.tex
echo "données : data/taille_A.dat, data/taille_B_*.dat, data/taille_tableau.tex"
