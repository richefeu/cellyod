#!/bin/bash
# Régénère les figures et les données du tutoriel (nécessite inkscape et ImageMagick).
# Usage : cd doc && ./make_figures.sh && latexmk -pdf tutorial_cellyod.tex
set -e
cd "$(dirname "$0")"
C=../cellyod
mkdir -p figures work
( cd .. && make -s cellyod cellyod-gui )

png() { inkscape "$1" --export-type=png --export-filename="$2" -w "$3" > /dev/null 2>&1; }

# vignette : nom, puis paires clé=valeur passées à cellyod
fig() {
  local name=$1; shift
  $C output=work/$name.txt svg=work/$name.svg "$@" > work/$name.log
  png work/$name.svg figures/$name.png ${W:-500}
}

# 1. Germes initiaux (sans itération)
for l in hex square random uniform; do
  fig lattice_$l Lx=14 Ly=14 lattice=$l lloydIterations=0
done

# 2. Régularisation de Lloyd seule
for n in 0 5 50; do
  fig lloyd_$n Lx=14 Ly=14 lattice=uniform lloydIterations=$n
done

# 3. Attraction vers les bundles : évolution avec les itérations
for n in 0 10 50 200; do
  fig iter_$n Lx=30 Ly=30 lattice=random bundleFraction=0.05 alpha=0.05 lloydIterations=$n svgSeeds=0
done

# 4. Espace des formes : alpha x K (50 itérations)
for a in 0 0.05 0.15; do
  for K in 1 100 1000; do
    fig shape_a${a}_K$K Lx=30 Ly=30 lattice=random bundleFraction=0.05 alpha=$a "anisotropy=$K" \
      lloydIterations=50 svgSeeds=0
  done
done
fig angle30 Lx=30 Ly=30 lattice=random bundleFraction=0.05 alpha=0.05 "anisotropy=1000 30" lloydIterations=50 svgSeeds=0
fig motion_fixed Lx=30 Ly=30 lattice=random bundleFraction=0.05 alpha=0.05 bundleMotion=fixed lloydIterations=50 svgSeeds=0
fig tomo Lx=30 Ly=30 lattice=square disorder=0 bundleFraction=0.05 alpha=0.02 anisotropy=10000 lloydIterations=5 svgSeeds=0

# 5. Formes de domaine
W=600 fig dom_disque Lx=40 Ly=40 "shape=r - 18" lattice=random bundleFraction=0.05 alpha=0.05 lloydIterations=50 svgSeeds=0
W=600 fig dom_anneau Lx=40 Ly=40 "shape=max(r - 18, 9 - r)" lattice=random bundleFraction=0.05 alpha=0.05 lloydIterations=50 svgSeeds=0
W=600 fig dom_ellipse Lx=40 Ly=40 "shape=((x-cx)/18)^2 + ((y-cy)/11)^2 - 1" "anisotropy=1000 0" \
  lattice=random bundleFraction=0.05 alpha=0.05 lloydIterations=50 svgSeeds=0
W=1200 fig rein Lx=46 Ly=46 "shape=r - 20 + 4.5*exp(-t^2/0.12)" lattice=random bundleFraction=0.05 alpha=0.05 lloydIterations=50
# zoom sur le bord échancré (à droite)
magick figures/rein.png -crop 300x300+890+450 +repage figures/rein_zoom.png
# même forme sans le poids beta (cellules du bord attirées comme les autres)
W=600 fig rein_beta0 Lx=46 Ly=46 "shape=r - 20 + 4.5*exp(-t^2/0.12)" beta=0 lattice=random bundleFraction=0.05 \
  alpha=0.05 lloydIterations=50 svgSeeds=0

# 6. Interface graphique (capture)
( cd work && rm -f cellyod-gui.ini && ../../cellyod-gui ../../examples/rein.txt --screenshot gui.png )
cp work/gui.png figures/gui.png
( cd work && ../../cellyod-gui ../../examples/tissu.txt --popup files --screenshot gui_files.png )
magick work/gui_files.png -crop 600x500+340+160 +repage figures/gui_files.png

# 7. Données : métriques en fonction des itérations et de alpha (pgfplots)
hdr() { sed -e '1s/^# //' "$1" > "$2"; }
for K in 1 1000; do
  rm -f work/iter_K$K.txt
  for n in 0 5 10 20 30 40 50 60 80 100 125 150 200; do
    $C Lx=30 Ly=30 lattice=random bundleFraction=0.05 alpha=0.05 anisotropy=$K lloydIterations=$n \
      output=/dev/null svg=none morphology=work/iter_K$K.txt > /dev/null
  done
  hdr work/iter_K$K.txt figures/iter_K$K.dat
  rm -f work/alpha_K$K.txt
  for a in 0 0.01 0.02 0.03 0.05 0.075 0.1 0.125 0.15 0.2; do
    $C Lx=30 Ly=30 lattice=random bundleFraction=0.05 alpha=$a anisotropy=$K lloydIterations=50 \
      output=/dev/null svg=none morphology=work/alpha_K$K.txt > /dev/null
  done
  hdr work/alpha_K$K.txt figures/alpha_K$K.dat
done

# 8. Extrait d'un fichier de cellules et sortie de la ligne de commande
$C Lx=6 Ly=6 lattice=hex lloydIterations=3 output=work/petit.txt svg=none > /dev/null
head -22 work/petit.txt > figures/petit_extrait.txt
( cd ../examples && ../cellyod tissu.txt output=../doc/work/tissu_cells.txt svg=none ) > figures/sortie_tissu.txt
$C -h > figures/aide.txt
echo "figures régénérées dans doc/figures"
