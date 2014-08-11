#! /bin/bash -e

gplot="gplot"
[ -z "$(type -p $gplot)" ] && gplot="../gplot"
gp2eps="gp2eps"
[ -z "$(type -p $gp2eps)" ] && gp2eps="../gp2eps"
gv="gv"
[ -z "$(type -p $gv)" ] && gv="kghostview"
[ -z "$(type -p $gv)" ] && gv="okular"
[ -z "$(type -p $gv)" ] && gv="echo gv"

pause () {
  echo "Press enter"
  read
}

echo "First we generate a dummy datafile"
pause
for ((i=1;i<=5;i++)); do
  echo $i $((2*$i)) $((3*$i))
done > datafile
echo "Done"

echo "Now datafile looks like this"
cat datafile
pause

echo "The simple plot can be done by 'gplot datafile'"
$gplot datafile
pause

echo "Any gnuplot plot options can be appended after the filename"
echo "'gplot datafile w lp t \"data\"' does (mind the \"!!)"
$gplot datafile w lp t \"data\"
pause

echo "or with eps output with 'gplot -o data.eps datafile'"
$gplot -o data.eps datafile
$gv data.eps
echo "Not nice, but fast. You are just a sec away from 100 eps plots ;-)"
pause

echo "useful options are -x,-y,-z,-t for x,y,z-label and title"
echo "and range [??:??] in front of the file"
echo "'gplot -o data2.eps -x x -y \"My data\" -t "Useful" [1:3] datafile' creates"
$gplot -o data2.eps -x x -y "My data" -t Useful [1:3] datafile w l
$gv data2.eps
echo "Nicer, but still creepy"
pause

echo "To see what gplot is acually doing add -p option"
$gplot -p -o data2.eps -x x -y "My data" -t Useful [1:3] datafile w l
pause

echo "Let us ignore the beauty issue for a moment and come to some nice gplot functions"
echo "gplot can read from stdin 'seq 1 10 | gplot'"
seq 1 10 | $gplot
pause

echo "Add a - if you want to add more options: 'seq 1 10 | gplot - w l'"
seq 1 10 | $gplot - w l
pause
echo "Close all gplot windows now"
pause

echo "gplot can automatically replot if file was modified 'gplot --replot datafile2'"
echo "Don't close the windows now !"
seq 1 10 > datafile2
$gplot --replot datafile2 w l & 
pid=$!
for ((i=20;i<50;i+=10)); do
  sleep 2
  echo "Changing datafile2"
  seq 1 $i >> datafile2
done
kill $pid
pause

echo "Generating nice plots for publications etc. you can use --gp2eps option"
echo "'gplot --gp2eps datafile w lp t \"Data\"' will generate a datafile.gp with content:"
$gplot --gp2eps datafile w lp t \"data\"
cat datafile.gp
pause

echo "Now we add some more stuff to the gp file to make the plot a lot nicer"
sed -i '2i set title "A title"' datafile.gp
sed -i '2i set xlabel "$X_\\\\theta$ [unit]"' datafile.gp
sed -i '2i set mxtics 5' datafile.gp
sed -i '2i set ylabel "$Y_i$ [unit]"' datafile.gp
sed -i '2i set mytics 5' datafile.gp
sed -i '2i set key bottom' datafile.gp
sed -i '2i set format x "%.1f"' datafile.gp
echo "Now datafile.gp looks like this:"
cat datafile.gp
echo "Running 'gp2eps datafile.gp' will generate datafile.eps"
pause
$gp2eps datafile.gp
$gv datafile.eps

echo "gp2eps has 4 scales (see the manpage of gp2eps)"
echo "Most important there is -s which is by default 0.5"
echo "different -s scales allows make multicolumn plots"
echo "with -s 0.5 include with \includegraphics[width=0.5\textwidth]{} (two columns)"
echo "with -s 1.0 include with \includegraphics[width=1\textwidth]{} (one column)"
echo "with -s 0.33 include with \includegraphics[width=0.33\textwidth]{} (3 columns)"
echo
echo "For more tricks see the manpage of gp2eps"
echo
echo "Have fun - feel free to add comments on gnuplotutils.googlecode.com"
