#! /bin/bash -e

[ -d wiki ] || hg clone https://code.google.com/p/gnuplotutils.wiki/ wiki
hg -R wiki pull
hg -R wiki update
make wiki
for i in *.wiki; do
  cp -v $i ./wiki
done
