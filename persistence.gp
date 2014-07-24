#!/usr/bin/gnuplot -p
#set term pdfcairo enhanced size 8cm,4.8cm

set fit quiet
set xlabel '{/Symbol t} (bases)'
set ylabel 'C ({/Symbol t})'


f(x)=exp(-x/lp)*(a+(1-a)*cos(2*pi*x/lambda))
lp=100
lambda=10
a=0.5

#set xrange [:50]

fit f(x) 'persistence' via a,lp,lambda

#i#t xrange [:45]

plot 'persistence', f(x)

show variable l
