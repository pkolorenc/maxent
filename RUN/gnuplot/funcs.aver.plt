set style data l
Ed = real(system("awk 'NR==1 {print $3}' coup.0001"))
set arrow 1 from Ed,0 to Ed,1000 lt 0 nohead
plot [1:2*Ed] 'func.aver.03_15' u ($1+Ed):($2-$3):($2+$3) with filledcurves lt 1, '' u ($1+Ed):2 lt 1, 'func.aver.05_20' u ($1+Ed):2 lt 2
pause -1
