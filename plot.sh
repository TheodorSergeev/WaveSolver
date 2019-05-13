GPI=plot.gpi

> $GPI

echo "reset
set encoding utf8

set terminal pngcairo size 1485,1050 enhanced dashed font \"Verdana,20\"

set ylabel \"u\" offset 0,0,0
set xlabel \"x\" offset 0,0,0

set style line 10 lw 2.5 lc rgb \"#228B22\" pt 7 ps 1.4 # forestgreen

set yrange [-0.006:0.006]
set xrange [0:0.75]

" >> $GPI

#echo "plot \\" >> $GPI
NUM=0

echo ""
echo "Generating gnuplot script..."

for FILE in ./data/*.txt; do
  filename=$(basename -- "$FILE")
  time="${filename%.*}"

  graph_num=$(printf "%06d" "$NUM")

  echo "" >> $GPI
  echo "set title \"t = $time\"" >> $GPI
  echo "set output \"./graph/graph_$graph_num.png\"" >> $GPI
  echo "plot \"$FILE\" using 1:2 smooth unique notitle  ls 10" >>$GPI
  ((NUM++))
done
#  echo "plot \"$FILE\" using 1:2 smooth unique notitle  ls 10, \"$FILE\" using 1:2 notitle ls 10" >>$GPI

echo ""
echo "Running gnuplot script..."
gnuplot $GPI

echo ""
echo "Creating animation..."
ffmpeg -r 15 -i ./graph/graph_%06d.png -c:v libx264 movie.mp4 -v 0

echo ""
echo "Finished."
