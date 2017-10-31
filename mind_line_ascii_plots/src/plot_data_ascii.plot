# gnuplot script to send ascii plot to terminal
# designed for dummy terminals without access to host 
# or other graphics

# set the terminal
set terminal dumb

# plot CTE
set title 'CTE'
set style fill solid 0.3
set style data lines
plot '../data/cte_vals.txt' notitle

# plot detla
set title 'Delta (Radians)'
plot '../data/delta_vals.txt' notitle

# plot Velocity
set title 'Velocity'
plot '../data/v_vals.txt' notitle
