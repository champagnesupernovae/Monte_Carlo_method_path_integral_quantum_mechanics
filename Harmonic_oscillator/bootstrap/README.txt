Con il bootstrap prendo per ogni N la catena della quantità fisica desiderata, chiamata X.
Fissata la lunghezza di correlazione, la ricampiono più volte e calcolo media delle medie e la std.
Salvo su file (bootstrap_mean_sigma_len_file.txt) la media delle medie, la std e la lunghezza di correlazione.
Ripeto per varie lunghezze di correlazione.

Una volta fatto questo, uso plot_sigma_vs_corr_len.py per capire qual è per X la lunghezza di correlazione giusta da usare (quella al plateau)
Salvo quindi il valore di X corrispondente su un qualche file.