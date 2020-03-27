import numpy as np
import matplotlib.pyplot as plt

# in cvs each thread outputs to a different file
# therefore I need to cycle through each thread and append them all together

ntuple_name = "output_nt_102"

n_threads = 4

energy = np.empty([0,1])
length = np.empty([0,1])

for i in range(0, n_threads) :
	
	file_name = ntuple_name + "_t" + str(i) + ".csv"
	
	energy_thread, length_thread = np.loadtxt(file_name, delimiter=',', unpack=True)
	
	energy = np.append(energy, energy_thread)
	length = np.append(length, length_thread)

mean_path_length = np.average(length)
	
y = energy / mean_path_length


# define the log binning

minimum = np.amin(y)
maximum = np.amax(y)

exp_start = np.floor(np.log10(minimum))
exp_end = np.ceil(np.log10(maximum))

n_decades = int(exp_end - exp_start)

bins_per_dec = 300	# change this depending on how many points you have

n_bins = n_decades * bins_per_dec

y_bins = np.zeros(n_bins)

y_bins[0] = 10**exp_start

for i in range(1, n_bins) :
	y_bins[i] = y_bins[i-1] * 10**( 1 / bins_per_dec )


# create the histogram

f = np.histogram( y, bins=y_bins ) [0]	# now f is a number of counts...

tot_counts = str(np.sum(f))

for i in range(n_bins-1) :	# f has dimension n_bins-1
		
		# ... so I turn it into a density
		bin_width = y_bins[i+1] - y_bins[i]	
		f[i] = f[i] / bin_width
		
		# I also shift the abscissa to the middle of the bin
		y_bins[i] = y_bins[i] * 10**( 1 / bins_per_dec / 2 )

y_bins = y_bins[0:-1]	#remove the last point,
							# so both arrays have the same size


# normalize the spectra

def normalize(f, y) :
	integral = np.trapz(f, y)
	f_norm = f / integral
	return f_norm

f = normalize(f, y_bins)
d = normalize(y_bins*f, y_bins)


# plot

fig, ax1 = plt.subplots()

color = 'tab:blue'

ax1.semilogx( y_bins, y_bins * f, linewidth=0.5, label="???", color=color  )

ax1.set_xlabel(r'$y \,\, [keV]$')
ax1.set_ylabel(r'$y \cdot f(y) $', color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'

ax2.semilogx( y_bins, y_bins * d, linewidth=0.5, label= r'???', color=color )

ax2.set_ylabel(r'$y \cdot d(y) $', color=color)
ax2.tick_params(axis='y', labelcolor=color)

title = tot_counts + " counts, " + str(bins_per_dec) + " bins per decade"
fig.suptitle(title)

ax1.legend(fancybox=True, framealpha=1, loc="best")
fig.tight_layout()

plt.subplots_adjust(top=0.92)
plt.show()
