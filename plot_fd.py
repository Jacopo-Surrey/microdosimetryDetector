import numpy as np
import matplotlib.pyplot as plt
import glob

def Normalize(f, x) :
	dx = x[1:] - x[:-1]
	integral = np.sum(f[:-1] * dx)
	f_norm = f[:-1] / integral
	return f_norm

# in cvs each thread outputs to a different file
# therefore I need to cycle through each thread and append them all together

energy = np.empty([0,1])
length = np.empty([0,1])

ntuple_name = "output_nt_2_deposited_energy"
output_name = ntuple_name + "_t" + "*" + ".csv"
output_list = glob.glob(output_name)

for this_file in output_list :	
	energy_thread, length_thread = np.loadtxt(this_file, delimiter=',', unpack=True)
	
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

bins_per_dec = 60	# change this depending on how many points you have

n_bins = n_decades * bins_per_dec

y_bins = np.zeros(n_bins)

y_bins[0] = 10**exp_start

for i in range(1, n_bins) :
	y_bins[i] = y_bins[i-1] * 10**( 1 / bins_per_dec )


# create the histogram

f = np.histogram( y, bins=y_bins ) [0]	# for now f is a number of counts...

tot_counts = np.sum(f)

# ... so now I turn f into a density
bin_width = y_bins[1:] - y_bins[:-1]
f = f / bin_width

f = np.append(f, 0.)	# give f and y_bins arrays the same size

# normalize the spectra
f = Normalize(f, y_bins)
d = Normalize(y_bins*f, y_bins)


# save to file

output_file = "analysed_spectra.csv"
header = "y[keV/um], f(y)[um/keV], d(y)[um/keV]"

np.savetxt( output_file, np.c_[ y_bins, f, d ], header=header, delimiter=',' )


# plot

fig, ax1 = plt.subplots()

color = 'tab:blue'

ax1.semilogx( y_bins, y_bins * f, linewidth=0.5, color=color  )

ax1.set_xlabel(r'$y \,\, [keV / \mu m]$')
ax1.set_ylabel(r'$y \cdot f(y) $', color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:red'

ax2.semilogx( y_bins, y_bins * d, linewidth=0.5, color=color )

ax2.set_ylabel(r'$y \cdot d(y) $', color=color)
ax2.tick_params(axis='y', labelcolor=color)

title = str(tot_counts) + " counts, " + str(bins_per_dec) + " bins per decade"
fig.suptitle(title)

fig.tight_layout()

plt.subplots_adjust(top=0.92)
plt.show()
