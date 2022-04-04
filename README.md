# quaryonic_eos
n-u-d version of quarkyonic matter EOS based on https://arxiv.org/abs/2004.08293

The code runs on Python 3. Download it here: https://www.anaconda.com/products/individual#linux. Parallel library 'joblib' is used to speed up program. It is inincluded in the folder 'joblib' here already. Check latest version of 'joblib' here https://github.com/joblib/joblib.

It should work by typing in terminal:

>python3 main_eos_table.py

A folder with tables of EOS will be generated.


Source files:

config.py: all parameters to be modified, parameter range and number of points.
main_eos_table.py: calculate EOS and generate tables.

eos_quarkyonic.py: define n-u-d version of quarkyonic matter EOS, and its interpolation. Class ‘EOS_Quarkyonic_Potential’ contains all the smooth EOS interpolation e.g. p(n), n(p), cs2(p).

eos_potential.py: define PNM potential form and fit to energy Sv+BE and slope L.

unitconvert: define all the physical constant and unit conversion.

sly.dat, togashi.dat: txt file of EOS data (used as crust EOS). ‘sly.dat’ comes from http://xtreme.as.arizona.edu/NeutronStars/ 

joblib: a library that used for Parallel processing, take advantage of multiple available cpu.


Output file:

EOS_table_regular_kFn: directory contains EOS table with regular neutron fermi momentum.

EOS_table_regular_nB: directory contains EOS table with regular baryon number density.

crust_eos.txt:  It comes from the crustal part of existed EOS table(Sly4 or Togashi). Four columns correspond to baryon number density nB (fm-3), energy density epsilon (MeV fm-3), pressure p (MeV fm-3) and sound speed square cs2/c2.

k0_kt_kmax.txt: neutron upper Fermi momentums (MeV) at: core-crust transition, quark drip transition, and 12*saturation density. Used to generate regular grid in neutron upper Fermi momentums. N1 points between k0 and kt, N2 points between kt and kmax.

core_nB_grid.txt: number density grid used to generate core EOS table. (fm-3)

core_neutron_fermi_momentum.txt, core_baryon_number_density.txt, core_energy_density.txt, core_pressure.txt, core_sound_speed_square.txt: These files contain 2D array. Column number corresponds to baryon number density grid. Row number stands for different EOS parameters. Three parameters are varied in the default table, symmetry energy slope L(I=0..N_i), shell thickness parameter Lambda (j=0..N_j) and quark drip density nt (k=0..N_k). Row index is i*N_j*N_k+j*N_k+k. In total, there should be N_i* N_j* N_k rows.
