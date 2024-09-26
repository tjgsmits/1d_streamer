# 1d_streamer
This code contains an 1D fluid model to simulate a drift-diffusion equation,
$\frac{\partial n}{\partial t} = \frac{\partial}{\partial x}(\mu n \textbf{E} + D\frac{\partial n}{\partial x}) + S$
which is coupled with 
$E = - \frac{\partial \phi}{\partial x}$
$\frac{\partial^2 \phi}{\partial x^2}$ = \frac{e}{\epsilon_0} (n_i - n_e)$
which is solved using a finit volum method (FVM). 

The code can be complide with
--> $ make
and used with
--> $ ./streamer_1d

The script plot_profiles.py can be used to plot the E, phi, n_electron and n_ion as function of position at one time slice.
Multiple time slices can be plotted. The script can be used with
--> $ python plot_profiles.py data_output step1 step2 ... stepN --colors color1 color2 ... colorN

Parameters of the simulation and input can be changed in the file streamer.f90
