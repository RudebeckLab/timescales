# timescales
Analysis code for Zeisler et al. (2024): Consistent hierarchies of single-neuron timescales in mice, macaques and humans

Inside the *matlab_timescale_extraction* folder are the MATLAB scripts used to estimate timescales from arrays of spiketimes; AC_spiketimes.m uses the spike-count autocorrelation method from Murray et al., 2014, and ISI_spiketimes.m uses the ISI-based approach introduced by Fontanier et al., 2022.

Stepping through the scripts in the top-level folder starting with 0-6 will process and plot the data shown in the main manuscript figures. Within each script, comments are used to indicate where each figure is plotted.

Inside the *method_comparison* folder are the scripts used to process the data toward Figure S2. As before, stepping through scripts 0-3 will generate each of the plots, with clear comments where each plot is generated.
