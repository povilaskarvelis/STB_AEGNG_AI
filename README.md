# STB_AEGNG_AI

This repository contains the code for reproducing the results in Karvelis, P., & Diaconescu, A. O. (2022). A Computational Model of Hopelessness and Active-Escape Bias in Suicidality. Computational Psychiatry, 6(1), 34–59. DOI: http://doi.org/10.5334/cpsy.80


To reproduce model simulations run _Dashboard_AI_GoNoGo_paper.m_ 

- _GoNoGo_EA.m_ specifies the generative model and the generative process strucutres
- _spm_MDP_VB_LC_EA.m_ implements active inference and learning
- _plot_GoNoGo_EA_stats.m_ produces plots for a single task run
- _plot_GoNoGo_EA_many_many.m_ produces plots of more robust results by averaging relevant variables across repeated task runs
