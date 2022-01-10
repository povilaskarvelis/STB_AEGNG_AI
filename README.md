# STB_AEGNG_AI

This repository contains the code for reproducing the results in "_A computational model of hopelessness and active-escape bias in suicidality_"


To reproduce model simulations run _Dashboard_AI_GoNoGo_paper.m_ 

- _GoNoGo_EA.m_ specifies the generative model and the generative process strucutres
- _spm_MDP_VB_LC_EA.m_ implements active inference and learning
- _plot_GoNoGo_EA_stats.m_ produces plots for a single task run
- _plot_GoNoGo_EA_many_many.m_ produces plots of more robust results by averaging relevant variables across repeated task runs
