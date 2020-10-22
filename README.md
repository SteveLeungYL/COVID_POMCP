# COVID_POMCP

This is the release of POMCP software used in the AAMAS 2021 paper

- **Liang, Y.**, & Yadav, A. (2021). Let the DOCTOR Decide Whom to Test: Adaptive Testing Strategies to Tackle the COVID-19 Pandemic. (Submitted to AAMAS-2021)

The POMCP software are originated from "Online Monte-Carlo Planning in Large POMDPs"
by David Silver and Joel Veness. 

## Plots from the paper:
Plots are stored in repository: https://github.com/SteveLeungYL/Plots_from_the_COVID_POMCP.git 

## Error information:
For unknown reason, when clone the repository to a new location, the file /src/Makefile will always turn out to be corrupted. ``make clean`` and ``make`` command would not work. The solution is to copy the raw file of /src/Makefile in the Github repository, paste it and overwrite the Makefile in the local file system, and rerun the ``make clean``, ``make`` again.