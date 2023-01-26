Constraining Dark matter models with PINT
====

Project objectives
------------------

* What
The aim of the project is to investigate the possibility of obtaining new most stringent constraints 
on dark matter coupling constants using a modern pulsar timing package called PINT.

* How
We do this by observing pulsar binaries, whose orbital periods can be time-dependent as a result of 
dark matter-binary pulsar interaction. This effect on orbital period can be incorporated into a 
timing model, and as a result the coupling constant of a direct interaction can be estimated by PINT 
when applying the least square minimization procedure to fit the timing model parameters to the data.

* Assumption
Although the interaction can be both gravitational and direct, the former one is expected to be too weak 
for having significant imprints in measured data. Therefore, we only include into our analysis direct interaction and 
aim to use data to set upper bounds on direct interaction constants.

* Current status and objectives
This project is a work in progress and this repository includes the last modifications of PINT's code,
so it accounts for the dark matter effects on binary systems. Current objectives are:

** To modify the BT binary model (Blandford and Teukolsky, 1976), so it accounts for time-dependent period (plus derived quntities, including the Kepler equation). Variations of other orbital parameters are not considered.

** To constraint a coupling constant of universal direct interaction between scalar DM-ordinary matter. You can check it `here <https://arxiv.org/abs/1612.06789/>`_ and `here <https://arxiv.org/abs/1910.08544/>`_ .





