# MacAyeal_glacial

This repo contains code relevant for the following manuscript submitted to Philosophical Transactions of the Royal Society A: Ajagun-Brauns J and Ditlevsen P (2025) A minimal conceptual model for glacial-interglacial cycles Phil. Trans. R. Soc.

Functions.py contains all custom functions used in the numerical integration scripts:

MacAyeal_model_extended.py - contains numerical integration of MacAyeal's model with the configuration described in Section 4 of the manuscript, for sinusoidal solar forcing, Q. This code was used to produce Fig 2.

MacAyeal_modulatedforcing.py - contains essentially the same as the previous code but rather than ramping down the timescale ratio, R at t=250kyr, the amplitude of the modelled obliquity solar forcing is ramped down to zero and the amplitude of a different solar forcing is ramped up - sinusoidal solar forcing with period T=20kyr (precessional timescales) and modulated by another sinusoidal forcing with period T=100kyr

