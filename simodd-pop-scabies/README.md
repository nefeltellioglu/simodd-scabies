# simodd-pop-scabies

Slightly modified version of simodd-pop, a computational simulation model of population with demography.

For further detail, please see Geard et al. [Synthetic population dynamics: a model of household demography](http://jasss.soc.surrey.ac.uk/16/1/8.html). JASSS 16(1):8.

## Index:

`population` contains python classes for individuals, populations, and a simulation object.

`observers` contains observers for storing data on a population during a simulation, written to a HD5 file.

`example_pop` contains parameters and data files (contemporary Australian population) for a simple testcase.

## Deprecated:

`pop_explore` directory contains various functions for analysing and plotting population model output; however, this is deprecated and will be removed.

## Example:

To run the test case, from `example_pop`, execute `python main.py`

A small population will be created and evolved over time, with output written to `output/population.hd5`
 





