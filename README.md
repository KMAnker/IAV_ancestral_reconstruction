tester


# ScattBO Benchmark - Bayesian optimisation for materials discovery

A self-driving laboratory (SDL) is an autonomous platform that conducts machine learning (ML) selected experiments to achieve a user-defined objective. An objective can be to synthesise a specific material.[1] Such an SDL will synthesise a material, evaluate if this is the target material and if necessary optimise the synthesis parameters for the next synthesis. One way to evaluate if the material is the target material is by measuring scattering data and comparing that to the scattering pattern of the target material. However, these types of SDLs can be expensive to run, which means that intelligent experimental planning is essential. At the same time, only a few people have access to an SDL for materials discovery. Therefore, it is currently challenging to benchmark Bayesian optimisation algorithms for experimental planning tasks in SDLs.

Here, we present a Python-based benchmark (ScattBO) that is an in silico simulation of an SDL for materials discovery. Based on a set of synthesis parameters, the benchmark ‘synthesises’ a structure, calculates the scattering pattern[2] and compares this to the scattering pattern of the target structure. Note: Scattering data may not be enough to conclusively validate that the target material has been synthesised.[3] The benchmark can include other types of data as long they can be simulated.

**Synthesis parameters**:
  - pH (float):       The pH value, which scales the size of the structure. Range: [0, 14]
  - pressure (float): The pressure value, which controls the lattice constant of the structure. Range: [0, 100]
  - solvent (int):    The solvent type, which determines the structure type for large clusters. 
                      0 for 'Ethanol', 1 for 'Methanol', 2 for 'Water', 3 for 'Others'

**Features**:
  - Continous parameters (pH and pressure) and discrete parameter (solvent).
  - Two sizes of benchmark is provided: 1) ScatterBO_small_benchmark and 2) ScatterBO_large_benchmark.
  - Possibility of two different objectives: 1) The scattering data in Q-space (Sq) or 2) the scattering data in r-space (Gr).
  - Possibility of multi-objective optimisation using both Sq and Gr.
  - Possibility of four target scattering patterns: 1) simulated Sq, 2) simulated Gr, 3) experimental Sq, and 4) experimental Gr.
  - OBS: The scattering pattern calculations can be slow on CPU. **It is recommended to use GPU.**

# Scoreboard

<p align="center"><i>These scoreboards represent the performance of different BO algorithms on various of the ScattBO benchmarks. If you have new scores to report, feel free to contact us.</i></p>

## Small benchmark 


<sup>1</sup> Criteria for convergence for simulated Sq data is defined as the number of steps until Rwp < 0.04.<br>
<sup>2</sup> Criteria for convergence for simulated Gr data is defined as the number of steps until Rwp < 0.04.<br>
<sup>3</sup> Criteria for convergence for experimental Sq data is defined as the number of steps until Rwp < 0.84.<br>
<sup>4</sup> Criteria for convergence for experimental Gr data is defined as the number of steps until Rwp < 0.80.<br>
<sup>5</sup> Criteria for convergence for multi-objective optimisation is defined for both above criteria.

# Usage
See (https://github.com/AndySAnker/ScattBO/tree/main/tools) for examples of single-objective optimisation with [Dragonfly](https://github.com/dragonfly/dragonfly/tree/master) or [skopt](https://scikit-optimize.github.io/stable/auto_examples/bayesian-optimization.html).

## Example usage with [Dragonfly](https://github.com/dragonfly/dragonfly/tree/master)


## Visualise results
```python

# Use the parameters that minimize the function as input to the ScatterBO_small_benchmark function
ScatterBO_small_benchmark((pH, pressure, solvent), plot=True, simulated_or_experimental='simulated', scatteringfunction='Gr')

```


## Reporting issues

If you encounter any issues or problems with our software, please report them by opening an issue on our GitHub repository. Please include as much detail as possible, including steps to reproduce the issue and any error messages you received.

## Seeking support

If you need help using our software, please reach out to us on our GitHub repository. We'll do our best to assist you and answer any questions you have.

# References

See also the [Twitter submission post](https://twitter.com/SodeAndy/status/1773474538631651769)
