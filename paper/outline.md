# xTB Paper Outline

## Introduction
* Literature review
* TB context
* Brief set-up


## Motivation
* Can mention constraints of DFTB parameterisation and molecular xTB parameterisation schemes.
* We should also acknowledge the more recent papers: GFN2, periodic xTB library (that links with DFTB+) and any issues
  they originally flagged up with the parameterisation.


## Theory:
* Periodic expression for H0

* Brief introduction/review to the approximations of the H2 term and MNOK potential history
    * Hit the key points of the approximation: Integral, monopole charges
    * Choice of parameterisation scheme for the potential
    * Known issues in the literature with MNOK 

* Our treatment of H2: 
  * Subtract the non-converging terms from the real-space sum and treat with an Ewald-like expansion in respective 
    short-range (real-space) and long-range (Fourier space), respectively. See Williams.
  * Issue w.r.t present of eta as a shell-resolved quantity leads to charged terms. Charge corrections not available for n=3 term.
  * Rui’s mathematical treatment for n=3. 


## Approximations
* Does the choice of cut-off in Rui’s approach affect the results? 
  I.e. we should show an example test case, then run for a range of systems with min and max cut-offs, and show its effect
  on the total energy
  * We should compare this in a few cases to CP2K, and their truncation of the potential

* Show for the simplest system possible, what goes wrong in the H2 potential. 
  * Show this mathematically - as done by Rui (or put in the theory section and reference back to it)
  * Give a toy/simple system that demonstrates what occurs in practice - maybe the charge is the quantity to focus on
    * There is somewhat of a question regarding how we show this
  

## Results
* For a large range of binary materials, show the trend away from agreement with experiment as systems become more ionic: 
  * As Fred said, hopefully see a trend as systems become more ionic. 
  * Demonstrate with a scatter plot
      * Quantify w.r.t. a) lattice constants, b) overall error in compression curve (see the Delta definition - add this)
    or c) band gaps - don't necessarily expect this one to correlate as strongly, but would be interesting to see
  * High throughput would be beneficial. 
  
* From this, can we make any general statements about what class of materials we expect the current parameterisation for the H2 term 
perform reasonably?
  * We have to define what is reasonable - perhaps 5% deviation in the equilibrium lattice constant w.r.t. experiment
    Or be cheeky and define a smaller cut-off w.r.t. DFT, which might already be ~ 3% off experiment.

* Define a class of binary materials for which xTB gives reasonable lattice constants, compression curves, and perhaps
band gaps

* Show compression curves for materials where we define it as being valid, along with quoted band gaps.
  * Important narrative here is to show all comparisons for materials, which can be compared to DFTB+, and DFT
  * xTB needs to be comparable to DFTB in accuracy, but offer a larger materials phase space.
  
* Band structures: If yes, restrict it to VBT and CBB, for direct gaps. 

* Show an initial result for the MOFS - Sets up the next paper nicely. 
