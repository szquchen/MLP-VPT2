# MLP-VPT2

## About
The software performs VPT2 calculation for the anharmonic vibrational frequencies and optionally IR and Raman intensities using machine-learned potentials.

## Usage
The VPT2 calculation involves two steps.
The first is to obtain the cubic and quartic force constants (2 options for this step: Fortran or Python MLPs), and the
second step is to compute VPT2 energies for fundamentals, as well as
IR and Raman intensities.

### (1A) MLPs written in Fortran: Get the QFF using analytical gradients (folder step1_fortran)

(a) Modify the template `pes_shell_template.f90` and rename it to `pes_shell.f90`. The three key subroutines/functions need to be implemented are `pes_init()`, `pot(x)`, and `gradient(x)`. An example for a PIP potential of aspirin is given. Then compile the program using the Makefile (with necessary modifications). **Warning: the aspirin example may take ~1 hour to compile.**

(b) Run the executable. For the aspirin example
```bash
./vpt2.x aspirin.inp
```
         
         The input file contains Cartesian coordinates of atoms in standard xyz format, in Angstrom. This executable not only computes QFF constants (written to XXX.out files) but also performs a normal VPT2          calculation, i.e., without considering any resonances. You may get negative (imaginary) frequencies for certain modes. To explicitly  treat resonances, copy XXX.out to step2 run DVPT2 and GVPT2, see below.

(1B) MLPs written in Python: Get the QFF using Python-based MLP (folder step1_python):

     (a) The sample force_const_hess.py has been set up to use MACE-OFF (the 23 large model). You need to install MACE first. To use other potential surface, modify force_const_hess.py, particularly, the ASE calculator.

     (b) Modify the script "run" to adjust the input arguments. The detailed explanations for these arguments can be found at the beginning of "main()" in force_const_hess.py

     (c) Then simply execute
         ./run

(2) Run DVPT2+GVPT2 and (optionally) compute IR/Raman intensities in folder step2:

    Copy the XXX.out files from step1 to step2 Modify the script "run" to adjust the input arguments (explained at the beginning of "main" in gvpt2.py. Then simply execute 
    ./run
    
    If the "--enable_intensity" argument is present, the double harmonic  intensities of IR and Raman are calculated using MACE-MDP model for dipole and polarizability. Otherwise, the intensity will be skipped and only the energies of fundamental transitions are calculated.

## Third-Party Software

The Python implementation of the finite-difference force-constant
assembly and GVPT2 resonance-handling algorithms is adapted from
PyVPT2:
https://github.com/philipmnel/pyvpt2

Copyright 2021–2024 Philip Nelson

Licensed under the BSD 3-Clause License.
See LICENSES/PyVPT2_LICENSE.md

The intensity calculations interface with MACE-MDP and use code
adapted from MACE-MDP repository:
https://github.com/Nilsgoe/MACE-MDP

Copyright [MACE-MDP] is © 2026, [Nils Gönnheimer]

MACE-MDP is distributed under the **Academic Software License v1.0 (ASL)** for
academic non-commercial use. See LICENSES/MACE-MDP_LICENSE.md

## Contributors

Chen Qu: Step1 Fortran

Saikiran Kotaru: Step1 Python and Step2

## Citation

If you use this software in your research, please cite:

Kotaru, S.; Qu, C.; Nandi, A.; Houston, P. L.; Bowman, J. M.
VPT2 Calculations of Vibrational Energies of CH3COOC6H4COOH Done in Seconds on a Laptop Using a Machine Learned Potential.
J. Phys. Chem. Lett. 2026, 17, 6580–6586.
DOI: https://doi.org/10.1021/acs.jpclett.6c01186

## Related References

This software builds upon established VPT2/GVPT2 methodologies and incorporates concepts from the following software packages and references.

Nelson, P. M.; Sherrill, C. D.
pyVPT2: Interoperable Software for Anharmonic Vibrational Frequency Calculations.
J. Chem. Phys. 2025, 162, 032501.
DOI: https://doi.org/10.1063/5.0251445

Barone, V.
Anharmonic Vibrational Properties by a Fully Automated Second-Order Perturbative Approach.
J. Chem. Phys. 2005, 122, 014108.
DOI: https://doi.org/10.1063/1.1824881

Franke, P. R.; Stanton, J. F.; Douberly, G. E.
How to VPT2: Accurate and Intuitive Simulations of CH Stretching Infrared Spectra Using VPT2+K with Large Effective Hamiltonian Resonance Treatments.
J. Phys. Chem. A 2021, 125, 1301–1324.
DOI: https://doi.org/10.1021/acs.jpca.0c09526

Gönnheimer, N.; Reuter, K.; Kapil, V.; Margraf, J. T.
MACE-MDP: A General Dipole and Polarizability Model for Organic Molecules and Materials.
ChemRxiv 2026.
DOI: https://doi.org/10.26434/chemrxiv.15000716/v2

