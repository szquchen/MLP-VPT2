# MLP-VPT2

## About
The software performs VPT2 calculation for the anharmonic vibrational frequencies using machine-learned potentials.

## Usage
The VPT2 calculation involves two steps.
The first is to obtain the cubic and quartic force constants (2 options for this step: Fortran or Python MLPs), and the
second step is to compute VPT2 energies for fundamentals, as well as
IR and Raman intensities.

(1A) MLPs written in Fortran: Get the QFF constants using PIP analytical gradients (folder step1):
     (a) Modify the function fq() and gq() in step1/normal.f90
         using your routine from your pes_shell.f90, and compile

     (b) Prepare the input file. An example is given in oxalate.inp.
         The file begins with Cartesian coordinates (in Angstrom,
         in principal axis frame if Coriolis coupling terms are
         included in the calculation. For large molecules, it is
         fine to ignore these terms, and this code does not include
         Coriolis terms), then 3N-6 normal mode eigenvectors
         corresponding to this Cartesian geometry.

     (c) Run the executable
         ./vpt2.x {MOLENAME} {N_MODES}
         For the aspirin example
         ./vpt2.x aspirin 57
         This executable not only computes QFF constants (written to
         MOLENAME_XXX.out files) but also performs a "standard" VPT2
         calculation, i.e., without considering any resonances. You may
         get negative (imaginary) frequencies for certain modes.

(1B) MLPs written in Python: Get the QFF constants using MACE (or other MLP) in Python
     (folder step1_python)
     (a) If you haven't installed MACE (or other MLP) yet, create a
         Python environment and follow the installation instructions
         from the corresponding MLP

     (b) The sample force_const_hess.py has been set up to use MACE-OFF
         (the 23 large model) and MACE-MDP for dipole and polarizability.
         To use other potential/dipole/polarizability surface, modify
         force_const_hess.py, particularly, the ASE calculator.

     (c) Modify the script "run" to adjust the input arguments. The detailed
         explanations for these arguments can be found at the beginning of
         "main()" in force_const_hess.py

     (d) Then simply execute
         ./run

(2) Run DVPT2+GVPT2 and compute IR/ Raman intensities in folder step2:
    Copy the *.out files from step1 to the folder "step2"
    Modify the script "run" to adjust the input arguments (explained at the
    beginning of "main" in gvpt2.py. Then simply execute 
    ./run

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

## Contribution

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

