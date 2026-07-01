#!/usr/bin/env python3
"""
Usage: python gvpt2.py --mol_name={MOLECULE} --num_modes={N_MODES}
                       (optional) --enable_intensity

MOLECULE:
    Prefix of the harmonic frequency, optimized geometry, cubic, and quartic files
    (e.g., C2417905, cholesterol, C583788).

N_MODES:
    Number of vibrational modes included in the GVPT2 calculation.
    This should match the number of modes retained during force-field
    generation after excluding translational and rotational modes.
    e.g. 3 modes for water (3N-6). 

If "--enable_intensity" flag is present, double harmonic IR/Raman intensities
will be calculated; otherwise only the energies of fundamentals are computed

Inputs:
  - All outputs from step 1: e.g.{MOLECULE}_nma, opt, cubics, quartics.out

Outputs: Harmonic and GVPT2 frequencies, IR/ Raman intensities: {MOLECULE}_intensity.out 

Computes: VPT2 fundamentals + GVPT2 (Fermi) resonance treatment using a small effective Hamiltonian + IR/Raman intensities.

Workflow:
  1. Compute double-harmonic IR and Raman intensities using MACE-MDP.
  2. Read harmonic frequencies and cubic/quartic force constants.
  3. Convert force constants from Q to dimensionless q coordinates.
  4. Identify Fermi-I and Fermi-II resonances.
  5. Compute deperturbed anharmonic constants (DVPT2).
  6. Build and diagonalize resonant polyads (GVPT2).
  7. Write harmonic/GVPT2 frequencies and IR/Raman intensities.

Conventions:
- omega_cm: cm^-1
- Phi3_Q, Phi4_Q read from files: assumed Eh in Q-units (mass-weighted normal coords Q)
- Convert Q -> dimensionless q:
    q_i = sqrt(omega_i) * Q_i   with omega in Eh
    phi3(q) = Phi3(Q)/sqrt(wi wj wk)
    phi4(q) = Phi4(Q)/sqrt(wi wj wk wl)
- Convert phi(q) from Eh -> cm^-1
- Compute chi (cm^-1) using omega (cm^-1) and phi (cm^-1)

Resonances:
- Identify Fermi-1: 2*omega_i ~ omega_k
- Identify Fermi-2: omega_i + omega_j ~ omega_k
- Deperturb chi denominators when resonance flagged (DVPT2)
- Build and diagonalize GVPT2 effective Hamiltonians for resonant polyads
- No Darling-Dennison resonance.
"""

import re
import math
import itertools
import sys
import argparse
from pathlib import Path
from typing import Dict, Tuple, List, Set
from scipy.optimize import linear_sum_assignment
from ase.io import read
from ase import units

import numpy as np
from dataclasses import dataclass

CM_PER_EH = 219474.6313705  # cm^-1 per Hartree
EH_PER_CM = 1.0 / CM_PER_EH
BOHR_TO_ANG = 0.529177210903
AMU_TO_ME = 1822.888486209

# =========================================================================
# The State, Interaction, Polyad, and fermi_solver implementations below
# are adapted from PyVPT2: https://github.com/philipmnel/pyvpt2
#
# Copyright 2021–2024 Philip Nelson
# Licensed under the BSD 3-Clause License.
# See LICENSES/PyVPT2_LICENSE.md.
@dataclass
class State:
    state: tuple
    nu: float

@dataclass
class Interaction:
    left: State
    right: State
    phi: float
    ftype: int  # 1 or 2

def fermi_solver(interaction_list: List[Interaction]) -> Dict[tuple, float]:
    polyad_list = []
    for interaction in interaction_list:
        flag = False
        for polyad in polyad_list:
            if interaction.left.state in polyad.state_list:
                polyad.add(interaction)
                flag = True
            elif interaction.right.state in polyad.state_list:
                polyad.add(interaction)
                flag = True

        if flag is False:
            polyad_list.append(Polyad(interaction))

    state_list: Dict[tuple, float] = {}
    for polyad in polyad_list:
        state_list.update(polyad.solve())

    return state_list

class Polyad:
    def __init__(self, interaction: Interaction):
        left = interaction.left
        right = interaction.right
        self.state_list = set([left.state, right.state])
        self.nu_list = {left.state: left.nu, right.state: right.nu}
        self.phi_list = {(left.state, right.state): (interaction.phi, interaction.ftype)}

    def add(self, interaction: Interaction):
        left = interaction.left
        right = interaction.right
        self.state_list.update([left.state, right.state])
        self.nu_list.update({left.state: left.nu})
        self.nu_list.update({right.state: right.nu})
        self.phi_list.update({(left.state, right.state): (interaction.phi, interaction.ftype)})

    def build_hamiltonian(self):
        self.state_list_enum = {state: i for i, state in enumerate(self.state_list)}
        dim = len(self.state_list_enum.keys())
        self.H = np.zeros((dim, dim), dtype=float)

        for state, i in self.state_list_enum.items():
            self.H[i, i] = self.nu_list[state]

        for states, interaction in self.phi_list.items():
            i = self.state_list_enum[states[0]]
            j = self.state_list_enum[states[1]]
            phi, ftype = interaction

            # Keep exactly your scaling rules
            if ftype == 1:
                self.H[i, j] = 1.0 / 4.0 * phi
                self.H[j, i] = self.H[i, j]
            elif ftype == 2:
                self.H[i, j] = 1.0 / (np.sqrt(2.0) * 2.0) * phi
                self.H[j, i] = self.H[i, j]

    def solve(self) -> Dict[tuple, float]:
        self.build_hamiltonian()
        evals, evecs = np.linalg.eigh(self.H)

        print("\n===================================================")
        print("GVPT2 Polyad Prints Eigen values and vecs")
        print("===================================================")

        # Basis states in matrix order
        basis_states = [None] * len(self.state_list_enum)
        for state, idx in self.state_list_enum.items():
            basis_states[idx] = state

        print("Basis states:")
        for i, st in enumerate(basis_states):
            print(f"  {i}: {st}")

        print("\nEigenvalues and CI coefficients:")
        for root in range(len(evals)):
            print(f"\nRoot {root+1}: {evals[root]:12.6f} cm^-1")
            print(f"                    CI coef      coef**2")
            for i, st in enumerate(basis_states):
                coef = evecs[i, root]
                print(f"   {st!s:12s}  {coef:10.6f}   {coef**2:10.6f}")
#        print("GVPT2 Polyad Printing ends")
        print("===================================================")
         
        # Assign eigenvalues to “most-overlapping” basis state
        # inds = [np.argmax(vec) for vec in np.square(evecs.T)]
        # eval_dict = dict(zip(inds, evals))
        # freqs = {state: eval_dict[ind] for (state, ind) in self.state_list_enum.items()}
        overlaps = np.square(evecs)   # basis x eigenstate
        cost = -overlaps # -ve sign to maximize the overlap
        row_ind, col_ind = linear_sum_assignment(cost)
        inv_map = {i: state for state, i in self.state_list_enum.items()}
        freqs = {inv_map[i]: evals[j] for i, j in zip(row_ind, col_ind)}

        return freqs
# End of adaptation from PyVPT2
# ============================================================


# ============================================================
#  Parsing helpers for reading normal modes, QFF. 
# ============================================================

def read_freqs_and_vecs(path: str, n_modes: int) -> List[float]:
    text = Path(path).read_text()

    # parse Frequencies block to get harmonic frequencies
    header_pat = re.compile(r"^\s*Frequencies\s*\(cm\^-1\)\s*:\s*$", re.MULTILINE)
    matches = list(header_pat.finditer(text))

    if not matches:
        raise ValueError(f'Could not find "Frequencies (cm^-1):" block in {path}')

    start = matches[-1].end()
    tail = text[start:]

    float_re = r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?"
    line_pat = re.compile(rf"^\s*(\d+)\s+({float_re})\s*$")

    freqs: List[float] = []
    for line in tail.splitlines():
        if not line.strip():
            if freqs:
                break
            continue
        m = line_pat.match(line)
        if not m:
            if freqs:
                break
            continue
        val = float(m.group(2).replace("D", "E").replace("d", "E"))
        freqs.append(val)

    if len(freqs) < n_modes:
        raise ValueError(f"Found only {len(freqs)} freqs in last block, need {n_modes}.")

    # parse normal mode vectors
    mode_pat = re.compile(r"^\s*mode\s+(\d+)\s*$", re.IGNORECASE)
    n_blocks = len(freqs)

    lines = text.splitlines()
    blocks = [None] * n_blocks

    i = 0
    while i < len(lines):
        m = mode_pat.match(lines[i])

        if m:
            mode_idx = int(m.group(1)) - 1

            block = []
            j = i + 1
    
            while j < len(lines):
    
                line = lines[j].strip()
    
                # next mode starts
                if mode_pat.match(line):
                    break
    
                if line:
                    vals = [
                        float(x.replace("D", "E").replace("d", "E"))
                        for x in line.split()
                    ]
    
                    if len(vals) != 3:
                        raise ValueError(
                            f"Expected 3 numbers in mode {mode_idx+1}, "
                            f"got {len(vals)}"
                        )
                    block.append(vals)
    
                j += 1
    
            blocks[mode_idx] = np.asarray(block, dtype=float)
            i = j
        else:
            i += 1

    return freqs[-n_modes:], np.array(blocks[-n_modes:])


def parse_cubic_tensor(path: str, n_modes: int) -> Dict[Tuple[int, int, int], float]:
    text = Path(path).read_text()
    phi: Dict[Tuple[int, int, int], float] = {}

    line_re = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+(?:[Ee][-+]?\d+)?)\s*$"
    )

    for line in text.splitlines():
        m = line_re.match(line)
        if not m:
            continue
        i, j, k = map(int, m.group(1, 2, 3))
        val = float(m.group(4))
        for p in set(itertools.permutations([i, j, k], 3)):
            phi[p] = val

    for i in range(1, n_modes + 1):
        for j in range(1, n_modes + 1):
            for k in range(1, n_modes + 1):
                phi.setdefault((i, j, k), 0.0)

    return phi


def parse_quartic_tensor(path: str, n_modes: int) -> Dict[Tuple[int, int], float]:
    text = Path(path).read_text()
    phi4: Dict[Tuple[int, int], float] = {}

    line_re = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([-+]?\d*\.\d+(?:[Ee][-+]?\d+)?)\s*$"
    )

    for line in text.splitlines():
        m = line_re.match(line)
        if not m:
            continue
        i, j, k, l = map(int, m.group(1, 2, 3, 4))
        val = float(m.group(5))
        for p in set(itertools.permutations([i, k], 2)):
            phi4[p] = val

    for i in range(1, n_modes + 1):
        for k in range(1, n_modes + 1):
            phi4.setdefault((i, k), 0.0)

    return phi4


#=======================================================
# Rotational constants and Coriolis coupling
#=======================================================
def rotational_constants(x, masses):
    """
    Compute rotational constants (A, B, C) in atomic units.

    Parameters
    ----------
    x : (natom, 3) array
        Cartesian coordinates in principal axis frame (bohr).
    masses : (natom,) array
        Atomic masses in electron masses (atomic units).

    Returns
    -------
    A, B, C : floats
        Rotational constants in Hartree (since ħ = 1).
    """

    x = np.asarray(x)
    m = np.asarray(masses)

    x2 = x[:, 0]
    y2 = x[:, 1]
    z2 = x[:, 2]

    I_x = np.sum(m * (y2**2 + z2**2))
    I_y = np.sum(m * (x2**2 + z2**2))
    I_z = np.sum(m * (x2**2 + y2**2))

    # rotational constants in atomic units (ħ = 1)
    A = 1.0 / (2.0 * I_x)
    B = 1.0 / (2.0 * I_y)
    C = 1.0 / (2.0 * I_z)

    return np.array([A, B, C])

def coriolis_zeta(vec):
    """
    Compute Coriolis coupling constants ζ^α_{kl}.

    Parameters
    ----------
    vec : ndarray, shape (nmode, natom, 3)
        Mass-weighted normal modes:
        vec[k, i, a] with a = 0(x),1(y),2(z)

    Returns
    -------
    zeta : ndarray, shape (3, nmode, nmode)
        zeta[0] = ζ^x, zeta[1] = ζ^y, zeta[2] = ζ^z
    """

    vec = np.asarray(vec)
    nmode, natom, _ = vec.shape

    zeta = np.zeros((3, nmode, nmode))

    # unpack Cartesian components for clarity
    vx = vec[:, :, 0]  # (k,i)
    vy = vec[:, :, 1]
    vz = vec[:, :, 2]

    for k in range(nmode):
        for l in range(nmode):
            # x-component
            zeta[0, k, l] = np.sum(vy[k, :] * vz[l, :] - vz[k, :] * vy[l, :])

            # y-component
            zeta[1, k, l] = np.sum(vz[k, :] * vx[l, :] - vx[k, :] * vz[l, :])

            # z-component
            zeta[2, k, l] = np.sum(vx[k, :] * vy[l, :] - vy[k, :] * vx[l, :])

    return zeta


#=======================================================
# Convert cubic and quartic force-constant tensors from
# mass-weighted normal coordinates (Q) to dimensionless
# normal coordinates (q) used in GVPT2 calculations.
# Q -> dimensionless 
#=======================================================

def convert_tensors_Q_to_dimensionless_q(
    omega_eh: List[float],
    phi3_Q: Dict[Tuple[int, int, int], float],
    phi4_Q: Dict[Tuple[int, int], float],
) -> Tuple[Dict[Tuple[int, int, int], float], Dict[Tuple[int, int], float]]:
    sqrtw = [math.sqrt(w) for w in omega_eh]

    phi3_q: Dict[Tuple[int, int, int], float] = {}
    for (i, j, k), val in phi3_Q.items():
        phi3_q[(i, j, k)] = val / (sqrtw[i - 1] * sqrtw[j - 1] * sqrtw[k - 1]) * CM_PER_EH

    phi4_q: Dict[Tuple[int, int], float] = {}
    for (i, k), val in phi4_Q.items():
        phi4_q[(i, k)] = val / (sqrtw[i - 1] ** 2 * sqrtw[k - 1] ** 2) * CM_PER_EH

    return phi3_q, phi4_q


# Resonance detection (heuristic) for Fermi type 1 and 2.
def identify_fermi_resonances(
    omega_cm: List[float],
    phi3_cm: Dict[Tuple[int, int, int], float],
    *,
    omega_thresh: float,
    K_thresh: float,
    v_ind: List[int],
) -> Set[Tuple[int, Tuple[int, int]]]:
    fermi: Set[Tuple[int, Tuple[int, int]]] = set()

    # Fermi-1: 2*i ≈ k
    for i, k in itertools.permutations(v_ind, 2):
        d_omega = abs(2.0 * omega_cm[i - 1] - omega_cm[k - 1])
        if d_omega <= omega_thresh:
            phi = phi3_cm[(i, i, k)]
            dK = (phi ** 4) / (256.0 * (d_omega ** 3 + 1e-30))
            if dK >= K_thresh:
                fermi.add((k, (i, i)))

    # Fermi-2: i + j ≈ k
    for i, j, k in itertools.permutations(v_ind, 3):
        d_omega = abs(omega_cm[i - 1] + omega_cm[j - 1] - omega_cm[k - 1])
        if d_omega <= omega_thresh:
            phi = phi3_cm[(i, j, k)]
            dK = (phi ** 4) / (64.0 * (d_omega ** 3 + 1e-30))
            if dK >= K_thresh:
                fermi.add((k, (i, j)))

    return fermi

# ===================================================================
# Deperturbation (DVPT2) Compute deperturbed χij anharmonic constants
# prior to resonance treatment.
# ===================================================================

def delta_ijk(wi: float, wj: float, wk: float) -> float:
    return ((wi + wj + wk) *
            (wi + wj - wk) *
            (wi - wj + wk) *
            (wi - wj - wk))


def compute_chi_deperturbed_cm(
    omega_cm: List[float],
    phi3_cm: Dict[Tuple[int, int, int], float],
    phi4_cm: Dict[Tuple[int, int], float],
    rotc: np.ndarray,
    zeta: np.ndarray,
    *,
    fermi_set: Set[Tuple[int, Tuple[int, int]]],
    v_ind: List[int],
    eps_div: float = 1e-30,
) -> List[List[float]]:
    n = len(omega_cm)
    chi = [[0.0] * n for _ in range(n)]

    # diagonal chi_ii
    for i in v_ind:
        wi = omega_cm[i - 1]
        term = phi4_cm[(i, i)]

        s = 0.0
        for k in v_ind:
            wk = omega_cm[k - 1]
            phi = phi3_cm[(i, i, k)]

            if (k, (i, i)) in fermi_set:
                # safe form (no 4 wi^2 - wk^2)
                s += -0.5 * (phi ** 2) * (1.0 / (2.0 * wi + wk + eps_div) + 4.0 / (wk + eps_div))
            else:
                num = (8.0 * wi ** 2 - 3.0 * wk ** 2) * (phi ** 2)
                den = wk * (4.0 * wi ** 2 - wk ** 2)
                s += num / (den + eps_div)

        chi[i - 1][i - 1] = (term - s) / 16.0

    # off-diagonal chi_ij
    for i in v_ind:
        wi = omega_cm[i - 1]
        for j in v_ind:
            if i == j:
                continue

            wj = omega_cm[j - 1]
            term = phi4_cm[(i, j)]

            # - Σ_k phi_iik phi_jjk / wk   (convention)
            s1 = 0.0
            for k in v_ind:
                wk = omega_cm[k - 1]
                s1 += (phi3_cm[(i, i, k)] * phi3_cm[(j, j, k)]) / (wk + eps_div)

            # + Σ_k phi_ijk^2 * delta_ij
            s2 = 0.0
            for k in v_ind:
                wk = omega_cm[k - 1]
                phi = phi3_cm[(i, j, k)]

                if (k, (i, j)) in fermi_set:
                    # i + j = k : drop 1/(wi+wj-wk)
                    delta_ij = (
                        1.0 / (wi + wj + wk + eps_div) +
                        1.0 / (-wi + wj + wk + eps_div) +
                        1.0 / (wi - wj + wk + eps_div)
                    ) / (-2.0)

                elif (i, (j, k)) in fermi_set:
                    # j + k = i : drop 1/(-wi+wj+wk)
                    delta_ij = (
                        1.0 / (wi + wj + wk + eps_div) +
                        1.0 / (wi + wj - wk + eps_div) +
                        1.0 / (wi - wj + wk + eps_div)
                    ) / (-2.0)

                elif (j, (i, k)) in fermi_set:
                    # i + k = j : drop 1/(wi-wj+wk)
                    delta_ij = (
                        1.0 / (wi + wj + wk + eps_div) +
                        1.0 / (wi + wj - wk + eps_div) +
                        1.0 / (-wi + wj + wk + eps_div)
                    ) / (-2.0)

                else:
                    D = delta_ijk(wi, wj, wk)
                    delta_ij = 2.0 * wk * (wi ** 2 + wj ** 2 - wk ** 2) / (D + eps_div)

                s2 += (phi ** 2) * delta_ij

            # + 4(wi**2+wj**2) / wiwj sum Be^t zeta^t_{ij}**2
            s3 = 0.0
            for k in range(3):
                s3 += 4.0*(wi**2+wj**2) / (wi*wj) * rotc[k] * zeta[k,i-1,j-1]**2
                
            chi[i - 1][j - 1] = (term - s1 + s2 + s3) / 4.0

    return chi


def fundamentals_from_chi_cm(omega_cm: List[float], chi_cm: List[List[float]]) -> List[float]:
    """
    Nondegenerate deperturbed fundamentals.
    """
    n = len(omega_cm)
    nu = []
    for i in range(n):
        offsum = sum(chi_cm[i][j] for j in range(n) if j != i)
        nu.append(omega_cm[i] + 2.0 * chi_cm[i][i] + 0.5 * offsum)
    return nu


# ==============================================================
#  Build |i>, |ii>, |ij> energies
# ==============================================================

def compute_chi0_cm(
    omega_cm: List[float],
    phi3_cm: Dict[Tuple[int, int, int], float],
    phi4_cm: Dict[Tuple[int, int], float],
    rotc: np.ndarray,
    zeta: np.ndarray,
    *,
    fermi_set: Set[Tuple[int, Tuple[int, int]]],
    v_ind: List[int],
    eps_div: float = 1e-30,
) -> float:
    """
    Compute chi0 contribution in cm^-1.
    nondegenerate part chi0:

      chi0 += phi_iijj[i,i]
      chi0 -= (7/9) * phi_iii^2 / omega_i
      chi0 += 3 * omega_i * phi_ijj^2 / (4 omega_j^2 - omega_i^2)
      chi0 += 2 * phi_ijk^2 * delta_0   for i < j < k
      chi0 -= 16 * B_alpha * (1 + 2 * sum zeta_alpha,ij^2)

    followed by:

      chi0 /= 64
    """
    chi0 = 0.0

    # One-mode and two-mode non-rotational terms.
    for i in v_ind:
        ii = i - 1
        wi = omega_cm[ii]

        chi0 += phi4_cm[(i, i)]
        chi0 -= (7.0 / 9.0) * (phi3_cm[(i, i, i)] ** 2) / (wi + eps_div)

        for j in v_ind:
            if j == i:
                continue
            jj = j - 1
            wj = omega_cm[jj]
            den = 4.0 * (wj ** 2) - (wi ** 2)
            chi0 += 3.0 * wi * (phi3_cm[(i, j, j)] ** 2) / (den + eps_div)

    # Three-mode terms, counted only once as i < j < k.
    for i, j, k in itertools.combinations(v_ind, 3):
        wi = omega_cm[i - 1]
        wj = omega_cm[j - 1]
        wk = omega_cm[k - 1]

        if (k, (i, j)) in fermi_set:
            # i + j = k: drop 1/(wi + wj - wk) from delta_0.
            delta_0 = (
                1.0 / (wi + wj + wk + eps_div)
                - 1.0 / (wi - wj + wk + eps_div)
                - 1.0 / (-wi + wj + wk + eps_div)
            )
        elif (i, (j, k)) in fermi_set:
            # j + k = i: drop 1/(-wi + wj + wk) from delta_0.
            delta_0 = (
                1.0 / (wi + wj + wk + eps_div)
                - 1.0 / (wi + wj - wk + eps_div)
                - 1.0 / (wi - wj + wk + eps_div)
            )
        elif (j, (i, k)) in fermi_set:
            # i + k = j: drop 1/(wi - wj + wk) from delta_0.
            delta_0 = (
                1.0 / (wi + wj + wk + eps_div)
                - 1.0 / (wi + wj - wk + eps_div)
                - 1.0 / (-wi + wj + wk + eps_div)
            )
        else:
            D = delta_ijk(wi, wj, wk)
            delta_0 = -8.0 * wi * wj * wk / (D + eps_div)

        chi0 += 2.0 * (phi3_cm[(i, j, k)] ** 2) * delta_0

    # rotational terms
    for k in range(3):
        sr = 0.0
        for [i,j] in itertools.combinations(v_ind, 2):
            sr += zeta[k, i-1, j-1]**2
        chi0 -= 16.0 * rotc[k] * (1.0 + 2.0*sr)
   
    chi0 /= 64.0

    return chi0


def compute_zpve_cm(
    omega_cm: List[float],
    chi_cm: List[List[float]],
    v_ind: List[int],
    *,
    chi0_cm: float = 0.0,
) -> Tuple[float, float, float]:
    """
    Compute harmonic and anharmonic ZPVE in cm^-1.

      ZPVE_harm = 1/2 * sum_i omega_i

      ZPVE_anh = chi0
                 + 1/2 * sum_i (omega_i + 1/2 * chi_ii)
                 + 1/4 * sum_{i<j} chi_ij

    Here chi0_cm can include the non-rotational chi0 contribution.
    It does not include the rotational/Coriolis part.

    Returns
    -------
    harmonic_zpve_cm, anharmonic_zpve_cm, correction_cm
    """
    harmonic_zpve = 0.0
    anharmonic_zpve = chi0_cm

    for i in v_ind:
        ii = i - 1
        harmonic_zpve += 0.5 * omega_cm[ii]
        anharmonic_zpve += 0.5 * (omega_cm[ii] + 0.5 * chi_cm[ii][ii])

    for i, j in itertools.combinations(v_ind, 2):
        ii = i - 1
        jj = j - 1
        anharmonic_zpve += 0.25 * chi_cm[ii][jj]

    correction = anharmonic_zpve - harmonic_zpve

    return harmonic_zpve, anharmonic_zpve, correction


def build_explicit_state_energies(
    omega_cm: List[float],
    chi_cm: List[List[float]],
    v_ind: List[int],
) -> Tuple[List[float], Dict[Tuple[int, int], float], Dict[Tuple[int, int], float]]:
    """
    Build deperturbed state energies in the nondegenerate, g=0 limit.
    With 1-based mode labels in v_ind and 0-based Python storage:

      nu_i      = omega_i + 2 chi_ii + 1/2 * sum_{j != i} chi_ij
      overtone  = 2 omega_i + 6 chi_ii +     sum_{j != i} chi_ij
      band_ij   = omega_i + omega_j + 2 chi_ii + 2 chi_jj + 2 chi_ij
                  + 1/2 * sum_{k != i,j} (chi_ik + chi_jk)
    Equations from: J. Phys. Chem. A 2021, 125, 1301−1324. 

    Returns
    -------
    Tuple containing:
      - fundamentals list (0-based storage)
      - overtone dict keyed by (i,i) with 1-based labels
      - combination-band dict keyed by (i,j) with i<j and 1-based labels
    """
    n = len(omega_cm)
    nu = [0.0] * n
    overtone: Dict[Tuple[int, int], float] = {}
    band: Dict[Tuple[int, int], float] = {}

    for i in v_ind:
        ii = i - 1
        nu_i = omega_cm[ii] + 2.0 * chi_cm[ii][ii]
        ov_i = 2.0 * omega_cm[ii] + 6.0 * chi_cm[ii][ii]
        for j in v_ind:
            if j == i:
                continue
            jj = j - 1
            nu_i += 0.5 * chi_cm[ii][jj]
            ov_i += chi_cm[ii][jj]
        nu[ii] = nu_i
        overtone[(i, i)] = ov_i

    for i, j in itertools.combinations(v_ind, 2):
        ii = i - 1
        jj = j - 1
        bij = (
            omega_cm[ii]
            + omega_cm[jj]
            + 2.0 * chi_cm[ii][ii]
            + 2.0 * chi_cm[jj][jj]
            + 2.0 * chi_cm[ii][jj]
        )
        for k in v_ind:
            if k == i or k == j:
                continue
            kk = k - 1
            bij += 0.5 * (chi_cm[ii][kk] + chi_cm[jj][kk])
        band[(i, j)] = bij

    return nu, overtone, band


# ============================================================
#  GVPT2: build Interaction list from fermi_set and solve
# ============================================================

def apply_gvpt2_fermi(
    omega_cm: List[float],
    nu_deperturbed: List[float],  # fundamentals ν_i (cm^-1), 0-based list
    chi_cm: List[List[float]],
    phi3_cm: Dict[Tuple[int, int, int], float],
    fermi_set: Set[Tuple[int, Tuple[int, int]]],
    v_ind: List[int],
) -> Tuple[List[float], Dict[tuple, float]]:
    """
    Build the resonant polyads (fundamental <-> overtone/combination),
    diagonalize each polyad Hamiltonian, update only fundamentals (|i>) and return.

    Returns:
      updated_nu (list, cm^-1)
      state_energy_map (all states solved in polyads)
    """
    _, overtone, comb = build_explicit_state_energies(omega_cm, chi_cm, v_ind)

    interactions: List[Interaction] = []

    for (k, parents) in sorted(fermi_set):
        # k is the fundamental mode index (1-based) being resonant with parents
        if parents[0] == parents[1]:
            # Fermi-1: 2*i ~ k
            i = parents[0]
            left = State(state=(k,), nu=nu_deperturbed[k - 1])
            right = State(state=(i, i), nu=overtone[(i, i)])
            phi = phi3_cm[(k, i, i)]
            interactions.append(Interaction(left=left, right=right, phi=phi, ftype=1))
        else:
            # Fermi-2: i + j ~ k
            i, j = parents
            if i == j:
                continue
            ii, jj = (i, j) if i < j else (j, i)
            left = State(state=(k,), nu=nu_deperturbed[k - 1])
            right = State(state=(ii, jj), nu=comb[(ii, jj)])
            phi = phi3_cm[(k, ii, jj)]
            interactions.append(Interaction(left=left, right=right, phi=phi, ftype=2))

    if not interactions:
        return nu_deperturbed, {}

    state_map = fermi_solver(interactions)

    # Update only fundamentals
    nu_updated = nu_deperturbed.copy()
    for state, energy in state_map.items():
        if len(state) == 1:
            nu_updated[state[0] - 1] = energy

    return nu_updated, state_map


# ============================================================
#  Debug helper
# ============================================================

def stats(name, d):
    vals = [abs(v) for v in d.values() if v != 0.0]
    if not vals:
        print(name, "no nonzero values")
        return
    vals_sorted = sorted(vals)
    print(
        name,
        "count =", len(vals),
        "min =", vals_sorted[0],
        "median =", vals_sorted[len(vals_sorted)//2],
        "max =", vals_sorted[-1],
    )


# ============================================================
#  Main
# ============================================================

def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--mol_name", required=True, help='Name of the molecule, i.e., the string before _ in *.out files')
    ap.add_argument("--num_modes", default=None, type=int, help='Number of normal modes, default 3*N_atoms-6')
    ap.add_argument("--enable_intensity", action='store_true', help='Calculate double harmonic intensities')
    args = ap.parse_args()

    xyz_file = args.mol_name + "_opt.out"
    monomer_log = args.mol_name + "_nma.out"
    cubic_file = args.mol_name + "_cubics.out"
    quartic_file = args.mol_name + "_quartics.out"

    # Step 0: setting up the molecule
    # 0.1) Read in Cartesian coordinates, harmonic freqs and normal modes
    atoms = read(xyz_file, index=0, format="xyz")
    n_atoms = len(atoms)
    if args.num_modes == None:
        n_modes = 3*n_atoms-6
    else:
        n_modes = args.num_modes

    omega_cm, modes = read_freqs_and_vecs(monomer_log, n_modes=n_modes)

    # 0.2) Compute rotational constants and Coriolis coupling
    x = atoms.get_positions() / BOHR_TO_ANG
    masses = atoms.get_masses() * AMU_TO_ME
    rotc = rotational_constants(x, masses) * CM_PER_EH
    print("Rotational constants (cm-1):")
    print(rotc)

    zeta = coriolis_zeta(modes)

    # Step 1 (optional): Double-harmonic IR/Raman intensity calculation
    if args.enable_intensity:
        from mace.calculators.mace import MACECalculator

        # =================================================================
        # This section below that computes double harmonic intensities is 
        # adapted from MACE-MDP: https://github.com/Nilsgoe/MACE-MDP
        # 
        # Copyright [MACE-MDP] is © 2026, [Nils Gönnheimer]
        # MACE-MDP is distributed under the Academic Software License v1.0
        # (ASL) for academic non-commercial use.
        # See LICENSES/MACE-MDP_LICENSE.md
        intens_out = args.mol_name + "_intensity.out"
        # 1.1) setup the calculator for dipole+polarizability
        mdp_path = "./MACE-MDP.model"
        mdp_calc = MACECalculator(
            model_paths=mdp_path,
            model_type="DipolePolarizabilityMACE",
            device='cpu',
            default_dtype="float64",
        )

        atoms = atoms.copy()

        # 1.2) Unweight normal mode vectors to get Cartesian displacements
        masses = atoms.get_masses()
        modes *= masses[np.newaxis, :, np.newaxis] ** -0.5

        # 1.3) IR/Raman stick intensity using dipole/polarizability derivatives
        ir_intens = []
        total_intens = [] # Raman total intensity
        iso_intens = [] # Raman isotropic intensity
        aniso_intens = [] # Raman anisotropic intensity
        
        pos0 = atoms.get_positions().copy()
        # Finite-difference displacement for dipole/polarizability derivatives
        step = 1.0e-3
        
        for mode_idx, mode in enumerate(modes):
            dpos = mode * step
        
            atoms.set_positions(pos0 + dpos)
            mu_plus = mdp_calc.get_property("dipole", atoms)
            alpha_plus = mdp_calc.get_property("polarizability", atoms)
        
            atoms.set_positions(pos0 - dpos)
            mu_minus = mdp_calc.get_property("dipole", atoms)
            alpha_minus = mdp_calc.get_property("polarizability", atoms)
        
            atoms.set_positions(pos0)
        
            dmu_dq = (mu_plus - mu_minus) / (2.0 * step)
            dalpha_dq = (alpha_plus - alpha_minus) / (2.0 * step)
        
            a_iso = np.trace(dalpha_dq) / 3.0
            a_ani_sq = 0.5 * (
                (dalpha_dq[0, 0] - dalpha_dq[1, 1]) ** 2
                + (dalpha_dq[1, 1] - dalpha_dq[2, 2]) ** 2
                + (dalpha_dq[2, 2] - dalpha_dq[0, 0]) ** 2
                + 6.0 * (dalpha_dq[0, 1] ** 2 + dalpha_dq[1, 2] ** 2 + dalpha_dq[2, 0] ** 2)
            )
        
            if omega_cm[mode_idx] < 5.0:
                i_ir = 0.0
                i_iso = 0.0
                i_ani = 0.0
                i_tot = 0.0
            else:
                i_ir = np.linalg.norm(dmu_dq) ** 2
                i_iso = 45.0 * (a_iso ** 2)
                i_ani = 7.0 * a_ani_sq
                i_tot = i_iso + i_ani
        
            ir_intens.append(i_ir)
            iso_intens.append(i_iso)
            aniso_intens.append(i_ani)
            total_intens.append(i_tot)
        
        ir_intens = np.array(ir_intens)
        iso_intens = np.array(iso_intens)
        aniso_intens = np.array(aniso_intens)
        total_intens = np.array(total_intens)
        # End of adaptation from MACE-MDP
        #==================================================================

    # ============================================================
    # Step 2: Anharmonic frequency calculation (DVPT2/GVPT2)
    # ============================================================
    # Resonance thresholds (tune)
    FERMI_OMEGA_THRESH = 100.0  # cm^-1
    FERMI_K_THRESH = 1.0

    # 2.1) Convert omega to Eh for Q->q scaling ONLY
    omega_eh = [w * EH_PER_CM for w in omega_cm]

    # 2.2) Parse Phi(Q) in Eh/Q^3, Eh/Q^4
    phi3_Q = parse_cubic_tensor(cubic_file, n_modes=n_modes)
    phi4_Q = parse_quartic_tensor(quartic_file, n_modes=n_modes)

    print("\nRAW TENSOR MAGNITUDES (before any scaling)")
    stats("Phi3_Q", phi3_Q)
    stats("Phi4_Q", phi4_Q)
    print("omega_cm min/max:", min(omega_cm), max(omega_cm))

    # 2.3) Q -> dimensionless q and in cm-1
    phi3_q_cm, phi4_q_cm = convert_tensors_Q_to_dimensionless_q(omega_eh, phi3_Q, phi4_Q)

    print("\nRAW TENSOR MAGNITUDES (after Q->q scaling, still Eh)")
    stats("Phi3_q_cm", phi3_q_cm)
    stats("Phi4_q_cm", phi4_q_cm)

    v_ind = list(range(1, n_modes + 1))
    # 2.4) Detect Fermi resonances (for deperturbation + GVPT2)
    fermi = identify_fermi_resonances(
        omega_cm,
        phi3_q_cm,
        omega_thresh=FERMI_OMEGA_THRESH,
        K_thresh=FERMI_K_THRESH,
        v_ind=v_ind,
    )

    print(f"\nDetected {len(fermi)} Fermi candidates:")
    for (k, parents) in sorted(fermi):
        if parents[0] == parents[1]:
            i = parents[0]
            print(f"  (k={k}, (i,i)=({i},{i}))  |2ω_i-ω_k|={abs(2*omega_cm[i-1]-omega_cm[k-1]):.3f} cm^-1")
        else:
            i, j = parents
            print(f"  (k={k}, (i,j)=({i},{j}))  |ω_i+ω_j-ω_k|={abs(omega_cm[i-1]+omega_cm[j-1]-omega_cm[k-1]):.3f} cm^-1")

    # 2.5a) Ordinary VPT2 chi: no deperturbation
    chi_vpt2_cm = compute_chi_deperturbed_cm(
        omega_cm, phi3_q_cm, phi4_q_cm,
        rotc, zeta,
        fermi_set=set(),   # no resonant terms removed
        v_ind=v_ind
    )

    nu_vpt2 = fundamentals_from_chi_cm(omega_cm, chi_vpt2_cm)
    
    # 2.5b) DVPT2 chi (cm^-1)
    chi_cm = compute_chi_deperturbed_cm(
        omega_cm, phi3_q_cm, phi4_q_cm,
        rotc, zeta,
        fermi_set=fermi, v_ind=v_ind
    )

    # 2.6) Deperturbed fundamentals (Explicit nondegenerate expression)
    nu_depert = fundamentals_from_chi_cm(omega_cm, chi_cm)

    # 2.7) Zero-point vibrational energy (ZPVE) with chi0
    chi0_cm = compute_chi0_cm(
        omega_cm=omega_cm,
        phi3_cm=phi3_q_cm,
        phi4_cm=phi4_q_cm,
        rotc=rotc,
        zeta=zeta,
        fermi_set=fermi,
        v_ind=v_ind,
    )

    harmonic_zpve_cm, anharmonic_zpve_cm, zpve_corr_cm = compute_zpve_cm(
        omega_cm=omega_cm,
        chi_cm=chi_cm,
        v_ind=v_ind,
        chi0_cm=chi0_cm,
    )
#    for freq in nu_depert:
#        print(freq)

    # 2.8) GVPT2: variational mixing in resonant polyads
    nu_gvpt2, state_map = apply_gvpt2_fermi(
        omega_cm=omega_cm,
        nu_deperturbed=nu_depert,
        chi_cm=chi_cm,
        phi3_cm=phi3_q_cm,
        fermi_set=fermi,
        v_ind=v_ind
    )

    # -----------------------------
    # Print summary
    # -----------------------------
    print("\nHarmonic ω (cm^-1):")
    for i, w in enumerate(omega_cm, start=1):
        print(f"  ω_{i:2d} = {w:12.6f}")

    print("\nDeperturbed χ_ii (cm^-1):")
    for i in range(n_modes):
        print(f"  χ_{i+1:02d}{i+1:02d} = {chi_cm[i][i]: .8f}")

    print("\nFundamentals (cm^-1):")
    print("  mode     Harmonic     VPT2      DVPT2      GVPT2       GVPT2-HOshift   GVPT2-VPT2shift")
    for i in range(1, n_modes + 1):
        dv = nu_depert[i - 1]
        vp = nu_vpt2[i - 1]
        gv = nu_gvpt2[i - 1]
        ho = omega_cm[i - 1]
        print(f"  {i:4d}  {ho:10.4f}  {vp:10.4f}  {dv:10.4f}  {gv:10.4f}  {gv-ho:10.4f}  {gv-vp:10.4f}")

    if state_map:
        print("\nGVPT2 solved state energies (cm^-1) in resonant polyads (includes combo/overtones):")
        for st, e in sorted(state_map.items(), key=lambda x: (len(x[0]), x[0])):
            print(f"  {st!s:12s}  {e:12.6f}")

    print("\nZero-Point Vibrational Energy:")
    print(f"  chi0 = {chi0_cm: .8f} cm^-1")
    print("  unit          harmonic ZPVE      correction       anharmonic ZPVE")

    unit_list = [
        ("cm^-1", 1.0),
        ("kcal/mol", 0.002859144),
        ("kJ/mol", 0.011962656),
    ]

    for unit, factor in unit_list:
        print(
            f"  {unit:8s}  "
            f"{factor * harmonic_zpve_cm:14.6f}  "
            f"{factor * zpve_corr_cm:14.6f}  "
            f"{factor * anharmonic_zpve_cm:14.6f}"
        )

    print("\nDone.\n")
    print("GVPT2 code run is now complete")

    if args.enable_intensity:
        print("HO/GVPT2 freq +IR/Raman intensities written to: ", intens_out, ".")    
        with open(intens_out, "w") as f:
            f.write("# HO_freq  GVPT2  IR_intens  Raman_intens  iso_intens  aniso_intens\n")
            for i in range(n_modes):
                # Vibrational modes begin after the first 6 external modes
                f.write(f"{omega_cm[i]:7.2f}  {nu_gvpt2[i]:7.2f}  {ir_intens[i]:.8f}  {total_intens[i]:.8f}  {iso_intens[i]:.8f}  {aniso_intens[i]:.8f}\n")

if __name__ == "__main__":
    main()
