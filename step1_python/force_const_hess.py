#!/usr/bin/env python3
"""
Usage:
  ./run

Executes the default force-constant calculation using the
settings defined in the run script.

Computes: cubic (Phi_ijk) and semi-diagonal quartic (Phi_iijj ≡ Phi_i i j j)
force constants in mass-weighted normal coordinates Q using finite differences
of Hessians from MLPs such as MACE-OFF/ PhysNet.

Conventions:
  Hessian is obtained in eV/Å^2 and converted to Eh/Bohr^2.

Inputs:
  - C2417905.inp (geometry file of the molecule to run VPT2 on)

Outputs:
  - optg_fromhessian.out: optimized geometry using MLP 
  - freqs_fromhessian.out: harmonic normal mode analysis performed on the optimized geometry
  - cubics_fromhessian.out  (unique i<=j<=k): cubic force constants
  - quartics_fromhessian.out (explicit i i j j lines for i<=j): quartic force constants
"""

import argparse
import itertools
from typing import Dict, List, Tuple

import numpy as np
from ase.io import read, write
from mace.calculators import mace_off
from ase.optimize import BFGS

# ---- constants ----
EV_PER_HARTREE = 27.211386246
BOHR_TO_ANG = 0.529177210903
AMU_TO_ME = 1822.888486209
AUCM = 219474.6313705

def principal_axis_frame(x, mass):
    """
    Place molecule in the principal axis frame.

    Parameters
    ----------
    x : (natom,3) ndarray
        Cartesian coordinates.
    mass : (natom,) ndarray
        Atomic masses.

    Returns
    -------
    x_pa : (natom,3) ndarray
        Coordinates in principal axis frame.
    eigvals : (3,) ndarray
        Principal moments of inertia.
    """

    mass = np.asarray(mass)
    x = np.asarray(x)

    # Center of mass
    com = np.sum(mass[:, None] * x, axis=0) / np.sum(mass)

    # Shift coordinates to COM
    xc = x - com

    # Inertia tensor
    I = np.zeros((3, 3))

    for m, r in zip(mass, xc):
        x_, y_, z_ = r

        I[0, 0] += m * (y_**2 + z_**2)
        I[1, 1] += m * (x_**2 + z_**2)
        I[2, 2] += m * (x_**2 + y_**2)

        I[0, 1] -= m * x_ * y_
        I[0, 2] -= m * x_ * z_
        I[1, 2] -= m * y_ * z_

    I[1, 0] = I[0, 1]
    I[2, 0] = I[0, 2]
    I[2, 1] = I[1, 2]

    # Diagonalize inertia tensor
    eigvals, eigvecs = np.linalg.eigh(I)

    # Rotate coordinates into principal-axis frame
    x_pa = xc @ eigvecs

    # Ensure right-handed coordinate system
    if np.linalg.det(eigvecs) < 0:
        eigvecs[:, 0] *= -1
        x_pa = xc @ eigvecs

    return x_pa, eigvals


def hessian_eh_per_bohr2(calc, atoms) -> np.ndarray:
    """
    Compute the Cartesian Hessian and convert it from
    eV/Å² to atomic units (Hartree/Bohr²).
    """

    H = calc.get_hessian(atoms)
    H = np.asarray(H, dtype=float).reshape((3 * len(atoms), 3 * len(atoms)))
    return (H / EV_PER_HARTREE) * (BOHR_TO_ANG ** 2)


def mass_vector_me(atoms) -> np.ndarray:
    """
    Construct the atomic mass vector in atomic units (electron masses).

    ASE provides atomic masses in atomic mass units (amu). For vibrational
    analysis, masses are converted to electron masses and are subsequently
    used to build the mass-weighted Hessian,

        Hmw = M^{-1/2} H M^{-1/2}

    and to transform between Cartesian and mass-weighted normal coordinates.
    """
    return np.asarray(atoms.get_masses(), dtype=float) * AMU_TO_ME


def diagonalize_mass_weighted_hessian(H_eh_bohr2: np.ndarray, m_me: np.ndarray):
    """
    Construct and diagonalize the mass-weighted Cartesian Hessian.

    The Cartesian Hessian H is transformed as

        Hmw = M^{-1/2} H M^{-1/2}

    Diagonalization of Hmw gives the harmonic normal modes and their
    corresponding squared angular frequencies,

        Hmw C = C Ω²

    where the columns of C are the mass-weighted normal-mode eigenvectors.

    The eigenvalues are converted to harmonic frequencies in cm^-1 using
    the atomic-unit conversion factor AUCM.
    """
    mcoord = np.repeat(m_me, 3)
    inv_sqrt_m = 1.0 / np.sqrt(mcoord)
    Hmw = inv_sqrt_m[:, None] * H_eh_bohr2 * inv_sqrt_m[None, :]

    omega2, eigvecs = np.linalg.eigh(Hmw)
    freq_cm = np.sign(omega2) * np.sqrt(np.abs(omega2)) * AUCM
    return omega2, eigvecs, freq_cm


def project_hessian_to_modes(Hx_eh_bohr2: np.ndarray, m_me: np.ndarray, Cmw: np.ndarray) -> np.ndarray:
    """
    Project the mass-weighted Cartesian Hessian into the
    harmonic normal-mode basis,    H_Q = Cmw^T * Hmw * Cmw
    where Hmw = M^{-1/2} Hx M^{-1/2}.
    """
    mcoord = np.repeat(m_me, 3)
    inv_sqrt_m = 1.0 / np.sqrt(mcoord)
    Hmw = inv_sqrt_m[:, None] * Hx_eh_bohr2 * inv_sqrt_m[None, :]
    return Cmw.T @ Hmw @ Cmw


class HQOracle:
    """
    Generates displaced geometries, evaluates MACE Hessians, projects
    them into the normal-mode basis, and caches results to avoid
    redundant Hessian evaluations.

    Q -> x -> Hessian -> projected H_Q(Q), cached by rounded Q.

      qmw = Cmw @ Q
      dx  = qmw / sqrt(mcoord)
      x_bohr = xref_bohr + dx
      x_ang  = x_bohr * BOHR_TO_ANG
    """

    def __init__(
        self,
        atoms,
        calc,
        xref_bohr: np.ndarray,
        Cmw: np.ndarray,
        m_me: np.ndarray,
        cache_digits: int = 12,
    ):
        self.atoms = atoms
        self.calc = calc
        self.xref_bohr = np.asarray(xref_bohr, dtype=float).reshape(-1)  # (3N,)
        self.Cmw = np.asarray(Cmw, dtype=float)  # (3N,nm)
        self.m_me = np.asarray(m_me, dtype=float)  # (N,)
        self.mcoord = np.repeat(self.m_me, 3)  # (3N,)
        self.cache_digits = cache_digits
        self.cache: Dict[bytes, np.ndarray] = {}

    def _key(self, Qv: np.ndarray) -> bytes:
        q = np.round(Qv.astype(float), self.cache_digits)
        return q.tobytes()

    def HQ(self, Qv: np.ndarray) -> np.ndarray:
        key = self._key(Qv)
        if key in self.cache:
            return self.cache[key]

        qmw = self.Cmw @ Qv
        dx = qmw / np.sqrt(self.mcoord)
        x_bohr = self.xref_bohr + dx
        x_ang = x_bohr.reshape(-1, 3) * BOHR_TO_ANG

        self.atoms.set_positions(x_ang)
        Hx = hessian_eh_per_bohr2(self.calc, self.atoms)
        HQ = project_hessian_to_modes(Hx, self.m_me, self.Cmw)

        self.cache[key] = HQ
        return HQ

# =========================================================================
# The cubic/quartic consistency checks and Hessian finite-difference
# force-constant assembly below are adapted from PyVPT2:
# https://github.com/philipmnel/pyvpt2
#
# Copyright 2021–2024 Philip Nelson
# Licensed under the BSD 3-Clause License.
# See LICENSES/PyVPT2_LICENSE.md.
def check_cubic(phi_ijk: np.ndarray, v_ind: List[int], tol: float = 1.0, verbose: bool = True) -> np.ndarray:
    """
    PyVPT2-style cubic sanity check: check permutation inconsistencies + symmetrize by averaging over 6 permutations.
    """
    no_inconsistency = True
    if verbose:
        print("Checking for numerical inconsistencies in cubic terms...")

    for unique_ijk in itertools.combinations_with_replacement(v_ind, 3):
        perms = list(set(itertools.permutations(unique_ijk, 3)))
        for a in range(len(perms)):
            for b in range(a + 1, len(perms)):
                ind1 = perms[a]
                ind2 = perms[b]
                diff = abs(phi_ijk[ind1] - phi_ijk[ind2])
                if diff >= tol:
                    if verbose:
                        print(ind1, ind2, diff)
                    no_inconsistency = False

    if verbose and no_inconsistency:
        print("No inconsistencies found")

    sym_phi = np.zeros_like(phi_ijk)
    for axes in itertools.permutations((0, 1, 2), 3):
        sym_phi += np.transpose(phi_ijk, axes)
    return sym_phi / 6.0


def check_quartic(phi_iijj: np.ndarray, v_ind: List[int], tol: float = 1.0, verbose: bool = True) -> np.ndarray:
    """
    PyVPT2-style quartic sanity check: enforce symmetry and optionally report asymmetry.
    """
    if verbose:
        print("Checking for numerical inconsistencies in quartic (iijj) terms...")

    no_inconsistency = True
    for i in v_ind:
        for j in v_ind:
            diff = abs(phi_iijj[i, j] - phi_iijj[j, i])
            if diff >= tol:
                if verbose:
                    print((i, j), (j, i), diff)
                no_inconsistency = False

    if verbose and no_inconsistency:
        print("No inconsistencies found")

    return 0.5 * (phi_iijj + phi_iijj.T)


def assemble_from_hessians(oracle: HQOracle, nm: int, h: float, v_ind: List[int]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Implement the Hessian-FD formula.
    Compute cubic and semi-diagonal quartic force constants from
    finite differences of projected Hessians in normal coordinates.
    Returns:
      phi_ijk : (nm,nm,nm)
      phi_iijj: (nm,nm)
    """
    phi_ijk = np.zeros((nm, nm, nm), dtype=float)
    phi_iijj = np.zeros((nm, nm), dtype=float)

    H0 = oracle.HQ(np.zeros(nm, dtype=float))

    Hp: Dict[int, np.ndarray] = {}
    Hn: Dict[int, np.ndarray] = {}

    for i in v_ind:
        Q = np.zeros(nm, dtype=float)
        Q[i] = +h
        Hp[i] = oracle.HQ(Q)
        Q[i] = -h
        Hn[i] = oracle.HQ(Q)

    # ---- cubics ----
    for i, j, k in itertools.product(v_ind, repeat=3):
        if i == j == k:
            phi_ijk[i, i, i] = (Hp[i][i, i] - Hn[i][i, i]) / (2.0 * h)
        else:
            val = (Hp[i][j, k] - Hn[i][j, k]) + (Hp[j][k, i] - Hn[j][k, i]) + (Hp[k][i, j] - Hn[k][i, j])
            phi_ijk[i, j, k] = val / (6.0 * h)

    # ---- quartics (semi-diagonal iijj) ----
    for i, j in itertools.product(v_ind, repeat=2):
        if i == j:
            phi_iijj[i, i] = (Hp[i][i, i] + Hn[i][i, i] - 2.0 * H0[i, i]) / (h**2)
        else:
            val = Hp[j][i, i] + Hn[j][i, i] + Hp[i][j, j] + Hn[i][j, j] - 2.0 * H0[i, i] - 2.0 * H0[j, j]
            phi_iijj[i, j] = val / (2.0 * (h**2))

    return phi_ijk, phi_iijj
# End of adaptation from PyVPT2
# =========================================================================


def main():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("-i", "--input", required=True, help="reference geometry XYZ (Angstrom)")
    ap.add_argument("--opt_out", default="optg_fromhessian.out", help="Optimized geometry file")
    ap.add_argument("--eps", type=float, default=0.2, help="step in Q-units")
    ap.add_argument("--freq_tol", type=float, default=5.0, help="drop modes with |freq| <= freq_tol (cm^-1)")
    ap.add_argument("--max_modes", type=int, default=None, help="optional cap on number of vibrational modes used")
    ap.add_argument("--cache_digits", type=int, default=12, help="rounding digits for Q cache keys")
    ap.add_argument("--freq_out", default="freqs_fromhessian.out", help="harmonic freq and eigenvector file")
    ap.add_argument("--cub_out", default="cubics_fromhessian.out", help="cubic output file")
    ap.add_argument("--quart_out", default="quartics_fromhessian.out", help="quartic output file (semi-diagonal iijj)")
    ap.add_argument("--cubic_tol", type=float, default=1.0, help="tolerance to report cubic permutation inconsistencies")
    ap.add_argument("--quartic_tol", type=float, default=1.0, help="tolerance to report quartic symmetry inconsistencies")
    args = ap.parse_args()
    
    # Load the reference molecular structure (Cartesian coordinates in Å).
    atoms = read(args.input, format="xyz")

    # Attach an ASE-compatible calculator. This can be replaced with any
    # ASE calculator providing energies, forces, and Hessians (e.g., PhysNet).
    calc = mace_off(model="large", default_dtype="float64", device='cpu')
    
    # Uncomment below to use the MACE-OFF23 medium model or any other ASE-compatible calculator of choice.
    
    # calc = mace_off(model="medium", default_dtype="float64", device='cpu')

    atoms.calc = calc

    # Optimize geometry using MACE-OFF
    optimizer = BFGS(atoms, logfile=None)
    optimizer.run(fmax=0.0001, steps=5000)
    x = atoms.get_positions() / BOHR_TO_ANG
    masses = atoms.get_masses() * AMU_TO_ME
    x_pa, p_inertia = principal_axis_frame(x, masses)
    atoms.set_positions(x_pa * BOHR_TO_ANG)
    print("Geometry optimization finished")
    write(args.opt_out, atoms, format="xyz")
    xref_bohr = atoms.get_positions().reshape(-1) / BOHR_TO_ANG

    # Reference Hessian for modes
    Hx0 = hessian_eh_per_bohr2(calc, atoms)
    m_me = mass_vector_me(atoms)
    _, eigvecs, freq_cm = diagonalize_mass_weighted_hessian(Hx0, m_me)

    with open(args.freq_out, "w") as f:
        f.write("Frequencies (cm^-1):\n")
        np.savetxt(
            f,
            np.column_stack([np.arange(1, len(freq_cm) + 1), freq_cm]),
            fmt=["%6d", "%16.8f"]
        )
        for i in range(1,len(freq_cm)+1):
            f.write("\n")
            f.write(" Mode " + str(i) + "\n")
            np.savetxt(
                f,
                eigvecs[:,i-1].reshape(len(atoms),3),
                fmt="%20.12f"
            )

    vib_idx = np.where(np.abs(freq_cm) > args.freq_tol)[0]
    if vib_idx.size == 0:
        raise RuntimeError(f"No modes with |freq| > {args.freq_tol} cm^-1. Try lowering --freq_tol.")
    if args.max_modes is not None:
        vib_idx = vib_idx[: args.max_modes]

    Cmw = eigvecs[:, vib_idx]  # (3N, nm)
    nm = Cmw.shape[1]
    v_ind = list(range(nm))

    oracle = HQOracle(atoms, calc, xref_bohr, Cmw, m_me, cache_digits=args.cache_digits)

    phi3, phi4 = assemble_from_hessians(oracle, nm=nm, h=args.eps, v_ind=v_ind)

    # Checks symmetrization
    phi3 = check_cubic(phi3, v_ind=v_ind, tol=args.cubic_tol, verbose=True)
    phi4 = check_quartic(phi4, v_ind=v_ind, tol=args.quartic_tol, verbose=True)

    # ---- write cubics unique i<=j<=k ----
    with open(args.cub_out, "w") as f:
        f.write("# cubic force constants from python MLPs via finite-difference of Hessians\n")
        f.write(f"# xyz={args.input}\n")
        f.write(f"# eps={args.eps} freq_tol={args.freq_tol} cm^-1\n")
        f.write(f"# nm (modes kept) = {nm}\n")
        f.write("#  i  j  k        Phi_ijk   (Hartree / Q^3)\n")
        for i in range(nm):
            for j in range(i, nm):
                for k in range(j, nm):
                    f.write(f"{i+1:4d} {j+1:4d} {k+1:4d}  {phi3[i, j, k]: .16e}\n")

    # ---- write quartics semi-diagonal explicitly as i i j j ----
    with open(args.quart_out, "w") as f:
        f.write("# semi-diagonal quartic force constants Phi_ijkl (i=i, k=l) from python MLPs via finite-difference of Hessians\n")
        f.write(f"# xyz={args.input}\n")
        f.write(f"# eps={args.eps} freq_tol={args.freq_tol} cm^-1\n")
        f.write(f"# nm (modes kept) = {nm}\n")
        f.write("#  i  j  k  l        Phi_ijkl   (Hartree / Q^4)\n")
        for i in range(nm):
            for j in range(i, nm):
                f.write(f"{i+1:4d} {i+1:4d} {j+1:4d} {j+1:4d}  {phi4[i, j]: .16e}\n")

    print(f"Wrote cubics:   {args.cub_out}")
    print(f"Wrote quartics: {args.quart_out}")
    print(f"Modes kept (nm) = {nm}")
    print(f"Cached HQ evaluations = {len(oracle.cache)}")


if __name__ == "__main__":
    main()

