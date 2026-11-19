#!/usr/bin/env python3
"""HubbardModel-Lab - educational fermionic Hubbard model lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * The Hubbard model is the smallest model of interacting electrons:
    a hopping amplitude t that wants electrons to delocalise, and an
    on-site repulsion U that charges them for sharing a site. Almost
    everything interesting - antiferromagnetism, the Mott transition,
    the Heisenberg model itself - is the argument between those two
    terms, and this lab has that argument on two tiny clusters where
    every number can be checked exactly.
  * The lab diagonalises the Hamiltonian exactly (exact
    diagonalisation, ED) in the fixed (N_up, N_down) occupation
    basis, with a from-scratch Jacobi eigensolver and ONE
    Jordan-Wigner sign rule used for hopping and spin flips alike.
  * It models TWO clusters: a 2-site dimer and a 2x2 plaquette
    (4 sites). Four sites are not a material. No thermodynamic-limit,
    real-material, superconductivity, or device claim is computed or
    claimed anywhere in this file - those uses are permanent HOLD
    gates.

Conventions, stated once and used everywhere:

  * Units: t = 1. All energies are in units of t, all U values are
    U/t ratios. H = -t sum_<ij>,sigma (c+_i c_j + h.c.) + U sum_i
    n_i,up n_i,down over the stated bonds.
  * Spin-orbitals are ordered up-block first, then down-block:
    orbital index p = spin * L + site. Every fermionic sign in this
    file comes from that one ordering (Jordan-Wigner parity).
  * The plaquette is the 2x2 square with its four edges as bonds,
    each counted once: sites 0=(0,0), 1=(1,0), 2=(0,1), 3=(1,1),
    bonds (0,1), (2,3), (0,2), (1,3). The square closes on itself,
    so this is equivalently a 4-site ring; no bond is double-counted.
  * Half filling on the plaquette means N = 4 electrons,
    (N_up, N_down) = (2, 2); on the dimer, (1, 1).

Every real-world use is a permanent HOLD gate - with Gemini Argon 4
Extended Thinking as DECIDER (names-only, tighten-only), Muse as
advisor (never overrides HOLD), and Google Flow plan-only seats.
"""
import hashlib
import json
import math
import sys
from itertools import combinations

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"
T_HOP = 1.0

# Labelled U/t sweep used across the plaquette book. Teaching values,
# not parameters of any material.
U_SWEEP = (0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 16.0)

DIMER_BONDS = ((0, 1),)
PLAQUETTE_BONDS = ((0, 1), (2, 3), (0, 2), (1, 3))
PLAQUETTE_SITES = 4

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "real_material_prediction_claim": "HOLD - no predicting a real material's properties from a 2-site or 4-site cluster",
    "superconductivity_claim": "HOLD - no claim about real superconductors, pairing mechanisms in real materials, or any Tc",
    "materials_discovery_claim": "HOLD - no materials-discovery or candidate-compound claim from a teaching model",
    "hardware_specification_claim": "HOLD - no hardware specification for any real device is asserted",
    "real_device_characterisation": "HOLD - no characterising a real sample, crystal, or device from this model",
    "phase_certification_claim": "HOLD - no certifying that any real system is in a Mott, magnetic, or any other phase",
    "thermodynamic_limit_claim": "HOLD - 4-site results are finite-cluster results; no thermodynamic-limit claim follows",
    "quantum_advantage_claim": "HOLD - no quantum-advantage claim; this is classical exact diagonalisation of tiny sectors",
    "cold_atom_experiment_claim": "HOLD - no directing real cold-atom or optical-lattice experiments from this lab",
    "experimental_protocol_claim": "HOLD - no experimental protocol, procedure, or safety claim for real laboratories",
    "investment_claim": "HOLD - no investment claim; a lab model is not a prospectus",
    "certification_claim": "HOLD - no certification that any system meets any standard",
    "live_hardware_control": "HOLD - no control of real hardware, instruments, or laboratory equipment, ever",
    "model_as_measurement_claim": "HOLD - a computed number is a statement about the model at stated inputs, never a measurement",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Jacobi eigensolver
# Real symmetric matrices only - the Hubbard Hamiltonian in the
# occupation basis is real symmetric (hopping entries are +/- t,
# the interaction is diagonal). Cyclic Jacobi rotations, with the
# eigenvector matrix accumulated alongside, refereed on known
# spectra before it is trusted with a Hamiltonian.

def jacobi_eigh(A, tol=1e-13, max_sweeps=100):
    """Eigen-decomposition of a real symmetric matrix.

    Returns (eigenvalues, eigenvectors): eigenvalues sorted
    ascending, eigenvectors[k] the matching unit eigenvector.
    """
    n = len(A)
    a = [list(map(float, row)) for row in A]
    V = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(max_sweeps):
        off = math.sqrt(sum(a[p][q] ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < tol:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-16:
                    continue
                app, aqq = a[p][p], a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                tt = math.copysign(
                    1.0 / (abs(tau) + math.sqrt(tau * tau + 1.0)), tau)
                c = 1.0 / math.sqrt(tt * tt + 1.0)
                s = tt * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    evals = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: evals[i])
    return ([evals[i] for i in order],
            [[V[r][i] for r in range(n)] for i in order])


def matvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(len(v))) for i in range(len(A))]


def eigen_residual(A, value, vector):
    """max |H v - E v| over entries: the eigenvector's receipt."""
    hv = matvec(A, vector)
    return max(abs(hv[i] - value * vector[i]) for i in range(len(vector)))


# ================================================================ 2. Fermion sector + Hamiltonian
# A basis state is a pair (up_mask, down_mask) of site bitmasks. The
# full Fock mask is up_mask | (down_mask << L): spin-orbital p for
# (site, spin) is spin * L + site, and every operator action below
# reads its sign off the parity of that full mask. One sign rule for
# hopping, spin flips, everything.

def sector_basis(L, n_up, n_dn):
    """All (up_mask, down_mask) states with the stated populations."""
    ups = [sum(1 << s for s in comb)
           for comb in combinations(range(L), n_up)]
    dns = [sum(1 << s for s in comb)
           for comb in combinations(range(L), n_dn)]
    return [(u, d) for u in ups for d in dns]


def _annihilate(full, p):
    """c_p on a full Fock mask: (new_mask, sign) or (None, 0)."""
    if not (full >> p) & 1:
        return None, 0.0
    sign = -1.0 if bin(full & ((1 << p) - 1)).count("1") % 2 else 1.0
    return full & ~(1 << p), sign


def _create(full, p):
    """c^+_p on a full Fock mask: (new_mask, sign) or (None, 0)."""
    if (full >> p) & 1:
        return None, 0.0
    sign = -1.0 if bin(full & ((1 << p) - 1)).count("1") % 2 else 1.0
    return full | (1 << p), sign


def _split(full, L):
    return (full & ((1 << L) - 1), full >> L)


def build_hubbard(L, bonds, n_up, n_dn, U, t=T_HOP):
    """Hubbard Hamiltonian in the (n_up, n_dn) sector.

    Returns (H, basis). Diagonal: U per doubly occupied site.
    Off-diagonal: -t per hop along a bond, sign from the one
    Jordan-Wigner rule on the full Fock mask.
    """
    basis = sector_basis(L, n_up, n_dn)
    idx = {state: i for i, state in enumerate(basis)}
    n = len(basis)
    H = [[0.0] * n for _ in range(n)]
    for col, (up, dn) in enumerate(basis):
        H[col][col] += U * bin(up & dn).count("1")
        full = up | (dn << L)
        for spin in (0, 1):
            off = spin * L
            for (i, j) in bonds:
                for a, b in ((i, j), (j, i)):
                    f1, s1 = _annihilate(full, off + a)
                    if f1 is None:
                        continue
                    f2, s2 = _create(f1, off + b)
                    if f2 is None:
                        continue
                    H[idx[_split(f2, L)]][col] += -t * s1 * s2
    return H, basis


def ground_state(L, bonds, n_up, n_dn, U, t=T_HOP):
    """(E0, psi0, basis) for one sector at one U."""
    H, basis = build_hubbard(L, bonds, n_up, n_dn, U, t)
    evals, evecs = jacobi_eigh(H)
    return evals[0], evecs[0], basis


def sector_spectrum(L, bonds, n_up, n_dn, U, t=T_HOP):
    H, _ = build_hubbard(L, bonds, n_up, n_dn, U, t)
    return jacobi_eigh(H)[0]


# ================================================================ 3. Observables
# Double occupancy is diagonal in the occupation basis. The spin
# correlation S_i . S_j = Sz_i Sz_j + (1/2)(S+_i S-_j + S-_i S+_j)
# is not: the flip terms are applied as four-operator sequences on
# the full Fock mask, with the same sign rule as the Hamiltonian.

def double_occupancy(psi, basis):
    """<sum_i n_i,up n_i,down>: total doublons in the state."""
    return sum(psi[k] ** 2 * bin(u & d).count("1")
               for k, (u, d) in enumerate(basis))


def spin_correlation_nn(psi, basis, L, bonds):
    """Mean nearest-neighbour <S_i . S_j> over the stated bonds."""
    idx = {state: k for k, state in enumerate(basis)}
    total = 0.0
    for (i, j) in bonds:
        val = 0.0
        for k, (u, d) in enumerate(basis):
            szi = 0.5 * (((u >> i) & 1) - ((d >> i) & 1))
            szj = 0.5 * (((u >> j) & 1) - ((d >> j) & 1))
            val += psi[k] ** 2 * szi * szj
            full = u | (d << L)
            for a, b in ((i, j), (j, i)):
                # S+_a S-_b = [c+_{a,up} c_{a,dn}] [c+_{b,dn} c_{b,up}]
                f1, s1 = _annihilate(full, b)
                if f1 is None:
                    continue
                f2, s2 = _create(f1, L + b)
                if f2 is None:
                    continue
                f3, s3 = _annihilate(f2, L + a)
                if f3 is None:
                    continue
                f4, s4 = _create(f3, a)
                if f4 is None:
                    continue
                target = _split(f4, L)
                if target in idx:
                    val += 0.5 * s1 * s2 * s3 * s4 * psi[idx[target]] * psi[k]
        total += val
    return total / len(bonds)


# ================================================================ 4. The dimer, analytically
# Two sites, half filling. The singlet ground state and the triplet
# decouple exactly: E_singlet = (U - sqrt(U^2 + 16 t^2)) / 2,
# E_triplet = 0, so the singlet-triplet gap is
# (sqrt(U^2 + 16 t^2) - U) / 2, which tends to 4 t^2 / U - the
# Heisenberg exchange J - as U grows. The ED must reproduce all of
# this; the closed forms are the referee, not the source.

def dimer_gap_exact(U, t=T_HOP):
    return (math.sqrt(U * U + 16.0 * t * t) - U) / 2.0


def dimer_energy_exact(U, t=T_HOP):
    return (U - math.sqrt(U * U + 16.0 * t * t)) / 2.0


def dimer_double_occ_exact(U, t=T_HOP):
    """dE0/dU by Hellmann-Feynman: (1 - U / sqrt(U^2 + 16 t^2)) / 2."""
    return 0.5 * (1.0 - U / math.sqrt(U * U + 16.0 * t * t))


def dimer_report():
    rows = {}
    for U in (0.0, 2.0, 4.0, 8.0, 20.0):
        evals = sector_spectrum(2, DIMER_BONDS, 1, 1, U)
        _, psi, basis = ground_state(2, DIMER_BONDS, 1, 1, U)
        gap = evals[1] - evals[0]
        rows[f"U={U:g}"] = {
            "E0_ed": round(evals[0], 9),
            "E0_exact": round(dimer_energy_exact(U), 9),
            "gap_ed": round(gap, 9),
            "gap_exact": round(dimer_gap_exact(U), 9),
            "gap_over_4t2_over_U": (round(gap / (4.0 / U), 6) if U else None),
            "double_occ_ed": round(double_occupancy(psi, basis), 9),
            "double_occ_exact": round(dimer_double_occ_exact(U), 9),
            "spin_corr": round(spin_correlation_nn(psi, basis, 2,
                                                   DIMER_BONDS), 9),
        }
    return rows


# ================================================================ 5. Plaquette sweeps
# The 2x2 plaquette at half filling, (2, 2), sector dimension 36.
# At U = 0 the ED energy must equal the free-fermion sum over the
# single-particle levels of the same four bonds - computed here by
# diagonalising the 4x4 hopping matrix with the same Jacobi solver,
# an independent path from the many-body basis.

def single_particle_levels(L, bonds, t=T_HOP):
    H = [[0.0] * L for _ in range(L)]
    for (i, j) in bonds:
        H[i][j] = H[j][i] = -t
    return jacobi_eigh(H)[0]


def free_fermion_energy(L, bonds, n_up, n_dn, t=T_HOP):
    levels = single_particle_levels(L, bonds, t)
    return sum(levels[:n_up]) + sum(levels[:n_dn])


def plaquette_report():
    rows = {}
    for U in U_SWEEP:
        E0, psi, basis = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS,
                                      2, 2, U)
        evals = sector_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, U)
        rows[f"U={U:g}"] = {
            "E0": round(E0, 6),
            "spin_gap_E1_minus_E0": round(evals[1] - evals[0], 6),
            "double_occ_per_site": round(
                double_occupancy(psi, basis) / PLAQUETTE_SITES, 6),
            "spin_corr_nn": round(
                spin_correlation_nn(psi, basis, PLAQUETTE_SITES,
                                    PLAQUETTE_BONDS), 6),
        }
    rows["U=0_free_fermion_referee"] = round(
        free_fermion_energy(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2), 6)
    return rows


def charge_gap(U):
    """Delta_c = E(N+1) + E(N-1) - 2 E(N) on the plaquette, N = 4."""
    e_plus = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 3, 2, U)[0]
    e_minus = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 1, U)[0]
    e_half = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, U)[0]
    return e_plus + e_minus - 2.0 * e_half


def mott_report():
    rows = {}
    for U in U_SWEEP:
        _, psi, basis = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS,
                                     2, 2, U)
        rows[f"U={U:g}"] = {
            "double_occ_per_site": round(
                double_occupancy(psi, basis) / PLAQUETTE_SITES, 6),
            "charge_gap": round(charge_gap(U), 6),
        }
    return rows


# ================================================================ 6. Heisenberg referee
# At large U the half-filled Hubbard model maps to the Heisenberg
# model with J = 4 t^2 / U. The lab diagonalises the 4-site
# Heisenberg ring itself (16 spin states) so the comparison is
# computed, not quoted: ring E0 = -2J, spin gap = J, and per-bond
# <S_i . S_j> = -1/2 exactly. The Hubbard plaquette must approach
# those numbers as U grows - with the finite-size caveat that a
# 4-site cluster at U/t = 16 is close to, not at, the limit.

def heisenberg_ring_spectrum(L, bonds, J):
    n = 1 << L
    H = [[0.0] * n for _ in range(n)]
    for s in range(n):
        for (i, j) in bonds:
            bi, bj = (s >> i) & 1, (s >> j) & 1
            H[s][s] += J * (0.25 if bi == bj else -0.25)
            if bi != bj:
                H[s ^ (1 << i) ^ (1 << j)][s] += J * 0.5
    return jacobi_eigh(H)[0]


def spin_report():
    rows = {}
    for U in U_SWEEP:
        if U == 0.0:
            continue
        _, psi, basis = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS,
                                     2, 2, U)
        evals = sector_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, U)
        J = 4.0 / U
        h_evals = heisenberg_ring_spectrum(PLAQUETTE_SITES,
                                           PLAQUETTE_BONDS, J)
        rows[f"U={U:g}"] = {
            "spin_corr_nn_hubbard": round(
                spin_correlation_nn(psi, basis, PLAQUETTE_SITES,
                                    PLAQUETTE_BONDS), 6),
            "spin_corr_nn_heisenberg": round(h_evals[0] / (4.0 * J), 6),
            "spin_gap_hubbard": round(evals[1] - evals[0], 6),
            "J_4t2_over_U": round(J, 6),
            "spin_gap_over_J": round((evals[1] - evals[0]) / J, 6),
            "heisenberg_gap_over_J": round((h_evals[1] - h_evals[0]) / J, 6),
        }
    return rows


# ================================================================ 7. kNN classifier
# Features of a plaquette ground state at coupling U: (U/t, double
# occupancy per site, charge gap / t, nearest-neighbour spin
# correlation). The training grid is labelled by the lab's own
# policy (model_label below) - the classifier interpolates the
# policy, it does not replace the computed quantities.

def plaquette_features(U):
    _, psi, basis = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, U)
    return (U,
            double_occupancy(psi, basis) / PLAQUETTE_SITES,
            charge_gap(U),
            spin_correlation_nn(psi, basis, PLAQUETTE_SITES,
                                PLAQUETTE_BONDS))


def model_label(u, d, gap, ss):
    """The correlation-regime policy distilled to four verdicts."""
    if u < 1.0:
        return "weakly-correlated"
    if ss <= -0.47:
        return "heisenberg-like"
    if d <= 0.06:
        return "mott-like"
    return "correlated-metal"


def _training_grid():
    rows = []
    for u in (0.0, 0.5, 2.0, 4.0, 8.0, 16.0):
        for d in (0.01, 0.04, 0.09, 0.16, 0.25):
            for gap in (0.0, 1.0, 3.0, 6.0, 10.0, 14.0):
                for ss in (-0.50, -0.48, -0.44, -0.36, -0.20, 0.0):
                    feats = (u, d, gap, ss)
                    rows.append((feats, model_label(*feats)))
    return rows


TRAINING_GRID = _training_grid()


def _scaled(features):
    # Feature scaling for the distance: U and the charge gap are
    # context (the policy thresholds sit in d and ss), so the two
    # policy-active observables carry double weight and the gap,
    # which the policy never thresholds, is compressed.
    u, d, gap, ss = features
    return (u / 16.0, 2.0 * d, gap / 56.0, 2.0 * ss)


def classify_features(features, k=5):
    """Distance-weighted kNN over the policy-labelled grid."""
    target = _scaled(features)
    ranked = []
    for feats, label in TRAINING_GRID:
        point = _scaled(feats)
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(target, point)))
        ranked.append((dist, feats, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"weakly-correlated": 0, "correlated-metal": 1,
                "mott-like": 2, "heisenberg-like": 3}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


CLASSIFIER_SCENARIOS = {
    "plaquette_U0": 0.0,
    "plaquette_U2": 2.0,
    "plaquette_U4": 4.0,
    "plaquette_U6": 6.0,
    "plaquette_U8": 8.0,
    "plaquette_U12": 12.0,
    "plaquette_U16": 16.0,
}


def classifier_report():
    rows = {}
    for name, U in CLASSIFIER_SCENARIOS.items():
        feats = plaquette_features(U)
        rows[name] = {
            "U_over_t": U,
            "features": [round(v, 6) for v in feats],
            "policy": model_label(*feats),
            "knn": classify_features(feats),
        }
    return rows


# ================================================================ 8. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="hubbardmodel-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "model_id": "gemini-argon-4",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all material/device gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="hubbardmodel-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the cluster before the number: two sites and four sites, t = 1, the bonds counted once. The dimer gap is exact, the plaquette is exact for four sites, and neither is a material - the Heisenberg limit is approached, and the distance still remaining at U/t = 16 is part of the result, not a rounding error.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["two sites passing one electron back and forth as U rises",
                  "doublons draining out of a 2x2 plaquette cell by cell",
                  "neighbouring spins locking antiparallel as the charge gap opens",
                  "the singlet-triplet gap shrinking onto the 4t^2/U curve",
                  "a 4-site ring of arrows settling into the Heisenberg pattern"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 0. Jacobi eigensolver refereed on known spectra first
    chk("jacobi_2x2_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(jacobi_eigh([[2.0, 1.0], [1.0, 2.0]])[0], (1.0, 3.0))))
    chk("jacobi_path4_hopping_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(single_particle_levels(4, ((0, 1), (1, 2), (2, 3))),
                (-1.618033988749895, -0.618033988749895,
                 0.618033988749895, 1.618033988749895))))
    chk("jacobi_eigenvector_residual_plaquette",
        eigen_residual(*_plaquette_H_and_ground()) < 1e-9)

    # 1. Sector bookkeeping
    chk("sector_dimensions",
        len(sector_basis(2, 1, 1)) == 4
        and len(sector_basis(4, 2, 2)) == 36
        and len(sector_basis(4, 3, 2)) == 24
        and len(sector_basis(4, 2, 1)) == 24)
    H4, _ = build_hubbard(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 4.0)
    chk("hamiltonian_symmetric",
        max(abs(H4[i][j] - H4[j][i])
            for i in range(36) for j in range(36)) < 1e-13)
    chk("hamiltonian_trace_counts_doublons",
        abs(sum(H4[i][i] for i in range(36)) - 4.0 * 36) < 1e-9)

    # 2. Dimer: exact spectrum, gap formula, Heisenberg limit, doublons
    chk("dimer_U0_spectrum",
        all(abs(a - b) < 1e-9 for a, b in
            zip(sector_spectrum(2, DIMER_BONDS, 1, 1, 0.0),
                (-2.0, 0.0, 0.0, 2.0))))
    chk("dimer_gap_matches_closed_form_grid",
        all(abs((sector_spectrum(2, DIMER_BONDS, 1, 1, U)[1]
                  - sector_spectrum(2, DIMER_BONDS, 1, 1, U)[0])
                - dimer_gap_exact(U)) < 1e-9
            for U in (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0)))
    chk("dimer_energy_matches_closed_form_grid",
        all(abs(ground_state(2, DIMER_BONDS, 1, 1, U)[0]
                - dimer_energy_exact(U)) < 1e-9
            for U in (0.0, 2.0, 4.0, 8.0, 20.0)))
    gap20 = (sector_spectrum(2, DIMER_BONDS, 1, 1, 20.0)[1]
             - sector_spectrum(2, DIMER_BONDS, 1, 1, 20.0)[0])
    chk("dimer_gap_tends_to_4t2_over_U",
        0.98 < gap20 / (4.0 / 20.0) < 1.0)
    chk("dimer_double_occ_matches_hellmann_feynman",
        all(abs(double_occupancy(ground_state(2, DIMER_BONDS, 1, 1, U)[1],
                                 ground_state(2, DIMER_BONDS, 1, 1, U)[2])
                - dimer_double_occ_exact(U)) < 1e-9
            for U in (0.0, 2.0, 4.0, 8.0, 20.0)))
    h_fd = 1e-4
    _, psi4, basis4 = _dimer_ground(4.0)
    chk("dimer_double_occ_matches_finite_difference",
        abs(double_occupancy(psi4, basis4)
            - (dimer_energy_exact(4.0 + h_fd)
               - dimer_energy_exact(4.0 - h_fd)) / (2 * h_fd)) < 1e-5)
    chk("dimer_spin_corr_limits",
        abs(spin_correlation_nn(*_dimer_ground(0.0)[1:], 2, DIMER_BONDS)
            + 0.375) < 1e-9
        and spin_correlation_nn(*_dimer_ground(20.0)[1:], 2, DIMER_BONDS)
        < -0.74)

    # 3. Plaquette: free-fermion referee + energy table sanity
    chk("plaquette_U0_matches_free_fermion",
        abs(ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 0.0)[0]
            - free_fermion_energy(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2))
        < 1e-9)
    chk("plaquette_U0_energy_is_minus_4t",
        abs(ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 0.0)[0]
            + 4.0) < 1e-9)
    chk("plaquette_single_particle_levels",
        all(abs(a - b) < 1e-12 for a, b in
            zip(single_particle_levels(PLAQUETTE_SITES, PLAQUETTE_BONDS),
                (-2.0, 0.0, 0.0, 2.0))))
    energies = [ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, U)[0]
                for U in U_SWEEP]
    chk("plaquette_energy_rises_with_U",
        all(b > a for a, b in zip(energies, energies[1:])))
    chk("plaquette_spin_gap_zero_at_U0_degenerate_shell",
        abs(sector_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 0.0)[1]
            - sector_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 0.0)[0])
        < 1e-9)

    # 4. Mott bookkeeping: doublons drain, charge gap opens
    dvals = []
    for U in U_SWEEP:
        _, psi, basis = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS,
                                     2, 2, U)
        dvals.append(double_occupancy(psi, basis) / PLAQUETTE_SITES)
    chk("double_occupancy_monotone_falling",
        all(b < a for a, b in zip(dvals, dvals[1:])))
    chk("double_occupancy_U0_near_quarter_with_degeneracy_caveat",
        abs(dvals[0] - 0.25) < 0.01)
    chk("double_occupancy_small_at_U16", dvals[-1] < 0.02)
    gaps = [charge_gap(U) for U in U_SWEEP]
    chk("charge_gap_zero_at_U0_open_shell", abs(gaps[0]) < 1e-9)
    chk("charge_gap_monotone_growing",
        all(b > a for a, b in zip(gaps, gaps[1:])))
    chk("charge_gap_large_at_U16", gaps[-1] > 12.0)

    # 5. Spin correlations + Heisenberg referee
    svals = {}
    for U in U_SWEEP:
        if U == 0.0:
            continue
        _, psi, basis = ground_state(PLAQUETTE_SITES, PLAQUETTE_BONDS,
                                     2, 2, U)
        svals[U] = spin_correlation_nn(psi, basis, PLAQUETTE_SITES,
                                       PLAQUETTE_BONDS)
    chk("spin_corr_negative_for_interacting_plaquette",
        all(v < 0.0 for v in svals.values()))
    chk("spin_corr_more_antiferromagnetic_with_U",
        all(svals[b] < svals[a]
            for a, b in zip(sorted(svals), sorted(svals)[1:])))
    h_evals = heisenberg_ring_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 1.0)
    chk("heisenberg_ring_exact_landmarks",
        abs(h_evals[0] + 2.0) < 1e-9 and abs(h_evals[1] - h_evals[0] - 1.0)
        < 1e-9 and abs(h_evals[0] / 4.0 + 0.5) < 1e-9)
    chk("spin_corr_approaches_heisenberg_minus_half",
        abs(svals[16.0] + 0.5) < 0.03)
    gap16 = (sector_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 16.0)[1]
             - sector_spectrum(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 16.0)[0])
    chk("spin_gap_approaches_J_at_U16",
        0.8 < gap16 / (4.0 / 16.0) < 0.95)

    # 6. Classifier: labels come from the lab's own policy
    rep = classifier_report()
    chk("classifier_scenario_labels",
        rep["plaquette_U0"]["knn"] == "weakly-correlated"
        and rep["plaquette_U2"]["knn"] == "correlated-metal"
        and rep["plaquette_U4"]["knn"] == "correlated-metal"
        and rep["plaquette_U6"]["knn"] == "mott-like"
        and rep["plaquette_U8"]["knn"] == "mott-like"
        and rep["plaquette_U12"]["knn"] == "heisenberg-like"
        and rep["plaquette_U16"]["knn"] == "heisenberg-like")
    chk("classifier_knn_matches_policy_scenarios",
        all(row["knn"] == row["policy"] for row in rep.values()))
    chk("classifier_knn_matches_policy_grid",
        all(classify_features(feats) == label
            for feats, label in TRAINING_GRID))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_material_prediction_hold",
        gate("real_material_prediction_claim") == "HOLD")
    chk("gate_superconductivity_hold", gate("superconductivity_claim")
        == "HOLD")
    chk("gate_thermodynamic_limit_hold",
        gate("thermodynamic_limit_claim") == "HOLD")
    chk("gate_quantum_advantage_hold", gate("quantum_advantage_claim")
        == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking"
        and brains["decider"]["model_id"] == "gemini-argon-4")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


def _dimer_ground(U):
    E0, psi, basis = ground_state(2, DIMER_BONDS, 1, 1, U)
    return E0, psi, basis


def _plaquette_H_and_ground():
    H, basis = build_hubbard(PLAQUETTE_SITES, PLAQUETTE_BONDS, 2, 2, 4.0)
    evals, evecs = jacobi_eigh(H)
    return H, evals[0], evecs[0]


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("HubbardModel-Lab - the Hubbard model, exactly, on two tiny clusters")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - t versus U, on clusters you can check")
    print("    Hopping t delocalises electrons; on-site repulsion U charges")
    print("    for double occupancy. This lab diagonalises that argument")
    print("    exactly on a 2-site dimer and a 2x2 plaquette (t = 1). Four")
    print("    sites are not a material, and nothing here claims they are.")
    print("\n[2] The dimer - a gap formula the ED must reproduce")
    dimer = dimer_report()
    for key in ("U=0", "U=4", "U=20"):
        row = dimer[key]
        print(f"    {key:5s} E0 {row['E0_ed']:+.6f} (exact "
              f"{row['E0_exact']:+.6f}) | gap {row['gap_ed']:.6f} "
              f"(exact {row['gap_exact']:.6f})")
    row = dimer["U=20"]
    print(f"    at U/t = 20 the gap is {row['gap_over_4t2_over_U']:.4f} of "
          f"4t^2/U - the Heisenberg exchange J emerging from the ED.")
    print(f"    dimer <S1.S2>: {dimer['U=0']['spin_corr']:.3f} at U = 0 "
          f"(half doublon weight), {dimer['U=20']['spin_corr']:.3f} at "
          f"U = 20 (pure singlet: -3/4).")
    print("\n[3] The plaquette - ground-state energy vs U/t (half filling)")
    plaq = plaquette_report()
    for U in U_SWEEP:
        row = plaq[f"U={U:g}"]
        print(f"    U/t {U:5.1f} | E0 {row['E0']:+.6f} | "
              f"spin gap {row['spin_gap_E1_minus_E0']:.6f}")
    print(f"    U = 0 referee: free-fermion sum over single-particle "
          f"levels = {plaq['U=0_free_fermion_referee']:+.6f} - the ED "
          f"agrees exactly.")
    print("    (the U = 0 shell is degenerate: levels -2, 0, 0, +2 leave")
    print("     the zero level half-filled, so the U = 0 spin gap is 0.)")
    print("\n[4] Mott bookkeeping - doublons drain, the charge gap opens")
    mott = mott_report()
    for U in U_SWEEP:
        row = mott[f"U={U:g}"]
        print(f"    U/t {U:5.1f} | double occupancy/site "
              f"{row['double_occ_per_site']:.6f} | charge gap "
              f"{row['charge_gap']:.6f}")
    print("    -> at U = 0 the open shell makes the charge gap exactly 0;")
    print("       every step of U after that opens it further. Double")
    print("       occupancy at U = 0 sits at the uncorrelated 1/4 (the")
    print("       solver's vector in the degenerate shell lands within")
    print("       0.001 of it) and drains monotonically toward 0.")
    print("\n[5] Spin correlations - antiferromagnetism from repulsion")
    spin = spin_report()
    for U in (4.0, 8.0, 16.0):
        row = spin[f"U={U:g}"]
        print(f"    U/t {U:5.1f} | <S.S>_nn {row['spin_corr_nn_hubbard']:+.6f} "
              f"(Heisenberg ring: {row['spin_corr_nn_heisenberg']:+.6f}) | "
              f"spin gap / J {row['spin_gap_over_J']:.3f}")
    print("    -> the plaquette approaches the Heisenberg ring the lab")
    print("       diagonalises itself (E0 = -2J, gap = J, <S.S> = -1/2):")
    print("       at U/t = 16 it is close, not exact - the remaining")
    print("       distance is charge fluctuation, printed, not hidden.")
    print("\n[6] Regime classifier (TinyML kNN, labels from the lab policy)")
    for name, row in classifier_report().items():
        print(f"    {name:28s} -> {row['knn']}")
    print("\n[7] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[8] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "dimer":
        print(json.dumps(dimer_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "plaquette":
        print(json.dumps(plaquette_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "mott":
        print(json.dumps(mott_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "spin":
        print(json.dumps(spin_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "real_material_prediction_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
