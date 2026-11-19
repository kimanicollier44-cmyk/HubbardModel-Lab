"""Tests for HubbardModel-Lab. Every test grades a computed number
against a closed form or an independent implementation in main.py."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as hm

TOL = 1e-9


def test_jacobi_known_spectra():
    evals, _ = hm.jacobi_eigh([[2.0, 1.0], [1.0, 2.0]])
    assert all(abs(a - b) < 1e-12 for a, b in zip(evals, (1.0, 3.0)))
    levels = hm.single_particle_levels(4, ((0, 1), (1, 2), (2, 3)))
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(levels, (-1.618033988749895, -0.618033988749895,
                            0.618033988749895, 1.618033988749895)))


def test_eigenvector_residual():
    H, basis = hm.build_hubbard(4, hm.PLAQUETTE_BONDS, 2, 2, 4.0)
    evals, evecs = hm.jacobi_eigh(H)
    assert hm.eigen_residual(H, evals[0], evecs[0]) < 1e-9


def test_sector_dimensions():
    assert len(hm.sector_basis(2, 1, 1)) == 4
    assert len(hm.sector_basis(4, 2, 2)) == 36
    assert len(hm.sector_basis(4, 3, 2)) == 24
    assert len(hm.sector_basis(4, 2, 1)) == 24


def test_hamiltonian_symmetric_and_trace():
    H, _ = hm.build_hubbard(4, hm.PLAQUETTE_BONDS, 2, 2, 4.0)
    assert max(abs(H[i][j] - H[j][i])
               for i in range(36) for j in range(36)) < 1e-13
    # mean doublon count over the (2,2) sector is exactly 1 per state
    assert abs(sum(H[i][i] for i in range(36)) - 4.0 * 36) < 1e-9


def test_dimer_U0_spectrum():
    evals = hm.sector_spectrum(2, hm.DIMER_BONDS, 1, 1, 0.0)
    assert all(abs(a - b) < TOL for a, b in zip(evals, (-2.0, 0.0, 0.0, 2.0)))


def test_dimer_gap_closed_form():
    for U in (0.0, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0):
        evals = hm.sector_spectrum(2, hm.DIMER_BONDS, 1, 1, U)
        assert abs((evals[1] - evals[0]) - hm.dimer_gap_exact(U)) < TOL
        assert abs(evals[0] - hm.dimer_energy_exact(U)) < TOL


def test_dimer_gap_heisenberg_limit():
    evals = hm.sector_spectrum(2, hm.DIMER_BONDS, 1, 1, 20.0)
    ratio = (evals[1] - evals[0]) / (4.0 / 20.0)
    assert 0.98 < ratio < 1.0


def test_dimer_double_occupancy():
    for U in (0.0, 2.0, 4.0, 8.0, 20.0):
        _, psi, basis = hm.ground_state(2, hm.DIMER_BONDS, 1, 1, U)
        assert abs(hm.double_occupancy(psi, basis)
                   - hm.dimer_double_occ_exact(U)) < TOL


def test_dimer_spin_correlation():
    _, psi, basis = hm.ground_state(2, hm.DIMER_BONDS, 1, 1, 0.0)
    assert abs(hm.spin_correlation_nn(psi, basis, 2, hm.DIMER_BONDS)
               + 0.375) < TOL
    _, psi, basis = hm.ground_state(2, hm.DIMER_BONDS, 1, 1, 20.0)
    assert hm.spin_correlation_nn(psi, basis, 2, hm.DIMER_BONDS) < -0.74


def test_plaquette_free_fermion_referee():
    E0 = hm.ground_state(4, hm.PLAQUETTE_BONDS, 2, 2, 0.0)[0]
    assert abs(E0 - hm.free_fermion_energy(4, hm.PLAQUETTE_BONDS, 2, 2)) < TOL
    assert abs(E0 + 4.0) < TOL
    levels = hm.single_particle_levels(4, hm.PLAQUETTE_BONDS)
    assert all(abs(a - b) < 1e-12 for a, b in zip(levels, (-2.0, 0.0, 0.0, 2.0)))


def test_plaquette_energy_rises_with_U():
    energies = [hm.ground_state(4, hm.PLAQUETTE_BONDS, 2, 2, U)[0]
                for U in hm.U_SWEEP]
    assert all(b > a for a, b in zip(energies, energies[1:]))


def test_plaquette_degenerate_shell_at_U0():
    evals = hm.sector_spectrum(4, hm.PLAQUETTE_BONDS, 2, 2, 0.0)
    assert abs(evals[1] - evals[0]) < TOL


def test_double_occupancy_monotone():
    vals = []
    for U in hm.U_SWEEP:
        _, psi, basis = hm.ground_state(4, hm.PLAQUETTE_BONDS, 2, 2, U)
        vals.append(hm.double_occupancy(psi, basis) / 4)
    assert all(b < a for a, b in zip(vals, vals[1:]))
    assert abs(vals[0] - 0.25) < 0.01  # degenerate shell: solver's vector
    assert vals[-1] < 0.02


def test_charge_gap_behaviour():
    gaps = [hm.charge_gap(U) for U in hm.U_SWEEP]
    assert abs(gaps[0]) < TOL  # open shell at U = 0: gap exactly 0
    assert all(b > a for a, b in zip(gaps, gaps[1:]))
    assert gaps[-1] > 12.0


def test_spin_correlation_trend():
    vals = {}
    for U in hm.U_SWEEP:
        if U == 0.0:
            continue
        _, psi, basis = hm.ground_state(4, hm.PLAQUETTE_BONDS, 2, 2, U)
        vals[U] = hm.spin_correlation_nn(psi, basis, 4, hm.PLAQUETTE_BONDS)
    assert all(v < 0.0 for v in vals.values())
    keys = sorted(vals)
    assert all(vals[b] < vals[a] for a, b in zip(keys, keys[1:]))
    assert abs(vals[16.0] + 0.5) < 0.03


def test_heisenberg_ring_landmarks():
    evals = hm.heisenberg_ring_spectrum(4, hm.PLAQUETTE_BONDS, 1.0)
    assert abs(evals[0] + 2.0) < TOL
    assert abs(evals[1] - evals[0] - 1.0) < TOL
    assert abs(evals[0] / 4.0 + 0.5) < TOL


def test_spin_gap_approaches_J():
    evals = hm.sector_spectrum(4, hm.PLAQUETTE_BONDS, 2, 2, 16.0)
    ratio = (evals[1] - evals[0]) / (4.0 / 16.0)
    assert 0.8 < ratio < 0.95


def test_classifier_scenarios():
    rep = hm.classifier_report()
    assert rep["plaquette_U0"]["knn"] == "weakly-correlated"
    assert rep["plaquette_U2"]["knn"] == "correlated-metal"
    assert rep["plaquette_U4"]["knn"] == "correlated-metal"
    assert rep["plaquette_U6"]["knn"] == "mott-like"
    assert rep["plaquette_U8"]["knn"] == "mott-like"
    assert rep["plaquette_U12"]["knn"] == "heisenberg-like"
    assert rep["plaquette_U16"]["knn"] == "heisenberg-like"
    assert all(row["knn"] == row["policy"] for row in rep.values())


def test_classifier_grid_consistency():
    assert all(hm.classify_features(feats) == label
               for feats, label in hm.TRAINING_GRID)


def test_hold_gates():
    assert len(hm.HOLD_GATES) == 14
    assert hm.gate("real_material_prediction_claim") == "HOLD"
    assert hm.gate("superconductivity_claim") == "HOLD"
    assert hm.gate("thermodynamic_limit_claim") == "HOLD"
    assert hm.gate("run_simulation") == "PROCEED"


def test_brains_seats():
    brains = hm.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["model_id"] == "gemini-argon-4"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = hm.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_hold_demo_exits_2():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "hold-demo"],
        capture_output=True, text=True)
    assert proc.returncode == 2
    assert "HOLD" in proc.stdout
