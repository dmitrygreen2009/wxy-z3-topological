"""Stationarity and known-sector energy gates, separate from phase inference."""
import math

def energy_quality_gates(record, *, variance_tolerance=1e-8, known_energy_tolerance=1e-9):
    variance=record.get('energy_variance')
    variance_gate=None if variance is None else math.isfinite(variance) and abs(variance)<variance_tolerance
    target=record.get('known_fixed_sector_ed_energy')
    error=None if target is None else record['energy']-target
    known_gate=None if error is None else math.isfinite(error) and abs(error)<known_energy_tolerance
    return dict(energy_variance=variance_gate,known_sector_energy=known_gate)
