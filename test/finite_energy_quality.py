"""A stable but nonstationary sector plateau must never pass convergence."""
import json,pathlib,sys,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from finite_energy_quality import energy_quality_gates

class EnergyQualityTest(unittest.TestCase):
    def test_preserved_cold_plateau_is_rejected(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        record=json.loads((root/'results/archive/winding_optimizer_plateau/zigzag_L2_w1_Nup9_winding0_cold_plateau.json').read_text())['records'][-1]
        self.assertLess(max(record['sweep_max_truncation_errors'][-4:]),1e-10)
        self.assertLess(abs(record['sweep_energies'][-1]-record['sweep_energies'][-2]),1e-9)
        self.assertFalse(energy_quality_gates(record)['energy_variance'])

    def test_eigenstate_quality_does_not_override_wrong_ground_energy(self):
        gates=energy_quality_gates(dict(energy_variance=0.,energy=-1.,known_fixed_sector_ed_energy=-2.))
        self.assertTrue(gates['energy_variance']);self.assertFalse(gates['known_sector_energy'])

    def test_missing_or_nonfinite_evidence_cannot_pass(self):
        self.assertIsNone(energy_quality_gates({})['energy_variance'])
        self.assertFalse(energy_quality_gates(dict(energy_variance=float('nan')))['energy_variance'])

    def test_validated_small_state_retains_energy_quality(self):
        root=pathlib.Path(__file__).resolve().parents[1]
        record=json.loads((root/'results/armchair_L2_w1_chi256_star_Nup8_variance.json').read_text())
        record['known_fixed_sector_ed_energy']=-8.037741645057547
        self.assertTrue(all(energy_quality_gates(record).values()))

if __name__=='__main__':unittest.main()
