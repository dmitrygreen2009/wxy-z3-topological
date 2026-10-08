"""Inexpensive regression of physical identities and salvaged result consistency."""
import json,pathlib,unittest,math
ROOT=pathlib.Path(__file__).resolve().parents[1]
def read(p):return json.loads((ROOT/p).read_text())
class SalvagePhysics(unittest.TestCase):
 def test_trial_vs_saved_exact_ground_energies(self):
  b=read('results/physics48_density_bound.json');exact=read('results/ed.json')
  for r in b['fixed_number_benchmark_expectations']:
   self.assertGreaterEqual(r['energy']+1e-11,min(exact[r['geometry']]['energies']))
   self.assertLess(r['error'],1e-12)
  self.assertAlmostEqual(b['analytic_density_window'][0],1/(3*math.sqrt(2)),14)
  self.assertLess(b['analytic_density_window'][0],b['cluster_density_window'][0])
 def test_operator_charge_conservation_and_covariance(self):
  r=read('results/physics48_strings.json')
  for x in r['open_CGS_unitaries']:
   self.assertEqual(len(x['endpoint_vertices_zero_based']),2)
   self.assertLess(x['interior_H_invariance_error'],1e-12)
   self.assertLess(max(x['CGS_commutator_errors'].values()),1e-12)
   self.assertLess(max(x['energy_costs_four_ground_states']),x['endpoint_operator_norm_bound'])
  paths=read('results/physics48_string_paths.json');self.assertLess(paths['ground_vector_relation_error'],1e-12);self.assertLess(max(paths['ground_charge_phase_errors']),1e-9)
  q=r['neutral_ladder_pair'];self.assertEqual(q['number_change'],0)
  self.assertEqual(q['final_charges'],[2,0]);self.assertLess(max(q['covariance_errors'].values()),1e-12)
  self.assertGreater(q['energy'],q['exact_sector_minimum_reused'])
 def test_tiling_trial_and_independent_cut_fixture(self):
  trial=read('results/physics48_tiled_trial.json');b=read('results/physics48_tiled_density_bound.json')
  self.assertEqual(trial['patch_physical_spin_count'],9*trial['primitive_cells_per_patch'])
  self.assertEqual(len(set(trial['kept_physical_spin_ids'])),trial['patch_physical_spin_count'])
  self.assertFalse(trial['optimization_performed']);self.assertLess(trial['bulk_trial_energy_per_primitive_cell'],-3)
  self.assertGreater(b['upward_rounded_trial_upper'],trial['bulk_trial_energy_per_primitive_cell'])
  pair=read('results/physics48_pair_density_bound.json');self.assertEqual(sum(r['dimension'] for r in pair['records']),2048);self.assertGreaterEqual(pair['density_window'][0],b['same_local_cluster_density_window'][0]-1e-9)
  for r in read('results/physics48_tiled_trial_validation.json')['records']:
   self.assertAlmostEqual(r['full_H_energy'],-1/math.sqrt(3),12)
   self.assertAlmostEqual(r['retained_H_energy'],r['expected_retained_energy'],12)
 def test_energy_inventory_ph_and_no_overstatement(self):
  r=read('results/physics48_filling_energies.json')
  for x in r['exact_particle_hole_checks']:self.assertLess(x['particle_hole_energy_error'],1e-8)
  for x in r['records']:
   self.assertAlmostEqual(x['energy_per_site'],x['energy']/x['spins'],13)
   self.assertAlmostEqual(x['filling'],x['nup']/x['spins'],13)
  ed={x['geometry']:x for x in r['minima'] if x['method']=='ED'}
  self.assertEqual(ed['armchair_L2_w1']['lowest_recorded_numbers'],[8,12])
  self.assertEqual(ed['four_star_torus']['lowest_recorded_numbers'],[9])
  self.assertEqual(read('results/physics48_summary.json')['accepted_gamma_fit_points'],0)
  r=read('results/physics48_infinite_energy.json');self.assertLess(max(r['state_consistency_errors'].values()),r['unchanged_consistency_tolerance']);self.assertTrue(r['excluded_as_unrestricted_global_ground_candidate']);self.assertGreater(r['energy_above_trial_per_vertex'],.06)
if __name__=='__main__':unittest.main()
