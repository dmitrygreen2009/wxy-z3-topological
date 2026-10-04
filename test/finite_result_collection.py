"""Regression for the omitted lower-energy number-sector scan candidates."""
import json,pathlib,sys,tempfile,unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/'scripts'))
from finite_results import finite_groups
class CollectionTest(unittest.TestCase):
 def test_sector_envelope_does_not_hide_candidates(self):
  with tempfile.TemporaryDirectory() as tmp:
   root=pathlib.Path(tmp)
   ordinary={'family':'armchair','length':4,'width':2,'spins':76,'nup':38,'records':[{'cap':256,'maxlinkdim':256,'energy':-30.5,'entropy':2.7}]}
   candidate={'family':'armchair','length':4,'width':2,'physical_spins':76,'nup':33,'cap':256,'bond_dimension':256,'energy':-30.6,'entropy':2.6,'seed':7198}
   (root/'ordinary.json').write_text(json.dumps(ordinary))
   (root/'scan.json').write_text(json.dumps({'geometry':{'family':'armchair'},'records':[candidate]}))
   mirror=dict(candidate,nup=43)
   (root/'mirror.json').write_text(json.dumps({'source_scan':'scan.json','records':[mirror]}))
   rows=list(finite_groups(str(root/'*.json')))
   self.assertEqual({r['nup'] for p,r,s in rows},{33,38,43})
   best=min(rows,key=lambda row:row[2][-1][1]['energy'])
   self.assertAlmostEqual(best[2][-1][1]['energy'],-30.6)
   self.assertEqual(best[2][-1][1]['maxlinkdim'],256)
   self.assertIn(best[1]['nup'],[33,43])
if __name__=='__main__':unittest.main()
