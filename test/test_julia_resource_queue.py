"""Scheduler regression tests: synthetic processes only; no Julia jobs launched."""
import importlib.util,json,pathlib,tempfile,unittest
from unittest.mock import patch,Mock
spec=importlib.util.spec_from_file_location('queue',pathlib.Path(__file__).parents[1]/'scripts/julia_resource_queue.py')
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)

class QueueTests(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.path=pathlib.Path(self.tmp.name)/'queue.json'
  self.proc=dict(pid=123,ppid=12,state='T',started='date',command='/julia --project=. scripts/extend_study.jl')
  self.job=dict(process=self.proc.copy(),purpose='protected chi512',protected=True)
  self.data=dict(julia_executable='/julia',critical_process=dict(self.proc,pid=999),resume_queue=[self.job])
 def write(self):q.save(self.path,self.data)
 def read(self):return json.loads(self.path.read_text())
 def gates(self):
  return patch.object(q.pathlib.Path,'exists',return_value=True),patch.object(q.pathlib.Path,'read_text',side_effect=lambda p: '')
 def run_tick(self,procs):
  original=q.pathlib.Path.read_text
  def read(p,*a,**kw):
   if 'winding' in str(p):return json.dumps({'records':[{'converged_within_fixed_number_and_loop_sector':True}]})
   return original(p,*a,**kw)
  with patch.object(q,'processes',return_value=procs),patch.object(q.pathlib.Path,'exists',return_value=True),patch.object(q.pathlib.Path,'read_text',read),patch.object(q.subprocess,'run'):
   q.tick(self.path,{})
 def test_protected_exact_process_resumes(self):
  self.write()
  with patch.object(q.os,'kill') as kill,patch.object(q.subprocess,'Popen') as launch:
   self.run_tick({123:self.proc});kill.assert_called_once_with(123,q.signal.SIGCONT);launch.assert_not_called()
  self.assertEqual(self.read()['queue_state'],'RUNNING')
 def test_surviving_active_prevents_duplicate(self):
  self.write()
  with patch.object(q.os,'kill') as kill,patch.object(q.subprocess,'Popen') as launch:
   self.run_tick({123:dict(self.proc,state='R')});kill.assert_not_called();launch.assert_not_called()
 def test_missing_protected_blocks_not_restart(self):
  self.write();self.run_tick({})
  self.assertEqual(self.read()['queue_state'],'BLOCKED');self.assertTrue((self.path.parent/'QUEUE_BLOCKED.txt').is_file())
 def test_exited_unverified_is_not_complete(self):
  self.data['resume_queue']=[dict(self.job,protected=False)];self.write();self.run_tick({})
  self.assertEqual(self.read()['queue_state'],'BLOCKED')
  self.assertIsNone(self.read()['resume_queue'][0]['exit_code'])
 def test_verified_completion_record(self):
  self.data['resume_queue']=[dict(self.job,execution_finished=True,scientific_completion_verified=True)];self.write();self.run_tick({})
  self.assertEqual(self.read()['queue_state'],'COMPLETE');self.assertTrue((self.path.parent/'QUEUE_COMPLETE.txt').is_file())
 def test_blocked_does_not_repeat_launch(self):
  self.data['queue_state']='BLOCKED';self.write()
  with patch.object(q.subprocess,'Popen') as launch,patch.object(q.os,'kill') as kill:
   self.run_tick({123:self.proc});launch.assert_not_called();kill.assert_not_called()
 def test_validated_next_checkpoint_launches_once(self):
  snapshot=self.path.parent/'saved.jls';snapshot.write_bytes(b'fixture')
  import hashlib
  self.data['resume_queue']=[dict(self.job,protected=False,restart_ready=True,restart_released=True,
   restart_validation={'validated':True},restart_log=str(self.path.parent/'solver.log'),
   restart_plan={'snapshot':str(snapshot),'snapshot_sha256':hashlib.sha256(b'fixture').hexdigest(),'target_cap':64})]
  self.write();child=Mock();child.pid=456
  with patch.object(q.subprocess,'Popen',return_value=child) as launch:
   self.run_tick({});launch.assert_called_once()
  self.assertEqual(self.read()['resume_queue'][0]['restart_pid'],456)
 def test_child_exit_status_recorded_and_next_resumed(self):
  old=dict(self.job,process=dict(self.proc,pid=111),protected=False,restart_launched=True,restart_pid=222)
  self.data['resume_queue']=[old,self.job];self.write();child=Mock();child.poll.return_value=7;child.returncode=7
  original=q.pathlib.Path.read_text
  def read(p,*a,**kw):
   if 'winding' in str(p):return json.dumps({'records':[{'converged_within_fixed_number_and_loop_sector':True}]})
   return original(p,*a,**kw)
  with patch.object(q,'processes',return_value={123:self.proc}),patch.object(q,'checkpoint_record',return_value={'log':'test'}),patch.object(q.pathlib.Path,'exists',return_value=True),patch.object(q.pathlib.Path,'read_text',read),patch.object(q.os,'kill') as kill:
   q.tick(self.path,{'111':child});kill.assert_called_once_with(123,q.signal.SIGCONT)
  self.assertEqual(self.read()['resume_queue'][0]['exit_code'],7)

if __name__=='__main__':unittest.main()
