"""Read-only Mac memory-pressure and swap follow-up; no pressure generation."""
from pathlib import Path
import datetime,json,re,subprocess,time

def read():
    raw={k:subprocess.check_output(cmd,text=True) for k,cmd in {
        'vm_stat':['vm_stat'],'memory_pressure':['memory_pressure','-Q'],
        'swap_usage':['sysctl','-n','vm.swapusage']}.items()}
    counts={k.strip():int(v) for k,v in re.findall(r'^([^\n:]+):\s+(\d+)\.',raw['vm_stat'],re.M)}
    return dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),**raw,counters=counts,
        compressed_physical_gib=counts['Pages occupied by compressor']*16384/2**30,
        compressed_logical_gib=counts['Pages stored in compressor']*16384/2**30)

if __name__=='__main__':
    a=read();time.sleep(5);b=read()
    record=dict(first=a,second=b,swapins_per_second=(b['counters']['Swapins']-a['counters']['Swapins'])/5,
        swapouts_per_second=(b['counters']['Swapouts']-a['counters']['Swapouts'])/5,concurrency_unchanged=True)
    Path('results/julia_checkpoint_ram_release_memory_followup.json').write_text(json.dumps(record,indent=2)+'\n')
    print('Compressed physical GiB:',b['compressed_physical_gib'],'stored logical GiB:',b['compressed_logical_gib'])
    print(b['swap_usage'].strip());print(b['memory_pressure'].splitlines()[-1])
    print('Swap traffic pages/s:',record['swapins_per_second'],record['swapouts_per_second'])
