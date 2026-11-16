"""Small deterministic file and provenance helpers."""
from yash_vaccine.paths import resolve_input
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path

def require_assertions(context):
    """Refuse assertion-based verification when Python has removed its checks."""
    if not __debug__:
        raise ValueError(context+' uses integrity checks that -O/PYTHONOPTIMIZE disables. '
                         'Rerun without -O/PYTHONOPTIMIZE so the evidence is checked.')

def child_env():
    """Run recorded verification workers with their assertions enabled."""
    return {name:value for name,value in os.environ.items() if name!='PYTHONOPTIMIZE'}

def open_text(path):
    path=Path(path)
    return gzip.open(resolve_input(path),'rt') if path.suffix=='.gz' else resolve_input(path).open()

def rows(path):
    with open_text(path) as f:yield from csv.DictReader(f,delimiter='\t')

def write_tsv(path,records):
    records=list(records);path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for row in records for k in row))
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(records)

def write_json(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,default=str)+'\n')

def sha256(path):
    h=hashlib.sha256()
    with resolve_input(Path(path)).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def ranked(records):
    return sorted(records,key=lambda r:(r['weighted_score_left']=='',-float(r['weighted_score_left'] or 0),r['event_id'],r['region_id']))
