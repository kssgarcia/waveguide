"""Execution bookkeeping only: no numerical criteria or Green modifications."""
from __future__ import annotations
import hashlib
import importlib
import json
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import redirect_stdout, redirect_stderr
from dataclasses import asdict, fields, is_dataclass, replace
from datetime import datetime
from pathlib import Path
from uuid import uuid4
from zoneinfo import ZoneInfo
import numpy as np

HARNESS_VERSION = 'theorem-2.3-execution-v1'


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f'.{path.name}-{uuid4().hex}.tmp')
    temporary.write_text(json.dumps(data, indent=2)+'\n')
    os.replace(temporary, path)


def encode(value):
    if is_dataclass(value):
        return {'record':type(value).__name__, 'fields':{f.name:encode(getattr(value,f.name)) for f in fields(value)}}
    if isinstance(value,np.ndarray):
        return {'array':str(value.dtype),'data':encode(value.tolist())}
    if isinstance(value,np.generic):return encode(value.item())
    if isinstance(value,complex):return {'complex':[value.real,value.imag]}
    if isinstance(value,tuple):return {'tuple':[encode(v) for v in value]}
    if isinstance(value,list):return [encode(v) for v in value]
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    return value


def decode(value,module):
    if isinstance(value,list):return [decode(v,module) for v in value]
    if not isinstance(value,dict):return value
    if 'record' in value:
        cls=getattr(module,value['record'])
        if not isinstance(cls,type) or not is_dataclass(cls):raise ValueError('Unknown checkpoint record')
        return cls(**{k:decode(v,module) for k,v in value['fields'].items()})
    if 'array' in value:
        if value['array'] not in ('float64','complex128','int64','bool'):raise ValueError('Unexpected dtype')
        return np.asarray(decode(value['data'],module),dtype=value['array'])
    if 'complex' in value:return complex(*value['complex'])
    if 'tuple' in value:return tuple(decode(v,module) for v in value['tuple'])
    return {k:decode(v,module) for k,v in value.items()}


def fingerprint(module,config,inputs=None):
    data={'harness_version':HARNESS_VERSION,'config':asdict(config),'inputs':encode(inputs),
          'script_sha256':hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
          'green_sha256':hashlib.sha256(Path(module.lattice.__file__).read_bytes()).hexdigest(),
          'harness_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()


def manifest(module,config,stage):
    directory=Path(config.output_directory)
    atomic_json(directory/'run_manifest.json',{'run_id':directory.name,'stage':stage,
        'harness_version':HARNESS_VERSION,'config':asdict(config),
        'fingerprint':fingerprint(module,config),
        'assembly_version':getattr(module,'ASSEMBLY_VERSION','unchanged-original-assembly')})


def allocate_run(config, purpose='full', write_manifest=None, guard=None):
    """Shared allocation with theorem-specific manifest and path protection."""
    root=Path(config.output_directory)
    if guard is not None:guard(root)
    root.mkdir(parents=True,exist_ok=True)
    stamp=datetime.now(ZoneInfo('America/Bogota')).strftime('%Y-%m-%d_%H-%M-%S')
    directory=root/f'{stamp}_{purpose}_{uuid4().hex[:12]}'
    directory.mkdir(exist_ok=False)
    active=replace(config,output_directory=str(directory.resolve()))
    if write_manifest is not None:write_manifest(directory,active,'created')
    atomic_json(root/'latest_run.json',{'run_id':directory.name,'directory':active.output_directory,'purpose':purpose})
    print(f'Results for this execution: {directory.resolve()}',flush=True)
    return active


def new_run(module,config,purpose='full'):
    return allocate_run(config,purpose,
                        lambda directory,active,stage:manifest(module,active,stage))


def checkpoint(directory, signature, resume, calculate, encoder, decoder, valid):
    """Common atomic, fail-closed checkpoint lifecycle."""
    directory=Path(directory)
    directory.mkdir(parents=True,exist_ok=True)
    path=directory/'completed.json'
    if resume and path.exists():
        try:
            saved=json.loads(path.read_text());payload=saved['payload']
            if (saved.get('schema')==1 and saved.get('stage')=='completed'
                    and saved.get('fingerprint')==signature
                    and saved.get('sha256')==hashlib.sha256(payload.encode()).hexdigest()):
                result=decoder(json.loads(payload))
                if valid(result):return result
        except (OSError,ValueError,KeyError,TypeError,AttributeError):pass
    with (directory/'run.log').open('w',buffering=1) as log,redirect_stdout(log),redirect_stderr(log):
        result=calculate()
    if not valid(result):raise ValueError('Invalid checkpoint result structure')
    payload=json.dumps(encoder(result),separators=(',',':'))
    atomic_json(path,{'schema':1,'stage':'completed','fingerprint':signature,
                     'sha256':hashlib.sha256(payload.encode()).hexdigest(),'payload':payload})
    return result


def case(module,config,key,args,kwargs=None):
    kwargs=kwargs or {}
    directory=Path(config.output_directory)/'epsilon_cases'/key
    return checkpoint(directory,fingerprint(module,config,(args,kwargs)),
                      config.execution_resume,
                      lambda:module.validate_epsilon(*args,config=config,**kwargs),
                      encode,lambda value:decode(value,module),
                      lambda result:isinstance(result,tuple))


def _worker(name,configuration,epsilon,mu,nu):
    module=importlib.import_module(name)
    config=module.Config(**configuration)
    return case(module,config,f'epsilon_{epsilon:.8f}',(epsilon,mu,nu))


def independent_cases(module,config,mu,nu):
    epsilons=config.epsilon_values
    results={}
    workers=min(len(epsilons),max(1,config.execution_workers),max(1,(os.cpu_count() or 1)//2))
    if workers>1:
        with ProcessPoolExecutor(max_workers=workers,mp_context=mp.get_context('spawn')) as pool:
            jobs={pool.submit(_worker,Path(module.__file__).stem,asdict(config),e,mu,nu):e for e in epsilons}
            for job in as_completed(jobs):
                e=jobs[job];results[e]=job.result()
                print(f'epsilon cases completed {len(results)}/{len(epsilons)}: {e:.4f}',flush=True)
    else:
        for e in epsilons:results[e]=case(module,config,f'epsilon_{e:.8f}',(e,mu,nu))
    return results


def run_entrypoint(module,config,operation):
    active=new_run(module,config)
    manifest(module,active,'running')
    try:operation(active)
    except BaseException as exc:
        manifest(module,active,'interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed')
        raise
    manifest(module,active,'completed')
