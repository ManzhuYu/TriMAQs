"""Resolve data paths and workflow settings from one configuration file."""
from pathlib import Path
import os,json

def load():
    filename=os.environ.get('TRIMAQS_CONFIG')
    if not filename:raise RuntimeError('Use run.py --config config/settings.json')
    p=Path(filename).expanduser().resolve();cfg=json.loads(p.read_text(encoding='utf-8-sig'))
    root=Path(os.environ.get('TRIMAQS_DATA_ROOT',cfg['data_root'])).expanduser()
    cfg['_root']=root if root.is_absolute() else (p.parent/root).resolve()
    return cfg

def path(key):
    cfg=load();value=Path(cfg['paths'][key]).expanduser()
    return str(value if value.is_absolute() else cfg['_root']/value)

def setting(key):return load()['settings'][key]

def relative_metadata(value):
    """Use dataset-relative file references in published metadata."""
    if isinstance(value,dict):return {key:relative_metadata(item) for key,item in value.items()}
    if isinstance(value,(list,tuple)):return [relative_metadata(item) for item in value]
    if isinstance(value,str):
        root=str(load()['_root'])
        if value==root:return '.'
        for separator in ['/',chr(92)]:
            if value.startswith(root+separator):return value[len(root)+1:].replace(chr(92),'/')
        p=Path(value)
        if p.is_absolute():return 'external/'+p.name
    return value
