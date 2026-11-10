"""Run one selected workflow stage."""
from pathlib import Path
import argparse,json,os,subprocess,sys,ast

def main():
    root=Path(__file__).resolve().parent;stages=json.loads((root/'stages.json').read_text())
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=sorted(stages));parser.add_argument('--config',required=True);parser.add_argument('--check',action='store_true')
    args=parser.parse_args();config=Path(args.config).expanduser().resolve();os.environ['TRIMAQS_CONFIG']=str(config)
    script=root/stages[args.stage];ast.parse(script.read_text(encoding='utf-8'))
    if args.check:
        from trimaqs_config import load,path
        cfg=load();print('Data root:',cfg['_root']);print('Stage:',args.stage)
        for key in cfg['paths']:
            if key.startswith(args.stage+'.'):print(key,'->',path(key))
        print('Configuration and Python syntax checked; no models were run.');return
    raise SystemExit(subprocess.run([sys.executable,str(script)],env=os.environ.copy(),check=False).returncode)

if __name__=='__main__':main()
