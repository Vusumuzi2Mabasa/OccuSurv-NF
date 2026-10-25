"""Keep the original generator and QC; perform no statistical fitting."""
import sys,json,hashlib
from pathlib import Path
import pandas as pd
from longitudinal import generate,validate
out=Path('longitudinal_input');out.mkdir(exist_ok=True)
generate(sys.argv[1],int(sys.argv[2]),out/'raw')
raw=out/'raw'
qc=validate(*[pd.read_csv(raw/f'{n}.csv') for n in ['intervals','subjects','history']])
qc['sha256']={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in raw.iterdir()}
(out/'qc.json').write_text(json.dumps(qc,indent=2))
