import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
for eq in json.loads((root/'docs/methods.json').read_text()):
 fig=plt.figure(figsize=(10,1));fig.text(.02,.5,'$'+eq['latex']+'$',fontsize=16,va='center');fig.savefig(root/'assets/equations'/(eq['id']+'.png'),dpi=250,bbox_inches='tight',pad_inches=.12);plt.close(fig)
