"""Synthetic notification sink for the Slurm channel-timing canary only."""
import argparse
from datetime import datetime,timezone
import json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--bam-ready',required=True)
a=p.parse_args();c=json.loads(Path(a.config).read_text())
Path(c['receipt']).write_text(json.dumps(dict(received_at=datetime.now(timezone.utc).isoformat(),
    metrics=a.bam_ready,completed_output_present=(Path(a.bam_ready)/'event.txt').is_file()))+'\n')
