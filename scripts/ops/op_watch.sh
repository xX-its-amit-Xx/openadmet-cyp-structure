#!/bin/sh
# Emit a line only when something actually MOVES. A monitor that reprints an unchanged
# count every cycle trains you to ignore it, which is how the 168-job batch got killed.
prev=""
while true; do
  c=$(python scripts/cofold/openprotein_cofold.py collect --engine protenix_v2 --tag op1 2>/dev/null \
      | python -c "import sys,json; d=json.load(sys.stdin); print(f\"folds {d['collected']} collected / {d['pending']} pending / {d['failed']} failed\")" 2>/dev/null)
  m=$(python scripts/cofold/p450_campaign.py msa-status 2>/dev/null \
      | python -c "import sys,json; d=json.load(sys.stdin); print('msa ' + ' '.join(f'{k}={v}' for k,v in sorted(d.items())))" 2>/dev/null)
  cur="$c | $m"
  if [ "$cur" != "$prev" ] && [ -n "$c" ]; then
    echo "$cur"
    prev="$cur"
  fi
  sleep 240
done
