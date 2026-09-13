#!/bin/sh
# Collect every engine, not just the first one. A watcher that polls one engine reports
# "0 pending" while three other campaigns sit uncollected on the server.
# Emits only when something MOVES: a monitor that reprints an unchanged count every cycle
# trains you to ignore it, which is how a healthy 168-job batch once got killed.
prev=""
while true; do
  line=""
  for e in protenix_v2 protenix esmfold2 rosettafold_3; do
    r=$(python scripts/cofold/openprotein_cofold.py collect --engine "$e" --tag op1 2>/dev/null \
        | python -c "import sys,json
try:
    d=json.load(sys.stdin); print(f\"{d['collected']}/{d['pending']}\")
except Exception: print('-')" 2>/dev/null)
    line="$line $e:$r"
  done
  p=$(python scripts/cofold/p450_campaign.py collect 2>/dev/null \
      | python -c "import sys,json
try:
    d=json.load(sys.stdin); print(f\"p450:{d['collected']}/{d['pending']}\")
except Exception: print('p450:-')" 2>/dev/null)
  m=$(python scripts/cofold/p450_campaign.py msa-status 2>/dev/null \
      | python -c "import sys,json
try:
    d=json.load(sys.stdin); print('msa:'+str(d.get('SUCCESS',0)))
except Exception: print('msa:-')" 2>/dev/null)
  cur="collected/pending -$line $p $m"
  if [ "$cur" != "$prev" ]; then
    echo "$cur"
    prev="$cur"
  fi
  sleep 300
done
