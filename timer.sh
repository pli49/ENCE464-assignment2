echo "Timing One"

{ time ./runner.sh; } 2>&1 1>/dev/null | awk '/real/ {print "First Experiment: " $2}'
