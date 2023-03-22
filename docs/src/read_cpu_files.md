# Read CPU log files

Gadget writes a number of CPU log files. We provide some parsers to work with their output.

## balance.txt

The balance file contains information on the wallclock time of each timestep and how many particles are active during that step.
You can get this information with

```@docs
parse_balance
```

This automatically skips over restarts and only outputs unique timesteps, i.e. if you start from a restartfile and Gadget has to re-do some of the timesteps, this function will only output the second iteration.