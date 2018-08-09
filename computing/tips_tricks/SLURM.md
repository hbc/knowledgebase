* Useful aliases
```bash
alias bjobs='sacct -u ${USER} --format="JobID,JobName%25,NodeList,State,ncpus,start,elapsed" -s PD,R'
alias bjobs_all='sacct -u ${USER} --format="JobID,JobName%25,NodeList,State,ncpus,AveCPU,AveRSS,MaxRSS,MaxRSSTask,start,elapsed"'
```
