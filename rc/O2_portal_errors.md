# O2 Portal - R

These are common errors found when running R on the O2 portal and ways to fix them.

## Issue 1 - You can make a session and open Rstudio on O2 but cannot actually type.

Potential solution:  Make a new session and put the following under "Slurm Custom Arguments":  
```
-x compute-f-17-[09-25]
```

## Issue 2 - Everything was fine but then you lost connection.

When you attempt to reload you see:

<p align = "center">
<img src="../img/r_taking_longer.png">
</p>

Potential solutions: Refresh your interactive sessions page first then refresh your R page. 
If that doesn't work close your R session and re-open from the interactive sessions page.
If that doesn't work wait 5-10 min then repeat.

## Issue 3 - You made a session but cannot connect

When you attempt to connect you see:

<p align = "center">
<img src="../img/can_not_connect.png">
</p>

Potential solutions: This error indicates that either you did not load a gcc module or you loaded the incorrect one for the version of R you are running.
Kill the current session and start a new one with the correct gcc loaded in the modules to be loaded tab.
