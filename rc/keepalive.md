---
title: Transfer files inside cluster
description: This code helps with transfer files inside cluster
category: computing
subcategory: tips_tricks
tags: [ssh, hpc]
---

Useful for file transfers on O2's new transfer cluster (transfer.rc.hms.harvard.edu).

The nohup command can be prepended to the bash command and the command will keep running after you logout (or have your connection interrupted).

From HMS RC:

`From one of the file transfer systems under transfer.rc.hms.harvard.edu , you can prefix your command with "nohup" to put it in the background and be able to log out without interrupting the process.`

`For example, after logging in to e.g. the transfer01 host, run your command:`

`nohup rsync -av /dir1 /dir2`

`and then log out. rsync will keep running.`

`To check in on the process later, just remember which machine you ran rsync and you can directly re-login to that system if you like.`

`For example:`

`1. ssh transfer.rc.hms.harvard.edu (let's say you land on transfer03), and then:`
`2. ssh transfer01`
`-- from there you can run the "ps" command or however you like to monitor the process.`


## Another option from John
If you run tmux from the login node before you ssh to the transfer node to xfer files, you can drop your connection and then re-attach to your tmux session later. It should still be running your transfer. 

**General steps**
1) Login to O2
2) write down what login node your are on (usually something like login0#)  
*at login node*
3) Start a new tmux session  
`tmux new -s myname`
4) SSH to the transfer node  
`ssh user@transfer.rc.hms.harvard.edu`  
*on transfer node*
5) start transfer with rsync, scp etc.
6) close terminal window without logging out  
*time passes*
7) Login to O2 again
8) ssh to the login node you wrote down above  
`ssh user@login0#`
9) Reattach to your tmux session  
`tmux a -t myname`
10) Profit

