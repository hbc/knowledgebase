---
title: IPython notebook on O2
description: How to open up an ipython notebook running on O2
category: computing
subcategory: tips_tricsk
tags: [python, ipython, singlecell]
---

1. First connect to O2 and open up an interactive session with all of the cores and memory you want to use. Here I'm connecting to the short queue so I can get more cores to use.

```bash
srun -n 8 --pty -p short --mem 64G -t 0-12:00 --x11 /bin/bash
```

2. Note the name of the compute node you are on:

```bash
uname --nodename
```

3. Start a jupyter notebook server on a specific port:

```bash
jupyter notebook --no-browser --port=1234
```

This command will open up a notebook server on port 1234. You might have to pick
a different port if 1234 is being used. Note the token it provides for you, you
will need this token to use your notebook server.

4. Create a SSH tunnel from your local machine to the jupyter notebook:

On your local machine do:

```
ssh -L 9999:localhost:9998 134.174.159.22 ssh -L 9998:localhost:1234 -N compute-e-16-238
```

This sets up a two SSH tunnels. The first one is connecting port 9999 on your laptop to port 9998 on `login02` on o2 (134.174.159.22). The second is connecting port 9998 on `login02` to port 1234 on `compute-e-16-238`.

5. Open a web browser and put `localhost:9999` as the address.

This should now connect you to the jupyter notebook server. It will ask you for
the token. If you put the token in, you can now log in and will be in your
home directory on O2.

You are now running a notebook server. This is just running using a single core now-- we want to hook up our computing that we reserved. We asked for 8 cores, so we'll set up a cluster with 8 cores. Click on "IPython Clusters", set the number of engines on
"default" to 8, and you will have your notebook connected to the 8 cores.

6. Start working! 

You can open up a terminal by going to the Files tab and clicking on new and opening 
the terminal. You can start a new notebook by going to the Files tab, clicking on 
new and opening a python notebook.
