this post is for those who are interested in running notebooks seemingly on O2. There is a well written documentation about running jupyter notebook on O2. you can find it here. https://wiki.rc.hms.harvard.edu/display/O2/Jupyter+on+O2. However, this involves multiple steps, opening bunch of terminals at times, and importantly finding an unused port every-time. I found it quite cumbersome and annoying, so i spent some time solving it. It took me a while to nail it down, with the help of FAC RC, but they suggested a simpler solution. If you wish to run jupyter/R notebook on O2 {where your data sits},
Here what you need to do : 

Install https://github.com/aaronkollasch/jupyter-o2 by `pip install jupyter-o2` on your local terminal. 

Run `jupyter-o2 --generate-config` on command line. 
This will generate the configuration file and will tell you where it is located. Un comment the fields which are not needed. Since configuration file is the key, find attach the a template for use, you need to change your credentials though.

You are all set to run notebook on O2 from your local machine now, without logging into server.  
Now at your local terminal Run `jupyter-o2 notebook` for python notebooks. Alternatively you can also do `jupyter-o2 lab` for R/python
This will ask you a paraphrase, you should enter your ecommons password as paraphrase. 
Boom!!! you are good to go! Happy Pythoning :):)  
If you wish you run R notebooks on O2, refer this. https://docs.anaconda.com/anaconda/navigator/tutorials/r-lang/
