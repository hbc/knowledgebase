## My intitial notes on using the FAS-RC Open on Demand virtual desktop system

### Steps to get going
- Download and install Cisco connect VPN

- Login to to VPN using your FAS-RC credentials and two-pass authentication code

- Navigate to https://vdi.rc.fas.harvard.edu/pun/sys/dashboard

- Login to page using FAS-RC user id and password

- Click on Interactive Apps pulldown and select the “Rstudio Server” under the SEr5ver heading. DO NOT select, “RStudio Server (bioconductor + tidyverse)”

Here most of the settings are self-explanatory.   
- I have tried out multiple cores (12) and while I am unsure it is using the full 48 cores parallel::detectCores finds, my simulation did run substantially faster (4 cores = 5% done when “12 cores” was at 75%. I would be  interested to hear other people’s experiences and how transparently/well it works with R.

- Maximum memory allocation is supposed to bey 120GHB. I haven’t tried asking for more. 

- I’ve been loading the R/4.02-fasrrc01 Core R versions.  I tried the  R/4.02-fasrc Core & gcc 9.3.0 first and ran into package compilation issues. 

- You can set the R_LIBS_USER folder to use which will  contain your library packages. Using this approach, I was able to install packages in session , delete the server and come back to the installed packages in a new session.  You could theoretically also switch between R versions using this and the version selector. 

- I haven’t tired executing a script before staring Rstudio, but theoretically, I could see using this to launch a condo environment. 

- I don’t know about reservations but they sound interesting for getting a high mem machine.
