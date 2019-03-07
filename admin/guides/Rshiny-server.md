---
title: How to set up new app in hbc rshiny
description: Information on where to connect and where to put apps in hbc rshiny server.
category: admin
subcategory: guides
tags: [visualization, report]
---

To get access to the server, ask to the Director or Associate Director of the core

The server is only reachable from HMS network. The typical way to connect is from the login node in O2. Type this:

```
ssh hbcreports.med.harvard.edu
```

Use the password given to you by HMS RC team.

Once inside, copy the folder that contains the code for you app inside `/srv/shiny-server/`, from that point, the app should be available at http://hbcreports.med.harvard.edu/`app_name`/. For instance, look at this: http://hbcreports.med.harvard.edu/heatmaps/


**installation of new packages affect to all current apps. This need to be coordinated with the responsible person in the core. Raise any need to install or update during group meeting.**
