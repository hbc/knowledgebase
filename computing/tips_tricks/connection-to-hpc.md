---
title: Connecting to hpc from local
description: This code helps with connecting to hpc computers
category: computing
subcategory: tips_tricks
tags: [ssh, hpc]
---


# osx

Use [Homebrew](http://brew.sh/) to get linux-like functionality on OSX

Use [XQuartz](https://www.xquartz.org/) for X11 window functionality in OSX.

# Odyssey with 2FA
Enter one time password into the current window (https://github.com/jwm/os-x-otp-token-paster)

# Fix 'Warning: No xauth data; using fake authentication data for X11 forwarding'
Add this to your ~/.ssh/config on your OSX machine:

```
Host *
  XAuthLocation /opt/X11/bin/xauth
```

# Use ssh keys on remote server
This will add your key to the OSX keychain, here your private key is assumed to be named "id_rsa":

```
ssh-add -K ~/.ssh/id_rsa
```

Now tell ssh to use the keychain. Add this to the ~/.ssh/config on your OSX machine:

```
Host *
     AddKeysToAgent yes
     UseKeychain yes
     IdentityFile ~/.ssh/id_rsa
     XAuthLocation /opt/X11/bin/xauth
```
