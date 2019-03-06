---
title: DS_Store tips
description: Don't leave .DS_Store files on network volumes
category: computing
subcategory: tips_tricks
tags: [osx]
---

# Don't leave .DS_Store files on network volumes
defaults write com.apple.desktopservices DSDontWriteNetworkStores true
