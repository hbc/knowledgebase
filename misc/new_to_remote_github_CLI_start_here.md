# Putting local git repos into the HBC Github organization remotely via command line interface

Has your experience with github primarily been through the browser? This document has the basics to begin turning your working directories into github repositories which can be pushed to the HBC location remotely via command line.

### Wait, back up, what do you mean by push?

Confused by push, pull, branch, main, commit? If you're not sure, it's worthwhile to familiarize yourself with the basics of git/github. There are some great resources and tutorials to learn from out there. Here's an interactive one (I could only get it to work in safari, not chrome):
https://learngitbranching.js.org

This won't teach you how to put your things into the HBC Github organization though.

## Set up/configuration

You only need to do these once!

### 1. Configure git locally to link to your github account
Open up a terminal and type

```bash
git config --global user.email EMAIL_YOU_USE_TO_SIGN_IN_TO_GITHUB
git config --global user.name YOUR_GITHUB_USERNAME
```

### 2. Make personal access token
Configuring your local git isn't enough, as github is moving away from passwords. You will need to make a personal access token through your github account. Follow the instructions here:
https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens

**Copy this personal access token and save it somewhere or keep the window open for now. You will be prompted to enter your personal access token the first time you type `push -u origin main`. You will not be able to access this token again!**

## Creating your git repo

### 1. Initialize git repo on the HBC Github via web browser

I have yet to find a way to do this remotely via CLI, but as far as I can tell this step is a necessary pre-requisite to pushing a local repo to the HBC github. Will update to add CLI if possible.

Go to https://github.com/hbctraining and click the green "New Repository" button. Initialize a new, empty repository.

Once you do this, there will be some basic code you can copy under `Quick Setup`, including the `https` location of your repo which can be used below
   
### 2. Create a local git repo and push it to the HBC Github via CLI

In your terminal, navigate to the folder you would like to turn into a github repo and type the following:

```bash
echo "# text to add to readme" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/hbctraining/NAME_OF_REPOT.git
git push -u origin main
```
**You will be prompted to enter your personal access token the first time you type `push -u origin main`.**

## Useful tips/tricks

If you are doing this in a directory with folders/data/files you don't necessarily want to put on github, you will want to pick/choose what you upload. Here are some tips and notes:

### Add all, but exclude some

**note: will not exclude if already pushed to HBC repo! Just untracks them!**

The best time to implement this is when you are making your first upload

```bash
git add .
git reset -- path/to/thing/to/exclude
git reset -- path/to/more/things/to/exclude
git commit -m "NAME_OF_COMMIT"
git push -u origin main
```

### To add specific files/folders:

```bash
git add path/to/files*
git commit -m "NAME_OF_COMMIT"
git push -u origin main
```

### Add all files/folder except this file/folder:

```bash
git add -- . ':!THING_TO_EXCLUDE'
git commit -m "NAME_OF_COMMIT"
git push -u origin main
```

### Remove a file you already pushed to Github

You might be tempted to just do this in the browser, but be warned! It will break your local repo until you pull from the HBC location. This could be a problem if that was important data you want to continue to store locally. Fortunately, you can "delete" files on the HBC Github without deleting them locally. Here's how:

```bash
git rm --cached NAME_OF_FILE
```

Or for a folder:
```bash
git rm -r --cached NAME_OF_FILE
```
**Side effects of resorting to `git rm -r --cached REALLY_IMPORTANT_DATA_DIRECTORY` include anxiety, sweating, heart palpitations, and appeals to spiritual beings**

### Check what changes have been made to the current commit

Very useful to see what will actually be added, removed, etc or if everything is up to date.
```bash
git status
```

### .gitignore: coming soon

