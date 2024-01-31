## For OSX

To have O2 accessible on your laptop/desktop as a folder, you need to use something called [`sshfs`](https://en.wikipedia.org/wiki/SSHFS) (ssh filesystem). This is a command that is not native to OSXand you need to go through several steps in order to get it. Once you have `sshfs`, then you need to set up ssh keys to connect O2 to your laptop without having to type in a password. 

### 1. Installing sshfs on OSX

Download OSXfuse from [https://github.com/osxfuse/osxfuse/releases](https://github.com/osxfuse/osxfuse/releases/download/osxfuse-3.6.0/osxfuse-3.6.0.dmg), and install it.

Download sshfs from [https://github.com/osxfuse/sshfs/releases](https://github.com/osxfuse/sshfs/releases/download/osxfuse-sshfs-2.5.0/sshfs-2.5.0.pkg), and install it.

> #### Use this only if the above option fails!
> 
> Step 1. Install [Xcode](https://developer.apple.com/xcode/)
> ```bash
> $ xcode-select --install
> ```
> 
> Step 2. Install Homebrew using ruby (from Xcode)
> ```bash
> $ /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
> 
> 	# Uninstall Homebrew
> 	# /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/uninstall)"
> ```
> 
> Step 2.1. Check to make sure that Homebrew is working properly
> ```bash
> $ brew doctor
> ```
> 
> Step 3. Install Cask from Homebrew's caskroom
> ```bash
> $ brew tap caskroom/cask
> ```
> 
> Step 4. Install OSXfuse using Cask
> ```bash
> $ brew cask install osxfuse
> ```
> 
> Step 5. Install sshfs from fuse
> ```bash
> $ brew install sshfs
> ```

### 2. Set up "ssh keys"

Once `sshfs` is installed, the next step is to connect O2 (or a remote server) to our laptops. To make this process seamless,  first set up ssh keys which can be used to connect to the server without having to type in a password everytime.

```bash
# set up ssh keys
$ ssh-keygen -t rsa -b 4096 -f ~/.ssh/id_rsa -C "ecommonsID"
$ ssh-add -K ~/.ssh/id_rsa
```

Arguments for `ssh-keygen`:
* `-t` = Specifies the type of key to create. The possible values are "rsa1" for protocol version 1 and "rsa" or "dsa" for protocol version 2. *We want rsa.*
* `-b` = Specifies the number of bits in the key to create. For RSA keys, the minimum size is 768 bits and the default is 2048 bits. *We want 4096*
* `-f` = name of output "keyfile"
* `-C` = Provides a new comment

Arguments for `ssh-add`:
* `-K` = Store passphrases in your keychain

```bash
# copy the contents of `id_rsa.pub` to ~/.ssh/authorized_keys on O2
$ cat ~/.ssh/id_rsa.pub | pbcopy
```

> `pbcopy` puts the output of `cat` into the clipboard (in other words, it is equivalent to copying with <kbd>ctrl + c</kbd>) so you can just paste it as usual with <kbd>ctrl + v</kbd>.


Log into O2 and use `vim` to open `~/.ssh/authorized_keys` and paste the contents copied from your computer to this file and save it. 


### 3. Mount O2 using sshfs

Now, let's set up for running `sshfs` on our laptops (local machines), by creating a folder with an intuitive name for your home directory on the cluster to be mounted in.

```bash
$ mkdir ~/O2_mount
```

Finally, let's run the `sshfs` command to have O2 mount as a folder in the above space.
```bash
$ sshfs ecommonsID@transfer.rc.hms.harvard.edu:. ~/O2 -o volname="O2" -o compression=no -o Cipher=arcfour -o follow_symlinks
```

Now we can browse through our home directory on O2 as though it was a folder on our laptop. 

> If you want to access your lab's directory in `/groups/` or your directory in `/n/scratch2`, you will need to create sym links to those in your home directory and you will be able to access those as well.

Once you are finished using O2 in its mounted form, you can cancel the connection using `umount` and the name of the folder.

```bash
$ umount ~/O2_mount 
```

### 4. Set up alias (optional)

It is optional to set shorter commands using `alias` for establishing and canceling `sshfs` connection. Use `vim` to create or open `~/.bashrc` and paste the following `alias` commands and save it.

```bash
$ alias mounto2='sshfs ecommonsID@transfer.rc.hms.harvard.edu:. ~/O2_mount -o volname="O2" -o follow_symlinks'
$ alias umounto2='umount ~/O2_mount'
```

> If your default shell is `zsh` instead of `bash`, use `vim` to create or open `~/.zshrc` and paste the `alias` commands.

Update changes in `.bashrc`

```bash
$ source .bashrc
```
Now we can type `mounto2` and `umounto2` to mount and unmount O2.

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
