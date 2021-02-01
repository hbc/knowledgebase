# O2 tips

## Making conda not slow down your login
If you have a complex base environment that gets loaded on login, you can end up having freezes of 30 seconds or more when
logging into O2. It is ultra annoying. You can fix this by not running the `_conda_setup` script in your .bashrc, like this:

```bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
#__conda_setup="$('/home/rdk4/local/share/bcbio/anaconda/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
#if [ $? -eq 0 ]; then
#    eval "$__conda_setup"
#else
if [ -f "/home/rdk4/local/share/bcbio/anaconda/etc/profile.d/conda.sh" ]; then
    . "/home/rdk4/local/share/bcbio/anaconda/etc/profile.d/conda.sh"
else
    export PATH="/home/rdk4/local/share/bcbio/anaconda/bin:$PATH"
fi
#fi
#unset __conda_setup
# <<< conda initialize <<<
```

## Interactive function to request memory and hours

Can be added to .bashrc (or if you don't want to clutter it, put it in .o2_aliases and then source it from .bashrc)

Defaults: 4G mem, 8 hours.
```
function interactive() {
        mem=${1:-4}
        hours=${2:-8}

        srun --pty -p interactive --mem ${mem}G -t 0-${hours}:00 /bin/bash
}
```

