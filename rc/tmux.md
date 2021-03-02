Tmux is a great way to work on the server as it allow syou to:  
1) keep your session alive.  
3) have multiple named sessions open to firewall tasks/projects.  
4) run multiple windows/command lines from a single login (O2 allows a maximum of 2-3 logins to their system).  
5) quickly spin off windows to do small commands (see 3).  

### Useful resources

#### [Tmux cheat sheet](https://tmuxcheatsheet.com/)


#### Tmux configuration
Tmux works great but can have some issues upon first use that make it challenging to use for those of us used to a GUI:  
a) the default command key is not great.  
b) it doesn't work well with a mouse.  
c) it doesn't let you copy text easily.  
d) it doesn't scroll your window easily.  
e) resizing the windows can be challenging.  

The confirguation file code below should make some of these issues easier:

    set -g default-terminal "screen-256color"
    set -g status-bg red
    set -g status-fg black
    
    # Use <C-a> instead of the default <C-b> as Tmux prefix
    set-option -g prefix C-a
    unbind-key C-b
    bind-key C-a send-prefix
    
    
    # Options enable mouse support in Tmux
    #set -g terminal-overrides 'xterm*:smcup@:rmcup@'
    # For Tmux >= 2.1
    #set -g mouse on
    # For Tmux <2.1
    # Make mouse useful in copy mode
    setw -g mode-mouse on
    #
    # # Allow mouse to select which pane to use
    set -g mouse-select-pane on
    #
    # # Allow mouse dragging to resize panes
    set -g mouse-resize-pane on
    #
    # # Allow mouse to select windows
    set -g mouse-select-window on
    
    
    # set colors for the active window
    # START:activewindowstatuscolor
    setw -g window-status-current-fg white 
    setw -g window-status-current-bg red 
    setw -g window-status-current-attr bright
    # END:activewindowstatuscolor
    
    
    ## Optional- act more like vim:
    #set-window-option -g mode-keys vi
    #bind h select-pane -L
    #bind j select-pane -D
    #bind k select-pane -U
    #bind l select-pane -R
    #unbind p
    #bind p paste-buffer
    #bind -t vi-copy v begin-selection
    #bind -t vi-copy y copy-selection 
    
    
    # moving between panes
    # START:paneselect
    bind h select-pane -L 
    bind j select-pane -D 
    bind k select-pane -U
    bind l select-pane -R    
    # END:paneselect
    
    
    # START:panecolors
    set -g pane-border-fg green
    set -g pane-border-bg black
    set -g pane-active-border-fg white
    set -g pane-active-border-bg yellow
    # END:panecolors
    
    # Command / message line
    # START:cmdlinecolors
    set -g message-fg white
    set -g message-bg black
    set -g message-attr bright
    # END:cmdlinecolorsd -g pane-active-border-bg yellow
    
    
    bind-key C-a last-window
    setw -g aggressive-resize on
