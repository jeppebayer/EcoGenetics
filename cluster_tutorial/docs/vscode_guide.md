# VScode Guide

VS Code is great for developing scripts and editing text files.

Once you have installed VScode, you should install the *Remote Development* extension. You do that by clicking the funny squares in the left bar of VScode and search for 'Remote Development'. Once installed, you can click the small square with two arrows (kind of like '><') in the lower-left corner to connect to the cluster. Select *Connect current window to host* then *Add new SSH host*, then type:

```cmd
ssh <username>@login.genome.au.dk
```

Now select the config file `.ssh/config`. Now you can click the small green square in the lower-left corner to connect to the cluster by selecting login.genome.au.dk. It may take a bit, but once it is done installing a remote server, you will have access to the files in your home folder on the cluster.
