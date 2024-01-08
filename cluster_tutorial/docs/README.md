# Welcome to cluster computing

This small repository is intended to help new user set up their own space on GenomeDK and get introduced to some basic use-cases.

- [Welcome to cluster computing](#welcome-to-cluster-computing)
  - [Python and Conda and the terminal](#python-and-conda-and-the-terminal)
  - [Logging on to the cluster](#logging-on-to-the-cluster)
  - [Your home on the cluster](#your-home-on-the-cluster)


## Python and Conda and the terminal

First things first.
If you do not already have python installed on your computer. I personally recommend just installing *miniconda*. It's a minimal version of something called *Anaconda*. It includes `python` along with the `conda` package manager and can easily be installed on Windows, Mac and Linux. To install miniconda go to this [page](https://docs.conda.io/projects/miniconda/en/latest/) and download the version fitting your computers operating system. (This is not strictly necessary if you do all your work on the cluster, but it is nice to have on your personal computer for any small programming tasks you might want to do)

Pretty much all software you will be using on the cluster is *command-line applications*, meaning programs which are executed through a *terminal* rather than a graphical user interface where you click on icons with a cursor. There are many options for what terminal to use, but all computers come with a built-in standard, which if you're not interested in anything with fancy bells and whistles, will function just fine. On Windows you have *PowerShell* and *CMD*, on Mac it is called *Terminal*, and for you cool kids with Linux it is also called *Terminal*.

When working on the cluster you need to install packages and programs for use in your analyses and pipelines, but sometimes the programs you need for one analysis may conflict with those needed for another. Luckily this is where `conda` helps saves the day! You see conda creates *environments* which work as isolated units, whatever is installed in one environment is completely isolated from whatever is going on in another environment. You can make as many environments as you want, they are easy to swithc between AND `conda` even helps make sure you have any necessary dependencies installed in any given environment!

## Logging on to the cluster

Now lets try logging on to the cluster. Open the terminal on your computer and write:

```cmd
ssh <USERNAME>@login.genome.au.dk
```

Where \<USERNAME> is your username on the cluster and press *ENTER*. You should now be prompted for your password, enter it and press *ENTER*.

```cmd
  _____                                ______ _   __
 |  __ \                               |  _  \ | / /
 | |  \/ ___ _ __   ___  _ __ ___   ___| | | | |/ /
 | | __ / _ \ '_ \ / _ \| '_ ` _ \ / _ \ | | |    \
 | |_\ \  __/ | | | (_) | | | | | |  __/ |/ /| |\  \
  \____/\___|_| |_|\___/|_| |_| |_|\___|___/ \_| \_/
  ```

As a new safety measure two-step authentication is now mandatory on the cluster, so if you haven't already you need to download the *Microsoft Authenticator* to your phone.

In the *Microsoft Authenticator* you should tap the '+' icon in the top right and choose 'other'. Now in your terminal, which is still logged on to the cluster, write:

```bash
gdk-auth-show-qr
```

This will show a QR code in your terminal, scan it with your phone.

To activate two-factor authentication write:

```bash
gdk-auth-activate-twofactor <TOKEN>
```

Where \<TOKEN> is the current token from the authenticator app. If everything went well you should now have activated two-factor authentication.

(This only applies if you received you password in an email. If you didn't skip to 'So currently, everytime...') Now, to be a little more annoying with all the safety pre-cautions, you should change your default password. In the terminal write:

```cmd
gdk-auth-change-password
```

Press *ENTER* and you will be prompted for your old password and what your new password should be.

So currently, everytime you log into the cluster you will have you write your password, BUT we can also just tell the cluster that your computer can be trusted, so let's do that!

First, to exit the cluster write:

```cmd
exit
```

When you press *ENTER* the terminal should return to your local computer. While on your local computer write:

```cmd
ssh-keygen -t rsa
```

This will generate a pair of authentication keys. When asked 'Enter file in whhich to save the key' just press *ENTER* and when it asks you to 'Enter a passphrase', again just press *ENTER*.

Now we will create a directory for our new public key on the cluster. Write the following:

```cmd
ssh <USERNAME>@login.genome.au.dk mkdir -p .ssh
```

Press *ENTER* and you will be prompted to enter your password, please do so.

Now we append the public `ssh` key on your local computer to the file `.ssh/authorized_keys` on the cluster. Write the following:

```cmd
cat ~/.ssh/id_rsa.pub | ssh username@login.genome.au.dk 'cat >> .ssh/authorized_keys'
```

Press *ENTER* and you will again be asked to enter you password (for the last time). From now on you can log into the cluster from your local computer without being prompted for a password.

## Your home on the cluster

Let's log back into the cluster, do you remember how? Just in case here is the command again:

```cmd
ssh <USERNAME>@login.genome.au.dk
```

Notice how you didn't have to write your password!

When you log into the cluster you start out in your private home folder. You can see this by the '~' on the prompt
![cluster_home](/images/cluster_home.png) If you want to have a look at what's inside your current folder your can write:

```bash
ls
```

For me it looks like this:
![folder_content](/images/content.png) You should see the text in different colors. The colors signify what you are looking at, like in my case where turqouis are shortcuts (or what called 'soft-links'), blue are folders and white are files.

A thing you might have noticed is that I have a *miniconda* folder and you don't, so now we need *miniconda* on the cluster. In the terminal write:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

This will download the latest version of *miniconda* to your current folder. when it's done write:

```bash
bash Miniconda3-latest-Linux-x86_64.sh
```

This runs the script and installs *miniconda*. Just follow the instructions and say 'yes' when it asks you if it should run `conda init`. When it's done you can check if there if a new folder with the `ls` command. You should also notice a small '(base)' next to your prompt. This indicates you current environment, which in this is the 'base' environment.
