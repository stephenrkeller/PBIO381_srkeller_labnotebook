# P/BIO381 Tutorials

## Intro to connecting to unix servers and navigating the command-line

### February 1, 2017



*What is the command-line?*

The command-line, also known as a "terminal" or "shell", is a way of interacting with your local computer or a remote server by means of typing commands or scripts, without using a graphical user interface (GUI).



*Why do I want to be doing this?*

At first, the command-line can seem a little intimidating. But after you get used to typing instead of pointing and clicking to issue your commands, you'll realize how powerful it is. For example, it's quite easy to copy, move, edit, and search within thousands of files in multiple directories with some simple command-line code. It would take forever to do this by dragging/dropping with a mouse. The command-line also allows you to work with very large data files without uncompressing them fully, or loading the entire file's contents into memory…something that standard GUI type applications aren't good at.



*So, let's get started…*

- The first step is to open a terminal *shell* on your local computer. For windows users, this would be "PuTTy". For MacOS users, this is called "Terminal".

- We'll connect to our remote server running Unix using the secure shell (ssh) protocol. Our server's name is *pbio381* and we can connect to it using our UVM netid username and password (as long as we're on-campus)

```bash
ip0af52fbf:papers srkeller$ ssh srkeller@pbio381.uvm.edu
srkeller@pbio381.uvm.edu's password: 
Last login: Tue Jan 31 10:51:15 2017 from ip040027.uvm.edu
[srkeller@pbio381 ~]$ 
```



- The log-in screen tells us some basic info on when we last logged in, and then gives us our current location in the filesystem (~) followed by the $ prompt, that tells us the computer is ready for our next command. 

  - NOTE: The tilda (~) is short-hand for your home directory in UNIX. This is your own personal little corner of the computer's hard drive space, and is the location that you should use to create folders and input/output/results files that are specific to your own work. No one has access to any files stored in your home directory but you.

- To see the full path to your current directory, use the **pwd** command:

  ```bash
  [srkeller@pbio381 ~]$ pwd
  /users/s/r/srkeller
  [srkeller@pbio381 ~]$ 
  ```

  ​

- The path shows the full directory address up from the "root" of the file structure, which is the most basal level (appealing to all you phylogeneticists here…). The root is symbolized as "/" and each subdirectory is separated by an additional "/". So, the full path to my working directory on the server is */users/s/r/srkeller/*


- Let's make a new folder (aka, directory) using the **mkdir** command. Let's name this folder "mydata"

```bash
[srkeller@pbio381 ~]$ mkdir mydata
```

- We can then use the **ll** command to show the current contents of any folders and files in our current location:

```bash
[srkeller@pbio381 ~]$ ll
total 0
drwxr-xr-x. 6 srkeller users 82 Jan 31 17:21 archive
drwxr-xr-x. 2 srkeller users  6 Jan 31 17:21 mydata
drwxr-xr-x. 2 srkeller users  6 Jan 31 17:15 scripts
[srkeller@pbio381 ~]$ 
```

- You'll notive that I've got some extra folders in my output from previous work, whereas you will probably only see the "scripts" folder you just made. 
  - NOTE: Each row shows a file or a folder (in this case, these are all folders) diplaying (from right to left) its name, when it was last edited, size, who it belongs to , and who has permission to read (r) write (w) and exectue (x) it. More on permissions later...
  - Try making your own folder named "scripts" and then use the **ll** command to list the folders again 
- We can change our current location within the directory structure using the **cd** command. Let's use **cd** to move inside the *mydata/* directory and **ll** to list its contents:

```bash
[srkeller@pbio381 ~]$ cd mydata/
[srkeller@pbio381 mydata]$ ll
total 0
[srkeller@pbio381 mydata]$ 
```

- Hah — nothing in there yet! Let's go get some data!
  - We've placed the text file containing all the metadata information on the seastar sampling under a shared space on the server. The path to this shared space is: 
    - */data/*   Try using **cd** to navigate over to this location. Then **ll** to show its contents. You should see something like this:

```bash
drwxr-xr-x.  5 root root       73 Jan 31 17:35 archive
drwxrwxr-x.  2 root pb381adm   40 Nov 30  2015 packages
drwxrwxr-x. 33 root pb381adm 4096 Nov 30  2015 popgen
drwxrwxr-x.  3 root pb381adm   42 Jan 30 09:08 project_data
drwxrwxr-x.  2 root pb381adm    6 Oct  2  2015 scripts
drwxr-xr-x. 18 root root     4096 Sep  2  2015 users
[srkeller@pbio381 data]$ 
```

- Now, **cd** into the folder called "project_data" and **ll**. Do you see this?

```bash
[srkeller@pbio381 data]$ cd project_data/
[srkeller@pbio381 project_data]$ ll
total 8
drwxr-xr-x. 12 srkeller users 4096 Jan 30 09:06 archive
-rw-r--r--.  1 srkeller users 1255 Jan 30 09:08 ssw_samples.txt
[srkeller@pbio381 project_data]$ 
```

- The file called "ssw_samples.txt" is the one with the seastar metadata. We don't want to open and make changes to this file in the shared space, because we don't want to have our edits affect the rest of the group. So, let's first make a copy of this file over to our home directory and put it inside the "mydata" folder. Use the **cp** command, followed by the filename, and the path to your destination (remember the ~ signals your home directory, and each subdirectory is then separated by a /):

```bash
[srkeller@pbio381 project_data]$ cp ssw_samples.txt ~/mydata/
```

- **cd** back to your *~/mydata/* directory and look inside. You should see your file...

```bash
[srkeller@pbio381 project_data]$ cd ~/mydata/
[srkeller@pbio381 mydata]$ ll
total 4
-rw-r--r--. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
[srkeller@pbio381 mydata]$ 
```

