# P/BIO381 Tutorials

## Intro to connecting to Unix servers and navigating the command-line

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

- Let's take a peek at this file with the **head** command, which prints the first 10 lines to screen.

```bash
[srkeller@pbio381 mydata]$ head ssw_samples.txt 
Individual	Trajectory	Location	Day3	Day6	Day9	Day12	Day15
10	HH	INT	10_5-08_H	10_5-11_H	10_5-14_H	10_5-17_H	10_5-20_H
24	HH	INT	24_5-08_H	24_5-11_H	24_5-14_H	24_5-17_H	24_5-20_H
27	HH	INT	27_5-08_H	27_5-11_H	27_5-14_H	27_5-17_H	27_5-20_H
08	HS	INT	08_5-08_H	08_5-11_S	08_5-14_S	08_5-17_S	08_5-20_S
09	HS	INT	09_5-08_H		09_5-14_S	09_5-17_S	09_5-20_S
15	HS	INT	15_5-08_H	15_5-11_H	15_5-14_H	15_5-17_S	15_5-20_S
19	HS	INT		19_5-11_H	19_5-14_H	19_5-17_H	19_5-20_S
20	HS	INT	20_5-08_H	20_5-11_H	20_5-14_H	20_5-17_H	20_5-20_S
03	SS	INT	03_5-08_S	03_5-11_S
```

- The **tail** command provides similar functionality, but prints just the last lines in the file. These features may not seem a big deal right now, but when you're dealing with files that are 20 Gb compressed, and feature hundreds of millions of lines of data, you and your computer will be happy to have tools to peek inside without having to open the whole file!
- What if we want to extract just the rows of data that correspond to Healthy (HH) individuals? We can use the search tool **grep** to search for a target query. Any line matching our search string will be printed to screen.

```bash
[srkeller@pbio381 mydata]$ grep 'HH' ssw_samples.txt 
10	HH	INT	10_5-08_H	10_5-11_H	10_5-14_H	10_5-17_H	10_5-20_H
24	HH	INT	24_5-08_H	24_5-11_H	24_5-14_H	24_5-17_H	24_5-20_H
27	HH	INT	27_5-08_H	27_5-11_H	27_5-14_H	27_5-17_H	27_5-20_H
31	HH	SUB	31_6-12_H	31_6-15_H	31_6-18_H	31_6-21_H	31_6-24_H
32	HH	SUB	32_6-12_H	32_6-15_H	32_6-18_H	32_6-21_H	
33	HH	SUB	33_6-12_H	33_6-15_H	33_6-18_H	33_6-21_H	33_6-24_H
34	HH	SUB	34_6-12_H	34_6-15_H	34_6-18_H	34_6-21_H	34_6-24_H
35	HH	SUB	35_6-12_H	35_6-15_H	35_6-18_H	35_6-21_H	
[srkeller@pbio381 mydata]$
```

- What if instead of printing it to screen, we want to save the output of our search to a new file? This is easy, just use the ">" symbol to redirect the results of any command to an output file with your choice of name.

```bash
[srkeller@pbio381 mydata]$ grep 'HH' ssw_samples.txt >ssw_HHonly.txt
[srkeller@pbio381 mydata]$ ll
total 8
-rw-r--r--. 1 srkeller users  462 Jan 31 20:46 ssw_HHonly.txt
-rwxrwxr-x. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
[srkeller@pbio381 mydata]$ 
```

- We can do the same routine for the "SS" samples. Here's a trick, when you're doing a similar task as a previous command, hit the up arrow on your keyboard at the $ prompt, and it will recall the last command you issued. Then you just have to switch the HH's for SS's.

```bash
[srkeller@pbio381 mydata]$ grep 'SS' ssw_samples.txt >ssw_SSonly.txt
[srkeller@pbio381 mydata]$ ll
total 12
-rw-r--r--. 1 srkeller users  462 Jan 31 20:46 ssw_HHonly.txt
-rwxrwxr-x. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
-rw-r--r--. 1 srkeller users  342 Jan 31 20:48 ssw_SSonly.txt
[srkeller@pbio381 mydata]$ 
```

- **Grep** is a useful search tool and has many additional features for sorting and output of the results. These kinds of search algorithms are called "regular expressions", or "regexp", and are one of the most powerful tools for wokring with large text files. If you want to learn more about **grep** and its regexp capabilities, you can look at the **"man"** page or manual. In fact, every UNIX command-line program has a built-in **man** page that you can call up to help you. Just type **man** and then the program name and it will give you the manual (small excerpt shown below).

```bash
[srkeller@pbio381 mydata]$ man grep


GREP(1)                            General Commands Manual                           GREP(1)

NAME
       grep, egrep, fgrep - print lines matching a pattern

SYNOPSIS
       grep [OPTIONS] PATTERN [FILE...]
       grep [OPTIONS] [-e PATTERN | -f FILE] [FILE...]

DESCRIPTION
       grep searches the named input FILEs (or standard input if no files are named, or if a
       single hyphen-minus (-) is given as file name) for lines containing a  match  to  the
       given PATTERN.  By default, grep prints the matching lines.

       In  addition,  two variant programs egrep and fgrep are available.  egrep is the same
       as grep -E.  fgrep is the same as grep -F.  Direct  invocation  as  either  egrep  or
       fgrep  is  deprecated,  but is provided to allow historical applications that rely on
       them to run unmodified.

OPTIONS
   Generic Program Information
       --help Print a usage message briefly summarizing these command-line options  and  the
              bug-reporting address, then exit.

       -V, --version
              Print  the version number of grep to the standard output stream.  This version
              number should be included in all bug reports (see below).

   Matcher Selection
       -E, --extended-regexp
              Interpret PATTERN as an extended regular expression (ERE, see below).  (-E  is
              specified by POSIX.)

       -F, --fixed-strings, --fixed-regexp
              Interpret  PATTERN  as  a list of fixed strings, separated by newlines, any of
              which is to be matched.  (-F is  specified  by  POSIX,  --fixed-regexp  is  an
              obsoleted alias, please do not use it in new scripts.)

       -G, --basic-regexp
              Interpret PATTERN as a basic regular expression (BRE, see below).  This is the
              default.

       -P, --perl-regexp
              Interpret PATTERN as a Perl regular expression.  This is  highly  experimental
              and grep -P may warn of unimplemented features.
```

- One of the most useful aspects of UNIX is the ability to take the output from one command and use it as standard input (termed 'stdin') into another command without having to store the intermediate files. Such a workflow is called "piping", and makes use of the pipe character (|) located above the return key to feed data between programs.
  - Example: Say we wanted to know how many samples come from the Intertidal. We can use **grep** to do the search, and pipe the results to the command **wc** which will tally up the number of lines, words, and characters in the file…voila!

```bash
[srkeller@pbio381 mydata]$ grep 'INT' ssw_samples.txt | wc
     16     106     762
[srkeller@pbio381 mydata]$ 
```

- Looks like 16 INT samples in the original data. See how quick it was to get a line count on this match, without actully opening a file or printing/saving the outputs? 
- Now, what if we want to move the files we created with just individuals of a particular disease status. There's a way to do this quickly using the wildcard character "*". With the wildcard, the "*\*" takes the place of any character, and in fact any length of characters. For example, make a new directory called *samples_by_disease/* inside the *mydata/* folder. Then move all files that contain the word "only" into the new directory using the **mv** command.

```bash
[srkeller@pbio381 mydata]$ mkdir sample_by_disease/
[srkeller@pbio381 mydata]$ ll
total 12
drwxr-xr-x. 2 srkeller users   10 Jan 31 21:12 sample_by_disease
-rw-r--r--. 1 srkeller users  462 Jan 31 20:46 ssw_HHonly.txt
-rwxrwxr-x. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
-rw-r--r--. 1 srkeller users  342 Jan 31 20:48 ssw_SSonly.txt
[srkeller@pbio381 mydata]$ mv *only* sample_by_disease/
[srkeller@pbio381 mydata]$ ll
total 4
drwxr-xr-x. 2 srkeller users   60 Jan 31 21:12 sample_by_disease
-rwxrwxr-x. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
[srkeller@pbio381 mydata]$ cd sample_by_disease/
[srkeller@pbio381 sample_by_disease]$ ll
total 8
-rw-r--r--. 1 srkeller users 462 Jan 31 20:46 ssw_HHonly.txt
-rw-r--r--. 1 srkeller users 342 Jan 31 20:48 ssw_SSonly.txt
[srkeller@pbio381 sample_by_disease]$ 
```

- OK, what about when we have files we don't want anymore? How do we clean up our workspace? You can remove files and folders with the **rm** command. However, in its default mode, UNIX will not ask if you really mean it before getting rid of it forever(!), so this can be dangerous if you're not paying attention. 
  - As an example, let's use our **grep** command to pull out he seastar samples that started healthy and then became sick. But perhaps we later decide we're not going to work with those samples, so we use **rm** to delete that file:

```bash
[srkeller@pbio381 mydata]$ ll
total 8
drwxr-xr-x. 2 srkeller users   60 Jan 31 21:12 sample_by_disease
-rw-r--r--. 1 srkeller users  282 Feb  1 05:35 ssw_HSonly.txt
-rwxrwxr-x. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
[srkeller@pbio381 mydata]$ rm ssw_HSonly.txt 
[srkeller@pbio381 mydata]$ ll
total 4
drwxr-xr-x. 2 srkeller users   60 Jan 31 21:12 sample_by_disease
-rwxrwxr-x. 1 srkeller users 1255 Jan 31 17:42 ssw_samples.txt
[srkeller@pbio381 mydata]$
```
- Gone! Forever! If that worries you, you can change your personal settings so that the server asks you to confirm deletion before it acts. To do this, we'll need to follow a couple of new steps:


1.    **cd** to your home directory (~/)
      2. list all the files, including "hidden" ones that aren't usually shown. To do this, use `ll -a`.
         3. Look for a file called ".bashrc" — this contains your settings for how you interact with the server when you log in.
         4. We're going to open this file and edit it to add a setting to request that **rm** confirms deletion with us. To edit text files on the fly in UNIX, you can use the built-in text editor, "vim": `vim .bashrc`
         5. You should see something that looks like this:

```bash
  # .bashrc

  # Source global definitions
  if [ -f /etc/bashrc ]; then
          . /etc/bashrc
  fi

  # Uncomment the following line if you don't like systemctl's auto-paging feature:
  # export SYSTEMD_PAGER=

  # User specific aliases and functions

```

6.   Use your arrow key to move your cursor down to the last line, below ""# User specific aliases and functions" — this is where we're going to insert our new function.

7. By defauly, vim is in read-only mode when it opens files. To go into edit mode, press your "i" key (for "insert"). You are now able to make changes to the file.

8. Add the following text on a new line directly below the "# User specific…" line:

       `alias rm='rm -i'`

9. Your file should now look like this:

```bash
  # .bashrc

  # Source global definitions
  if [ -f /etc/bashrc ]; then
          . /etc/bashrc
  fi

  # Uncomment the following line if you don't like systemctl's auto-paging feature:
  # export SYSTEMD_PAGER=

  # User specific aliases and functions

  alias rm='rm -i'
```

10.    You're now ready to get out of edit mode (hit the `escape key`), save your changes (type `:w`), and exit vim (type `:q`).

11. These changes won't take effect until you log out (type `exit` to log out of the server). But from now on, every time you log in, the server will remember that you want a reminder before deleting any of your work.

      ​

##Let's review what we've learned so far…##

- Logging in to the server: `ssh netid@pbio381.uvm.edu`
- Finding what directory you're in: `pwd`
- Listing files in your current directory, or changing to a new directory: `ll`, `cd`
- Making a new folder: `mkdir foldername`
- Location of shared space, data, and programs on our class server:

```
[srkeller@pbio381 ~]$ cd /data/
[srkeller@pbio381 data]$ ll
total 8
drwxr-xr-x.  5 root root       73 Jan 31 17:35 archive
drwxrwxr-x.  2 root pb381adm   40 Nov 30  2015 packages
drwxrwxr-x. 33 root pb381adm 4096 Nov 30  2015 popgen
drwxrwxr-x.  3 root pb381adm   42 Jan 30 09:08 project_data
drwxrwxr-x.  2 root pb381adm    6 Oct  2  2015 scripts
drwxr-xr-x. 18 root root     4096 Sep  2  2015 users
[srkeller@pbio381 data]$ 
```

- Copying or moving files from one location to another: `cp filename destinationpath/` or `mv filename destinationpath/` 
- Peeking into the first or last few lines of a file: `head filename`, `tail filename`
- Searching within a file for a match: `grep 'search string' filename`
- Outputing the results of a command to a new file: `grep 'search string' filename >outputfilename`
- Using wildcards to work on multiple files at the same time: `mv *.txt ~/newfolder`
- Using the "pipe" to send the output of one command to the input of another: `grep 'INT' filename | wc `
- Removing files or folders: `rm`
- Editing text files on the server: `vim filename`       

******************

### Handy [UNIX cheat sheet](https://files.fosswire.com/2007/08/fwunixref.pdf) for helping to remember some of these commonly used commands (and others)

### Here's another useful [UNIX cheatsheet](http://cheatsheetworld.com/programming/unix-linux-cheat-sheet/)



