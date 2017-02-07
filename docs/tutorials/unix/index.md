<style src="../../includes/media/tute.css" ></style>
<style>em {font-style: normal; font-family: courier new;}</style>

# Introduction to Unix

A hands-on-workshop covering the basics of the Unix/Linux command line interface

## How to use this workshop

The workshop is broken up into a number of *Topics* each focusing on a particular aspect of Unix.  You should take a short break between 
each to refresh and relax before tackling the next.

*Topic*s may start with some background followed by a number of *exercises*.  Each *exercise* begins with a *question*, then 
sometimes a *hint* (or two) and finishes with the suggested *answer*.

### Question

An example question looks like:

\showable{What is the Answer to Life?}{question}

\endshowable

### Hint

Depending on how much of a challenge you like, you may choose to use hints.  Even if you work out the answer without hints, its a good 
idea to read the hints afterwards because they contain extra information that is good to know.

Note: *hint*s may be staged, that is, there may be a *more* section within a hint for further hints

\showable{Hint}{hint} &lt;- click here to reveal hint

What is the answer to everything?

As featured in "The Hitchhiker's Guide to the Galaxy"

\showable{More|Less} &lt;- and here to show more

It is probably a two digit number

\endshowable

\endshowable

### Answer

Once you have worked out the answer to the question expand the Answer section to check if you got it correct.

\showable{Answer}{answer} &lt;- click here to reveal answer

**Answer**: 42

Ref: [Number 42 (Wikipedia)](http://en.wikipedia.org/wiki/42_%28number%29)

\endshowable

### Usage Style

This workshop attempts to cater for two usage styles:

1. **Problem solver**: for those who like a challenge and learn best be trying to solve the problems by-them-selves (hints optional):
	* Attempt to answer the question by yourself.
	* Use hints when you get stuck.
	* Once solved, reveal the answer and read through our suggested solution.
	* Its a good idea to read the hints and answer description as they often contain extra useful information.
2. **By example**:  for those who learn by following examples:  [Expand](?exp) all sections
	* Expand the Answer section at the start of each question and follow along with the commands that are shown and check you get the
	  same (or similar) answers.
	* Its a good idea to read the hints and answer description as they often contain extra useful information.






## Topic 1: Remote log in

In this topic we will learn how to connect to a *Unix* computer via a method called *SSH* and run a few basic commands.


### Connecting to a Unix computer

To begin this workshop you will need to connect to an HPC.  Today we will use the LIMS-HPC.  The computer called 
*lims-hpc-m* (m is for master which is another name for head node) is the one that coordinates all the HPCs tasks.

**Server details**:

* **host**: lims-hpc-m.latrobe.edu.au
* **port**: 6022 
* **username**: trainingXX (where XX is a two digit number, provided at workshop)
* **password**: (provided at workshop) 

{!docs/includes/connecting.md!}



**Note**: for security reasons ssh will not display any characters when you enter your password. This 
can be confusing because it appears as if your typing is not recognised by the computer. Don’t be 
alarmed; type your password in and press return at the end.

LIMS-HPC is a high performance computer for La Trobe Users.  Logging in connects your local computer 
(e.g. laptop) to LIMS-HPC, and allows you to type commands into the Unix prompt which are run on 
the HPC, and have the results displayed on your local screen.

You will be allocated a training account on LIMS-HPC for the duration of the workshop. Your 
username and password will be supplied at the start of the workshop.

Log out of LIMS-HPC, and log back in again (to make sure you can repeat the process).

All the remaining parts assume that you are logged into LIMS-HPC over ssh.

### Exercises

\showable+{1.1) When you’ve logged into LIMS-HPC run the following commands and see what they do:}{question}

```sh
who
whoami
date
cal
hostname
/home/group/common/training/Intro_to_Unix/hi
```

\endshowable

\showable{Answer}{answer}

* **who**: displays a list of the users who are currently using this Unix computer.
* **whoami**: displays your username (i.e. they person currently logged in).
* **date**: displays the current date and time.
* **cal**: displays a calendar on the terminal.  It can be configured to display more than just 
the current month.
* **hostname**: displays the name of the computer we are logged in to.
* **/home/group/common/training/Intro_to_Unix/hi**: displays the text "Hello World"

\endshowable






## Topic 2: Exploring your home directory

In this topic we will learn how to "look" at the filesystem and further expand our repertoire of Unix commands. 

**Duration**: 20 minutes. <!-- See "The shell and the command line" and "The file system" section of the workshop notes. -->

**Relevant commands**: *ls*, *pwd*, *echo*, *man*

Your home directory contains your own private working space.  Your *current working directory* is automatically set 
to your *home* directory when you log into a Unix computer.

\showable{2.1) Use the *ls* command to list the files in your *home* directory.  How many files are there?}{question}\endshowable

\showable{Hint}{hint}

Literally, type *ls* and press the *ENTER* key.

\endshowable

\showable{Answer}{answer}

```sh
$ ls
exp01  file01  muscle.fq
```

When running the *ls* command with no options it will list files in your current working directory.  The place 
where you start when you first login is your *HOME* directory.

**Answer**: 3 (exp01, file01 and muscle.fq)

\endshowable

---

The above answer is not quite correct.  There are a number of *hidden* files in your home directory as well.

\showable{2.2) What *flag* might you use to display *all* files with the *ls* command?  How many files are really there?}{question}\endshowable

\showable{Hint}{hint}

Take the *all* quite literally.

\showable{More|Less}

Type *ls --all* and press the *ENTER* key.

\endshowable

\endshowable

\showable{Answer}{answer}

**Answer 1**: *--all* (or *-a*) flag

Now you should see several files in your home directory whose names all begin with a dot. All these files are 
created automatically for your user account. They are mostly configuration options for various programs including 
the shell. It is safe to ignore them for the moment.

```sh
$ ls --all
.              .bash_logout    exp01    .lesshst
..             .bash_profile   file01   muscle.fq
.bash_history  .bashrc         .kshrc   .viminfo
```

There are two trick files here; namely *.* and *..* which are not real files but instead, shortcuts.  *.* is a shortcut
for the current directory and *..* a shortcut for the directory above the current one.

**Answer 2**: 10 files (don't count *.* and *..*)

\endshowable

---

\showable{2.3) What is the full path name of your *home* directory?}{question}\endshowable

\showable{Hint}{hint}

Remember your *Current Working Directory* start's in your *home* directory (and the hint from the slides).

\showable{More|Less}

Try a shortened version of *print working directory*

\endshowable

\endshowable

\showable{Answer}{answer}

You can find out the full path name of the current working directory with the *pwd* command. Your home directory 
will look something like this:

```sh
$ pwd
/home/trainingXY
```

**Answer**: */home/trainingXY*

where *XY* is replaced by some 2 digit sequence.

**Alternate method**:
You can also find out the name of your home directory by printing the value of the *$HOME* shell variable:

```sh
echo $HOME
```

\endshowable

---

\showable{2.4) Run *ls* using the long flag (*-l*), how did the output change?}{question}\endshowable

\showable{Hint}{hint}

Run *ls -l*

\endshowable

\showable{Answer}{answer}

**Answer**: it changed the output to place 1 file/directory per line.  It also added some extra information
about each.

```sh
$ ls -l
total 32
drwxr-x--- 2 training01 training 2048 Jun 14 11:28 exp01
-rw-r----- 1 training01 training   97 Jun 14 11:28 file01
-rw-r----- 1 training01 training 2461 Jun 14 11:28 muscle.fq
```

**Details**:

```sh
drwxr-x--- 2 training01 training 2048 Jun 14 11:28 exp01
\--------/ ^ \--------/ \------/ \--/ \----------/ \---/
permission |  username   group   size    date       name
       /---^---\
       linkcount
```

Where:

* **permissions**: 4 parts, file type, user perms, group perms and other perms
	* *filetype*: 1 character, *d* = directory and *-* regular file
	* *user* permissions: 3 characters, *r* = read, *w* = write, *x* = execute and *-* no permission
	* *group* permissions: same as user except for users within the owner group
	* *other* permissions: same as user except for users that are not in either user *or* *group*
* **username**: the user who *owns* this file/directory
* **group**: the group name who *owns* this file/directory
* **size**: the number of bytes this file/directory takes to store on disk
* **date**: the date and time when this file/directory was *last edited*
* **name**: name of the file
* **linkcount**: technical detail which represents the number of links this file has in the file system (safe to ignore)

\endshowable

---

\showable{2.5) What type of file is *exp01* and *muscle.fq*?}{question}\endshowable

\showable{Hint}{hint}

Check the output from the *ls -l*.

\endshowable

\showable{Answer}{answer}

**Answer**:

* *exp01*: Directory (given the 'd' as the first letter of its permissions)
* *muscle.fq*: Regular File (given the '-')

\endshowable

---


\showable{2.6) Who has permission to *read*, *write* and *execute* your *home* directory?}{question}\endshowable

\showable{Hint}{hint}

You can also give *ls* a filename as the first option.

\showable{More|Less}

*ls -l* will show you the contents of the *CWD*; how might you see the contents of the *parent* directory? (remember
the slides)

\endshowable

\endshowable

\showable{Answer}{answer}

If you pass the *-l* flag to ls it will display a "long" listing of file information including file permissions.

There are various ways you could find out the permissions on your home directory.


**Method 1**: given we know the *CWD* is our home directory.

```sh
$ ls -l ..
...
drwxr-x--- 4 trainingXY training  512 Feb  9 14:18 trainingXY
...
```
The *..* refers to the parent directory.


**Method 2**: using $HOME.  This works no matter what our *CWD* is set to. 

You could list the permissions of all files and directories in the parent directory of your home:

```sh
$ ls -l $HOME/..
...
drwxr-x--- 4 trainingXY training  512 Feb  9 14:18 trainingXY
...
```

In this case we use the shell variable to refer to our home directory.


**Method 3**: using *~* (tilde) shortcut

You may also refer to your home directory using the *~* (tilde) character:

```sh
$ ls -l ~/..
...
drwxr-x--- 4 trainingXY training  512 Feb  9 14:18 trainingXY
...
```

All 3 of the methods above mean the same thing.

You will see a list of files and directories in the parent directory of your home directory. One of them will 
be the name of your home directory, something like *trainingXY*.  Where *XY* is replaced by a two digit string

**Altername**: using the *-a* flag and looking at the *.* (dot) special file.

```sh
$ ls -la
...
drwxr-x--- 4 trainingXY training  512 Feb  9 14:18 .
...
```


**Answer**: *drwxr-x---*

* **You**: read (see filenames), write (add, delete files), execute (change your CWD to this directory).
* **Training users**: read, execute
* **Everyone else**: No access

**Discussion on Permissions**:

The permission string is *"drwxr-x---"*. The *d* means it is a directory. The *rwx* means that the owner of the directory 
(your user account) can *read*, *write* and *execute* the directory. Execute permissions on a directory means that you 
can *cd* into the directory. The *r-x* means that anyone in the same user group as *training* can read or execute the 
directory. The *---* means that nobody else (other users on the system) can do anything with the directory.

\endshowable

---

<div class="info">
<p><b><em>man</em> is for manual</b>: and it will be your best friend!</p>
<p>Manual pages include a lot of detail about a command and its available flags/options.  It should be your first (or second) 
port of call when you are trying to work out what a command or option does.</p>
<p>You can scroll <em>up</em> and <em>down</em> in the man page using the <em>arrow</em> keys.</p>
<p>You can search in the man page using the forward 
slash followed by the search text followed by the <em>ENTER</em> key. e.g. 
type <em>/hello</em> and press <em>ENTER</em> to search for the word <em>hello</em>.  Press <em>n</em> key to find next 
occurance of <em>hello</em> etc.</p>
<p>You can <em>quit</em> the man page by pressing <em>q</em>.</p>
</div>

---

\showable{2.7) Use the *man* command to find out what the *-h* flag does for *ls*}{question}\endshowable

\showable{Hint}{hint}

Give *ls* as an option to *man* command.

\showable{More|Less}

*man ls*

\endshowable

\endshowable

\showable{Answer}{answer}

Use the following command to view the *man* page for *ls*:

```sh
$ man ls
```

**Answer**: You should discover that the *-h* option prints file sizes in human readable format

```sh
-h, --human-readable
              with -l, print sizes in human readable format (e.g., 1K 234M 2G)
```

\endshowable


---

\showable{2.8) Use the *-h*, how did the output change of *muscle.fq*?}{question}\endshowable

\showable{Hint}{hint}

Don't forget the *-l* option too.

\showable{More|Less}

Run *ls -lh*

\endshowable

\endshowable

\showable{Answer}{answer}

```sh
$ ls -lh
...
-rw-r----- 1 training01 training 2.5K Jun 14 11:28 muscle.fq
```

**Answer**: it changed the output so the *filesize* of *muscle.fq* is now *2.5K* instead of *2461*

\endshowable









## Topic 3: Exploring the file system

In this topic we will learn how to move around the filesystem and see what is there.

**Duration**: 30 minutes. <!-- See "The file system" section of the workshop notes. -->

**Relevant commands**: *pwd*, *cd*, *ls*, *file*

\showable{3.1) Print the value of your current working directory.}{question}\endshowable

\showable{Answer}{answer}

The *pwd* command prints the value of your current working directory.

```sh
$ pwd
/home/training01
```

\endshowable

---

\showable{3.2) List the contents of the root directory, called '*/*' (forward 
slash).}{question}\endshowable

\showable{Hint}{hint}

*ls* expects a single option which is the directory to change too.

\endshowable

\showable{Answer}{answer}

```sh
$ ls /
applications-merged  etc         media    root     tmp
bin                  home        mnt      sbin     usr
boot                 lib         oldhome  selinux  var
data                 lib64       opt      srv
dev                  lost+found  proc     sys
```

Here we see that *ls* can take a filepath as its argument, which allows you to list the contents of directories 
other than your current working directory.

\endshowable

---

\showable{3.3) Use the *cd* command to change your working directory to the root directory.  Did your prompt 
change?}{question}\endshowable

\showable{Hint}{hint}

*cd* expects a single option which is the directory to change to

\endshowable

\showable{Answer}{answer}

The *cd* command changes the value of your current working directory. To change to the root directory use the 
following command:

```sh
$ cd /
```

**Answer**: Yes, it now says the CWD is */* instead of *~*.

Some people imagine that changing the working directory is akin to moving your focus within the file system. 
So people often say "move to", "go to" or "charge directory to" when they want to change the working directory.

The root directory is special in Unix. It is the topmost directory in the whole file system.

\endshowable

---

<div class="info">
<b>Output on ERROR only</b>: Many Unix commands will not produce any output if everything went well; <em>cd</em> is one
such command.  However, it will get grumpy if something went wrong by way of an error message on-screen.
</div>

---

\showable{3.4) List the contents of the CWD and verify it matches the list in 3.2}{question}\endshowable

\showable{Hint}{hint}

*ls*

\endshowable

\showable{Answer}{answer}

Assuming you have changed to the root directory then this can be achieved with *ls*, or *ls -a* (for all files) or 
*ls -la* for a long listing of all files.

If you are not currently in the root directory then you can list its contents by passing it as an argument to ls:

```sh
$ ls
applications-merged  etc         media    root     tmp
bin                  home        mnt      sbin     usr
boot                 lib         oldhome  selinux  var
data                 lib64       opt      srv
dev                  lost+found  proc     sys
```

**Answer**: Yes, we got the same output as exercise 3.2



\endshowable

---

\showable{3.5) Change your current working directory back to your home directory. What is the simplest Unix command that 
will get you back to your home directory from anywhere else in the file system?}{question}\endshowable

\showable{Hint}{hint}

The answer to exercise 2.6 might give some hints on how to get back to the home directory

\showable{More|Less}

*$HOME*, *~*, */home/trainingXY* are all methods to name your home directory.  Yet there is a simpler method; the answer
is buried in *man cd* however *cd* doesn't its own manpage so you will need to search for it.

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *cd* command to change your working directory to your home directory. There are a number of ways to refer 
to your home directory:

```sh
cd $HOME
```

is equivalent to:

```sh
cd ~
```

The simplest way to change your current working directory to your home directory is to run the *cd* command with 
no arguments:

**Answer**: the simplest for is cd with NO options.

```sh
cd
```

This is a special-case behaviour which is built into *cd* for convenience.

\endshowable

---

\showable{3.6) Change your working directory to */home/group/common/training/Intro_to_Unix/*}{question}\endshowable

\showable{Answer}{answer}

```sh
cd /home/group/common/training/Intro_to_Unix/
```

\endshowable

---

\showable{3.7) List the contents of that directory. How many files does it contain?}{question}\endshowable

\showable{Hint}{hint}

*ls*

\endshowable

\showable{Answer}{answer}

You can do this with *ls*

```sh
$ ls
expectations.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

**Answer**: 7 files (expectations.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy)

\endshowable

---

\showable{3.8) What kind of *file* is */home/group/common/training/Intro_to_Unix/sleepy*?}{question}\endshowable

\showable{Hint}{hint}

Take the word *file* quite literally.

\showable{More|Less}

*file sleepy*

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *file* command to get extra information about the contents of a file:

Assuming your current working directory is */home/group/common/training/Intro_to_Unix/*

```sh
$ file sleepy
Bourne-Again shell script text executable
```

Otherwise specify the full path of sleepy:

```sh
$ file /home/group/common/training/Intro_to_Unix/sleepy
Bourne-Again shell script text executable
```

**Answer**: Bourne-Again shell script text executable

The "Bourne-Again shell" is more commonly known as BASH. The *file* command is telling us that sleepy 
is (probably) a shell script written in the language of BASH.

The file command uses various heuristics to guess the "type" of a file. If you want to know how it works 
then read the Unix manual page like so:

```sh
man file
```

\endshowable

---

\showable{3.9) What kind of *file* is */home/group/common/training/Intro_to_Unix/hi*?}{question}\endshowable

\showable{Hint}{hint}

Take the word *file* quite literally.

\endshowable

\showable{Answer}{answer}

Use the file command again. If you are in the same directory as *hi* then:

```sh
$ file hi
ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked (uses shared libs), for GNU/Linux 
2.6.9, not stripped
```

**Answer**: ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked (uses shared libs), for GNU/Linux 

This rather complicated output is roughly saying that the file called *hi* contains a binary executable 
program (raw instructions that the computer can execute directly).

\endshowable

---

\showable{3.10) What are the file permissions of */home/group/common/training/Intro_to_Unix/sleepy*? 
What do they mean?}{question}\endshowable

\showable{Hint}{hint}

Remember the *ls* command, and don't forget the *-l* flag

\endshowable

\showable{Answer}{answer}

You can find the permissions of *sleepy* using the *ls* command with the *-l* flag. If you are in the same 
directory as *sleepy* then:

```sh
$ ls -l sleepy
-rw-r--r-- 1 arobinson common 183 Feb  9 16:36 sleepy
```

**Answer**: We can see that this particular instance of sleepy is owned by the user arobinson, and is part of the common 
user group. It is 183 bytes in size, and was last modified on the 9th of February at 4:36pm. The file is 
readable to everyone, and writeable only to training01.  The digit '1' between the file permission string and 
the owner indicates that there is one link to the file. The Unix file system allows files to be referred to 
by multiple "links". When you create a file it is referred to by one link, but you may add others later. For 
future reference: links are created with the *ln* command.

\endshowable

---

\showable{3.11) Change your working directory back to your home directory ready for the next topic.}{question}\endshowable

\showable{Hint}{hint}

*cd*

\endshowable

\showable{Answer}{answer}

You should know how to do this with the cd command:

```sh
cd
```

\endshowable









## Topic 4: Working with files and directories

In this topic we will start to read, create, edit and delete files and directories.

**Duration**: 50 minutes.  <!-- See "Working with files" from the workshop notes. -->

**Relevant commands**: *mkdir*, *cp*, *ls*, *diff*, *wc*, *nano*, *mv*, *rm*, *rmdir*, *head*, *tail*, *grep*, *gzip*, *gunzip*

<div class="info">
<b>Hint</b>: Look at the commands above; you will need them roughly in order for this topic.  Use the <em>man</em>
command find out what they do, in particular the NAME, SYNOPSIS and DESCRIPTION sections.
</div>

---


\showable{4.1) In your home directory make a sub-directory called test.}{question}\endshowable

\showable{Hint}{hint}

You are trying to *make a directory*, which of the above commands looks like a shortened version of this?

\showable{More|Less}

*mkdir*

\endshowable

\endshowable

\showable{Answer}{answer}

Make sure you are in your home directory first. If not *cd* to your home directory.

Use the *mkdir* command to make new directories:

```sh
$ mkdir test
```

Use the *ls* command to check that the new directory was created.

```sh
$ ls
exp01  file01  muscle.fq  test
```

\endshowable

---

\showable{4.2) Copy all the files from */home/group/common/training/Intro_to_Unix/* into the newly created 
test directory.}{question}\endshowable

\showable{Hint}{hint}

You are trying to *copy*, which of the above commands looks like a shortened version of this?

\showable{More|Less}

```sh
$ man cp
...
SYNOPSIS
       cp [OPTION]... [-T] SOURCE DEST
...
DESCRIPTION
       Copy SOURCE to DEST, or multiple SOURCE(s) to DIRECTORY.
```

which means *cp* expects zero or more flags, a SOURCE file followed by a DEST file or directory 

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *cp* command to copy files. 

<div class="info"><b>Wildcards</b>: You could copy them one-by-one, but that would be tedious, so use 
the <em>*</em> wildcard to specify that you want to copy all the files.
</div>

There are a number of ways you could do this depending on how you specify the source and destination 
paths to *cp*. You only need to perform one of these ways, but we show multiple ones for your reference.

**Answer 1**: From your home directory:

```sh
$ cp /home/group/common/training/Intro_to_Unix/* test
```

**Answer 2**: Change to the test directory and then copy (assuming you started in your home directory):

```sh
$ cd test
$ cp /home/group/common/training/Intro_to_Unix/* .
```

In the example above the '*.*' (dot) character refers to the current working directory. It should be 
the test subdirectory of your home directory.

**Answer 3**: Change to the /home/group/common/training/Intro_to_Unix/ directory and then copy:

```sh
cd /home/group/common/training/Intro_to_Unix/
cp * ~/test
```

Remember that ~ is a shortcut reference to your home directory.

\endshowable

---

**Note**: This exercise assumes that the copy command from the previous exercise was successful. 

\showable{4.3) Check that the file size of *expectations.txt* is the same in both the directory that you copied 
it from and the directory that you copied it to.}{question}\endshowable

\showable{Hint}{hint}

Remember *ls* can show you the file size (with one of its flags)

\showable{More|Less}

*ls -l*

\endshowable

\endshowable

\showable{Answer}{answer}

Use *ls -l* to check the size of files.

You could do this in many ways depending on the value of your working directory. We just show one possible 
way for each file:

```sh
$ ls -l /home/group/common/training/Intro_to_Unix/expectations.txt

$ ls -l ~/test/expectations.txt
```

From the output of the above commands you should be able to see the size of each file and check that they 
are the same. 

**Answer**: They should each be *1033773* bytes

**Alternate**: Sometimes it is useful to get file sizes reported in more "human friendly" units than bytes. If this is 
true then try the *-h* option for ls:

```sh
$ ls -lh /home/group/common/training/Intro_to_Unix/expectations.txt
-rw-r--r-- 1 arobinson common 1010K Mar 26  2012 /home/group/common/training/Intro_to_Unix/expectations.txt
```

In this case the size is reported in kilobytes as *1010K*. Larger files are reported in megabytes, gigabytes 
etcetera.

\endshowable

---

**Note**: this exercise assumes your working directory is *~/test*; if not run *cd ~/test*

\showable{4.4) Check that the contents of expectations.txt are the same in both the directory that you copied 
it from and the directory that you copied it to.}{question}\endshowable

\showable{Hint}{hint}

What is the opposite of *same*?

\showable{More|Less}

*diff*erence

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *diff* command to compare the contents of two files.

```sh
$ diff /home/group/common/training/Intro_to_Unix/expectations.txt expectations.txt
```

If the two files are identical the *diff* command will NOT produce any output)

**Answer**: Yes, they are the same since no output was given.

\endshowable

---

\showable{4.5) How many lines, words and characters are in expectations.txt?}{question}\endshowable

\showable{Hint}{hint}

Initialisms are key

\showable{More|Less}

*w*ord *c*ount

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *wc* (for "word count") to count the number of characters, lines and words in a file:

```sh
$ wc expectations.txt
  20415  187465 1033773 expectations.txt
```

**Answer**: There are *20415* lines, *187465* words and *1033773* characters in expectations.txt.

To get just the line, word or character count:

```sh
$ wc -l expectations.txt
20415 expectations.txt
$ wc -w expectations.txt
187465 expectations.txt
$ wc -c expectations.txt
1033773 expectations.txt
```

\endshowable

---

\showable{4.6) Open *~/test/expectations.txt* in the *nano* text editor, delete the first line of text, and 
save your changes to the file. Exit *nano*.}{question}\endshowable

\showable{Hint}{hint}

*nano FILENAME*

Once *nano* is open it displays some command hints along the bottom of the screen.

\showable{More|Less}

*^O* means hold the *Control* (or CTRL) key while pressing the *o*.  Dispite what it displays, you need to type 
the lower-case letter that follows the *^* character.

WriteOut is another name for Save.

\endshowable

\endshowable

\showable{Answer}{answer}

Take some time to play around with the *nano* text editor.

*Nano* is a very simple text editor which is easy to use but limited in features. More powerful 
editors exist such as *vim* and *emacs*, however they take a substantial amount of time to learn.

\endshowable

---

\showable{4.7) Did the changes you made to *~/test/expectations.txt* have any effect on 
*/home/group/common/training/Intro_to_Unix/expectations.txt*? How can you tell if two files are the 
same or different in their contents?}{question}\endshowable

\showable{Hint}{hint}

Remember exercise 4.4

\showable{More|Less}

Use *diff*

\endshowable

\endshowable

\showable{Answer}{answer}

Use *diff* to check that the two files are different after you have made the change to the copy of 
*expectations.txt* in your *~/test* directory.

```sh
diff ~/test/expectations.txt \
/home/group/common/training/Intro_to_Unix/expectations.txt
```

You could also use *ls* to check that the files have different sizes.

\endshowable

---

\showable{4.8) In your *test* subdirectory, rename *expectations.txt* to *foo.txt*.}{question}\endshowable

\showable{Hint}{hint}

Another way to think of it is *moving* it from *expectations.txt* to *foo.txt*

\showable{More|Less}

*mv*

Use *man mv* if you need to work out how to use it.

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *mv* command to rename the file:

```sh
$ mv expectations.txt foo.txt
$ ls
foo.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

\endshowable

---

\showable{4.9) Rename foo.txt back to expectations.txt.}{question}\endshowable

\showable{Answer}{answer}

Use the *mv* command to rename the file:

```sh
$ mv foo.txt expectations.txt
$ ls
expectations.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

Use *ls* to check that the file is in fact renamed.

\endshowable

---

\showable{4.10) Remove the file *expectations.txt* from your *test* directory.}{question}\endshowable

\showable{Hint}{hint}

We are trying to *remove* a file, check the commands at the top of this topic.

\showable{More|Less}

*rm*

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *rm* command to remove files (carefully):

```sh
$ rm expectations.txt
$ ls
hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

\endshowable

---

\showable{4.11) Remove the entire *test* directory and all the files within it.}{question}\endshowable

\showable{Hint}{hint}

We are trying to *remove a directory*.

\showable{More|Less}

You could use *rmdir* but there is an easier way using just *rm* and a flag.

\endshowable

\endshowable

\showable{Answer}{answer}

You could use the *rm* command to remove each file individually, and then use the *rmdir* command 
to remove the directory. Note that *rmdir* will only remove directories that are empty (i.e. do not 
contain files or subdirectories).

A faster way is to pass the *-r* (for recursive) flag to *rm* to remove all the files and the 
directory in one go:

**Logical Answer**:
```sh
cd ~
rm test/*
rmdir test
```

**Easier Answer**:
```sh
cd ~
rm -r test
```

<div class="error"><b>Warning</b>: Be very careful with <em>rm -r</em>, it will remove all files 
and all subdirectories underneath the specified directory. This could be catastrophic if you do it 
in the wrong location! Now is a good moment to pause and think about file backup strategies.</div>

\endshowable

---

\showable{4.12) Recreate the test directory in your home directory and copy all the files from 
*/home/group/common/training/Intro_to_Unix/* back into the test directory.}{question}\endshowable

\showable{Hint}{hint}

See exercises 4.1 and 4.2

\endshowable

\showable{Answer}{answer}

Repeat exercises 4.1 and 4.2.

```sh
$ cd ~
$ mkdir test
$ cp /home/group/common/training/Intro_to_Unix/* test
```

\endshowable

---

\showable{4.13) Change directories to *~/test* and use the *cat* command to display the entire contents 
of the file *hello.c*}{question}\endshowable

\showable{Hint}{hint}

Use *man* if you can't guess how it might work.

\endshowable

\showable{Answer}{answer}

```sh
$ cd ~/test
$ cat hello.c
#include <stdio.h>
int main(void) {
    printf ("Hello World\n");
    return 0;
}
```

*hello.c* contains the source code of a C program. The compiled executable version of this code 
is in the file called *hi*, which you can run like so:

```sh
$ ./hi
Hello World
```

\endshowable

---

\showable{4.14) Use the *head* command to view the first *20* lines of the file *sample_1.fastq*}{question}\endshowable

\showable{Hint}{hint}

Remember your *best* friend!

\showable{More|Less}

Use *man* to find out what option you need to add to display a given number of *lines*.

\endshowable

\endshowable

\showable{Answer}{answer}

```sh
$ head -20 sample_1.fastq
@IRIS:7:1:17:394#0/1
GTCAGGACAAGAAAGACAANTCCAATTNACATTATG
+IRIS:7:1:17:394#0/1
aaabaa`]baaaaa_aab]D^^`b`aYDW]abaa`^
@IRIS:7:1:17:800#0/1
GGAAACACTACTTAGGCTTATAAGATCNGGTTGCGG
+IRIS:7:1:17:800#0/1
ababbaaabaaaaa`]`ba`]`aaaaYD\\_a``XT
@IRIS:7:1:17:1757#0/1
TTTTCTCGACGATTTCCACTCCTGGTCNACGAATCC
+IRIS:7:1:17:1757#0/1
aaaaaa``aaa`aaaa_^a```]][Z[DY^XYV^_Y
@IRIS:7:1:17:1479#0/1
CATATTGTAGGGTGGATCTCGAAAGATATGAAAGAT
+IRIS:7:1:17:1479#0/1
abaaaaa`a```^aaaaa`_]aaa`aaa__a_X]``
@IRIS:7:1:17:150#0/1
TGATGTACTATGCATATGAACTTGTATGCAAAGTGG
+IRIS:7:1:17:150#0/1
abaabaa`aaaaaaa^ba_]]aaa^aaaaa_^][aa
```

\endshowable

---

\showable{4.15) Use the *tail* command to view the last *8* lines of the file *sample_1.fastq*}{question}\endshowable

\showable{Hint}{hint}

Its very much like *head*.

\endshowable

\showable{Answer}{answer}

```sh
tail -8 sample_1.fastq
@IRIS:7:32:731:717#0/1
TAATAATTGGAGCCAAATCATGAATCAAAGGACATA
+IRIS:7:32:731:717#0/1
ababbababbab]abbaa`babaaabbb`bbbabbb
@IRIS:7:32:731:1228#0/1
CTGATGCCGAGGCACGCCGTTAGGCGCGTGCTGCAG
+IRIS:7:32:731:1228#0/1
`aaaaa``aaa`a``a`^a`a`a_[a_a`a`aa`__
```

\endshowable

---

\showable{4.16) Use the *grep* command to find out all the lines in *moby.txt* that contain the word 
"Ahab"}{question}\endshowable

\showable{Hint}{hint}

One might say we are 'looking for the *pattern* "Ahab"'

\showable{More|Less}

```sh
$ man grep
...
SYNOPSIS
       grep [OPTIONS] PATTERN [FILE...]
...
```

\endshowable

\endshowable

\showable{Answer}{answer}

```sh
$ grep Ahab moby.txt
"Want to see what whaling is, eh? Have ye clapped eye on Captain Ahab?"
"Who is Captain Ahab, sir?"
"Aye, aye, I thought so. Captain Ahab is the Captain of this ship."
... AND MUCH MUCH MORE ...
```

If you want to know how many lines are in the output of the above command you can "pipe" it 
into the *wc -l* command:

```sh
$ grep Ahab moby.txt | wc -l
491
```

which shows that there are *491* lines in *moby.txt* that contain the word Ahab.

\endshowable

---

\showable{4.17) Use the *grep* command to find out all the lines in *expectations.txt* that contain the 
word "the" with a case insensitive search (it should count "the" "The" "THE" "tHe" etcetera)
.}{question}\endshowable

\showable{Hint}{hint}

One might say we are *ignoring case*.

\showable{More|Less}

```sh
$ man grep
...
       -i, --ignore-case
              Ignore case distinctions in both the PATTERN and the input files.  (-i is specified by POSIX.)
...
```

\endshowable

\endshowable

\showable{Answer}{answer}

Use the *-i* flag to *grep* to make it perform case insensitive search:

```sh
$ grep -i the expectations.txt
The Project Gutenberg EBook of Great Expectations, by Charles Dickens
This eBook is for the use of anyone anywhere at no cost and with
re-use it under the terms of the Project Gutenberg License included
[Project Gutenberg Editor's Note: There is also another version of
... AND MUCH MUCH MORE ...
```

Again, "pipe" the output to *wc -l* to count the number of lines:

```sh
$ grep -i the expectations.txt  | wc -l
8165
```

\endshowable

---

\showable{4.18) Use the *gzip* command to compress the file *sample_1.fastq*. Use *gunzip* to decompress it 
back to the original contents.}{question}\endshowable

\showable{Hint}{hint}

Use the above commands along with *man* and *ls* to see what happens to the file.

\endshowable

\showable{Answer}{answer}

Check the file size of sample_1.fastq before compressing it:

```sh
# check filesize
$ ls -l sample_1.fastq
-rw-r--r-- 1 training01 training 90849644 Jun 14 20:03 sample_1.fastq

# compress it (takes a few seconds)
$ gzip sample_1.fastq

# check filesize (Note: its name changed)
$ ls -l sample_1.fastq.gz
-rw-r--r-- 1 training01 training 26997595 Jun 14 20:03 sample_1.fastq.gz

# decompress it
$ gunzip sample_1.fastq.gz

$ ls -l sample_1.fastq
-rw-r--r-- 1 training01 training 90849644 Jun 14 20:03 sample_1.fastq
```

You will see that when it was compressed it is *26997595* bytes in size, making it about *0.3* times the size of the 
original file.

**Note**: in the above section the lines starting with *#* are comments so don't need to be copied but if you
do then they wont do anything.

\endshowable





## Topic 5: Pipes, output redirection and shell scripts

In this section we will cover a lot of the more advanced Unix concepts; it is here where you will start to see
the power of Unix.  I say *start* because this is only the "tip of the iceberg".

**Duration**: 50 minutes. <!-- See "Processes" from the workshop notes. -->

**Relevant commands**: *wc*, *paste*, *grep*, *sort*, *uniq*, *nano*, *cut*




\showable{5.1) How many *reads* are contained in the file *sample_1.fastq*?}{question}\endshowable

\showable{Hint}{hint}

Examine some of the file to work out how many lines each *read* takes up.

\showable{More|Less}

Count the number of lines

\endshowable

\endshowable

\showable{Answer}{answer}

We can answer this question by counting the number of lines in the file and dividing by 4:

```sh
$ wc -l sample_1.fastq
3000000
```

**Answer**: There are *3000000* lines in the file representing *750000* reads.

If you want to do simple arithmetic at the command line then you can use the "basic calculator" 
called *bc*:

```sh
$ echo "3000000 / 4" | bc
750000
```

<div class="info"><b>Note</b>: that the vertical bar character "|" is the Unix pipe (and is often 
called the "pipe symbol"). It is used for connecting the output of one command into the input of 
another command. We'll see more examples soon.</div>

*bc* is suitable for small calculations, but it becomes cumbersome for more complex examples. If 
you want to do more sophisticated calculations then we recommend to use a more general purpose 
programming language (such as Python etcetera).

\endshowable

---

\showable{5.2) How many reads in *sample_1.fastq* contain the sequence *GATTACA*?}{question}\endshowable

\showable{Hint}{hint}

Check out exercise 4.16

\endshowable

\showable{Answer}{answer}

Use *grep* to find all the lines that contain *GATTACA* and "pipe" the output to *wc -l* to count them:

```sh
$ grep GATTACA sample_1.fastq | wc -l
1119
```

**Answer**: *1119*

If you are unsure about the possibility of upper and lower case characters then consider using 
the *-i* (ignore case option for grep).

\endshowable

---

\showable{5.3) On what line numbers do the sequences containing *GATTACA* occur?}{question}\endshowable

\showable{Hint}{hint}

We are looking for the *line numbers*.

\showable{More|Less}

Check out the manpage for *grep* and/or *nl*

\endshowable

\endshowable

\showable{Answer}{answer}

You can use the *-n* flag to grep to make it prefix each line with a line number:

**Answer 1**: 
```sh
$ grep -n GATTACA sample_1.fastq
5078:AGGAAGATTACAACTCCAAGACACCAAACAAATTCC
7170:AACTACAAAGGTCAGGATTACAAGCTCTTGCCCTTC
8238:ATAGTTTTTTCGATTACATGGATTATATCTGTTTGC
... AND MUCH MUCH MORE ...
```

**Answer 2**: Or you can use the *nl* command to number each line of sample_1.fastq and then search for *GATTACA* 
in the numbered lines:

```sh
$ nl sample_1.fastq | grep GATTACA
  5078	AGGAAGATTACAACTCCAAGACACCAAACAAATTCC
  7170	AACTACAAAGGTCAGGATTACAAGCTCTTGCCCTTC
  8238	ATAGTTTTTTCGATTACATGGATTATATCTGTTTGC
... AND MUCH MUCH MORE ...
```

**Just the line numbers**:

If you just want to see the line numbers then you can "pipe" the output of the above command into 
*cut -f 1*:

```sh
$ nl sample_1.fastq | grep GATTACA | cut -f 1
  5078
  7170
  8238
... AND MUCH MUCH MORE ...
```

*cut* will remove certain columns from the input; in this case it will remove all except column 1
(a.k.a. field 1, hence the *-f 1* option)

```sh
$ grep -n GATTACA sample_1.fastq | cut -d: -f 1
5078
7170
8238
... AND MUCH MUCH MORE ...
```

\endshowable

---

\showable{5.4) Use the *nl* command to print each line of *sample_1.fastq* with its corresponding line 
number at the beginning.}{question}\endshowable

\showable{Hint}{hint}

Check answer to 5.3.

\endshowable

\showable{Answer}{answer}

```sh
$ nl sample_1.fastq
     1	@IRIS:7:1:17:394#0/1
     2	GTCAGGACAAGAAAGACAANTCCAATTNACATTATG
     3	+IRIS:7:1:17:394#0/1
     4	aaabaa`]baaaaa_aab]D^^`b`aYDW]abaa`^
     5	@IRIS:7:1:17:800#0/1
     6	GGAAACACTACTTAGGCTTATAAGATCNGGTTGCGG
     7	+IRIS:7:1:17:800#0/1
     8	ababbaaabaaaaa`]`ba`]`aaaaYD\\_a``XT
... AND MUCH MUCH MORE ...
```

There are a lot of lines in that file so this command might take a while to print all its output. 
If you get tired of looking at the output you can kill the command with *control-c* (hold the 
*control* key down and simultaneously press the "*c*" character).

\endshowable

---

\showable{5.5) Redirect the output of the previous command to a file called *sample_1.fastq.nl*. Check 
the first *20* lines of *sample_1.fastq.nl* with the *head* command. Use the *less* command to 
interactively view the contents of *sample_1.fastq.nl* (use the arrow keys to navigate up and down, 
*q* to quit and '*/*' to search). Use the search facility in less to find occurrences of 
*GATTACA*.}{question}\endshowable

\showable{Hint}{hint}

Ok that one was tough, *> FILENAME* is how you do it if you didn't break out an internet search for 
"redirect the output in Unix"

\endshowable

\showable{Answer}{answer}

```sh
$ nl sample_1.fastq > sample_1.fastq.nl
```

The greater-than sign "*>*" is the file redirection operator. It causes the standard output of the 
command on the left-hand-side to be written to the file on the right-hand-side.

You should notice that the above command is much faster than printing the output to the screen. 
This is because writing to disk can be performed much more quickly than rendering the output on 
a terminal.

To check that the first 20 lines of the file look reasonable you can use the *head* command like so:

```sh
$ head -20 sample_1.fastq.nl
     1	@IRIS:7:1:17:394#0/1
     2	GTCAGGACAAGAAAGACAANTCCAATTNACATTATG
     3	+IRIS:7:1:17:394#0/1
     4	aaabaa`]baaaaa_aab]D^^`b`aYDW]abaa`^
     5	@IRIS:7:1:17:800#0/1
     6	GGAAACACTACTTAGGCTTATAAGATCNGGTTGCGG
     7	+IRIS:7:1:17:800#0/1
     8	ababbaaabaaaaa`]`ba`]`aaaaYD\\_a``XT
...
```

The *less* command allows you to interactively view a file. The arrow keys move the page up and 
down. You can search using the '*/*' followed by the search term. You can quit by pressing "*q*". Note 
that the *less* command is used by default to display man pages.

```sh
$ less sample_1.fastq.nl
```

\endshowable

---

\showable+{5.6) The four-lines-per-read format of FASTQ is cumbersome to deal with. Often it would be 
preferable if we could convert it to tab-separated-value (TSV) format, such that each read appears 
on a single line with each of its fields separated by tabs. Use the following command to convert 
sample_1.fastq.nl into TSV format:}{question}

```sh
$ cat sample_1.fastq | paste - - - - > sample_1.tsv
```
\endshowable

\showable{Answer}{answer}

The *'-'* (dash) character has a special meaning when used in place of a file; it means use the standard
input instead of a real file.  Note: while it is fairly common in most Unix programs, not all wil support it.

The *paste* command is useful for merging multiple files together line-by-line, such that the *Nth* 
line from each file is joined together into one line in the output, separated by default with a 
*tab* character. In the above example we give paste 4 copies of the contents of *sample_1.fastq*, 
which causes it to join consecutive groups of 4 lines from the file into one line of output.

\endshowable

---

\showable+{5.7) Do you expect the output of the following command to produce the same output as above? and why?}{question}

```sh
$ paste sample_1.fastq sample_1.fastq sample_1.fastq sample_1.fastq > sample_1b.tsv
```

Try it, see what ends up in sample_1b.tsv (maybe use *less*)

\endshowable

\showable{Hint}{hint}

Use *less* to examine it.

\endshowable

\showable{Answer}{answer}

**Answer**: No, in the second instance we get 4 copies of each line.

**Why**: In the first command *paste* will use the input file (standard input) 4 times since the *cat* 
command will only give one copy of the file to *paste*, where as, in the second command *paste* will open 
the file 4 times.  Note: this is quite confusing and is not necessory to remember; its just an interesting
side point.

\endshowable

---

\showable{5.8) Check that *sample_1.tsv* has the correct number of lines. Use the *head* command to view 
the first *20* lines of the file.}{question}\endshowable

\showable{Hint}{hint}

Remember the *wc* command.

\endshowable

\showable{Answer}{answer}

We can count the number of lines in *sample_1.tsv* using *wc*:

```sh
$ wc -l sample_1.tsv
```

The output should be *750000* as expected (1/4 of the number of lines in sample_1.fastq).

To view the first *20* lines of *sample_1.tsv* use the *head* command:

```sh
$ head -20 sample_1.tsv
```

\endshowable

---

\showable{5.9) Use the *cut* command to print out the second column of *sample_1.tsv*. Redirect the 
output to a file called *sample_1.dna.txt*.}{question}\endshowable

\showable{Hint}{hint}

See exercise 5.3 (for cut) and 5.5 (redirection)

\endshowable

\showable{Answer}{answer}

The file sample_1.tsv is in column format. The cut command can be used to select certain columns 
from the file. The DNA sequences appear in column 2, we select that column using the -f 2 flag 
(the f stands for "field").

```sh
cut -f 2 sample_1.tsv > sample_1.dna.txt
```

Check that the output file looks reasonable using *head* or *less*.

\endshowable

---

\showable{5.10) Use the *sort* command to sort the lines of *sample_1.dna.txt* and redirect the output to 
*sample_1.dna.sorted.txt*. Use *head* to look at the first few lines of the output file. You should 
see a lot of repeated sequences of As.}{question}\endshowable

\showable{Hint}{hint}

Use *man* (sort) and see exercise 5.5 (redirection)

\endshowable

\showable{Answer}{answer}

```sh
$ sort sample_1.dna.txt > sample_1.dna.sorted.txt
```

Running *head* on the output file reveals that there are duplicate DNA sequences in the input FASTQ 
file.

\endshowable

---

\showable{5.11) Use the *uniq* command to remove duplicate consecutive lines from *sample_1.dna.sorted.txt*, 
redirect the result to *sample_1.dna.uniq.txt*. Compare the number of lines in sample1_dna.txt to 
the number of lines in *sample_1.dna.uniq.txt*.}{question}\endshowable

\showable{Hint}{hint}

I am pretty sure you have already used *man* (or just guessed how to use *uniq*).  You're also a gun at 
redirection now.

\endshowable

\showable{Answer}{answer}

```sh
$ uniq sample_1.dna.sorted.txt > sample_1.dna.uniq.txt
```

Compare the outputs of:

```sh
$ wc -l sample_1.dna.sorted.txt
750000
$ wc -l sample_1.dna.uniq.txt
614490
```

View the contents of *sample_1.dna.uniq.txt* to check that the duplicate DNA sequences have been 
removed.

\endshowable

---

\showable{5.12) Can you modify the command from above to produce *only* those sequences of DNA which were 
duplicated in *sample_1.dna.sorted.txt*?}{question}\endshowable

\showable{Hint}{hint}

Checkout the *uniq* manpage

\endshowable

\showable{Hint}{hint}

Look at the man page for uniq.

\endshowable

\showable{Answer}{answer}

Use the *-d* flag to *uniq* to print out only the duplicated lines from the file:

```sh
$ uniq -d sample_1.dna.sorted.txt > sample_1.dna.dup.txt
```

\endshowable

---

\showable{5.13) Write a *shell pipeline* which will print the number of duplicated DNA sequences in 
sample_1.fastq.}{question}\endshowable

\showable{Hint}{hint}

That is, *piping* most of the commands you used above instead of redirecting to file

\showable{More|Less}

I.e. 6 commands (*cat*, *paste*, *cut*, *sort*, *uniq*, *wc*)

\endshowable

\endshowable

\showable{Answer}{answer}

Finally we can 'pipe' all the pieces together into a sophisticated pipeline which starts with a 
FASTQ file and ends with a list of duplicated DNA sequences:

**Answer**:
```sh
$ cat sample_1.fastq | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
56079
```

The output file should have *56079* lines.

\endshowable

---

\showable{5.14) (Advanced) Write a shell script which will print the number of duplicated DNA sequences 
in sample_1.fastq.}{question}\endshowable

\showable{Hint}{hint}

Check out the *sleepy* file (with *cat* or *nano*); there is a bit of magic on the first line that you will need. 

You also need to tell bash that this file can be executed (check out *chmod* command).

\endshowable

\showable{Answer}{answer}

Put the answer to *5.13* into a file called *sample_1_dups.sh* (or whatever you want). Use *nano* to 
create the file. 

**Answer**: the contents of the file will look like this:

```sh
#!/bin/bash

cat sample_1.fastq | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
```

<div class="info"><b>Note</b>: the first line has special meaning.  If it starts with '<em>#!</em>' (Hash 
then exclamation mark) then it tells bash this file is a script that can be interpreted.  The command 
(including full path) used to intepret the script is placed right after the magic code.</div>

Give everyone execute permissions on the file with chmod:

```sh
$ chmod +x sample_1_dups.sh 
```

You can run the script like so:

```sh
$ ./sample_1_dups.sh
```

If all goes well the script should behave in exactly the same way as the answer to 5.13.

\endshowable

---

\showable{5.15) (Advanced) Modify your shell script so that it accepts the name of the input FASTQ file 
as a command line parameter.}{question}\endshowable

\showable{Hint}{hint}

Shell scripts can refer to command line arguments by their position using special variables called 
*$0*, *$1*, *$2* and so on. 

\showable{More|Less}

*$0* refers to the name of the script as it was called on the command line. 
*$1* refers to the first command line argument, and so on.

\endshowable

\endshowable


\showable{Answer}{answer}

Copy the shell script from *5.14* into a new file:

```sh
$ cp sample_1_dups.sh fastq_dups.sh
```

Edit the new shell script file and change it to use the command line parameters:

```sh
#!/bin/bash

cat $1 | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
```

You can run the new script like so:

```sh
$ ./fastq_dups.sh sample_1.fastq
```

In the above example the script takes *sample_1.fastq* as input and prints the number of duplicated 
sequences as output.

**A better Answer**:

Ideally we would write our shell script to be more robust. At the moment it just assumes there 
will be at least one command line argument. However, it would be better to check and produce an 
error message if insufficient arguments were given:

```sh
#!/bin/bash
if [ $# -eq 1 ]; then
    cat $1 | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
else
    echo "Usage: $0 <fastq_filename>"
    exit 1
fi
```

The '*if ...; then*' line means: do the following line(s) ONLY if the *...* (called condition) bit is true.

The '*else*' line means: otherwise do the following line(s) instead.  Note: it is optional.

The '*fi*' line means: this marks the end of the current *if* or *else* section.

The '*[ $# -eq 1 ]*' part is the condition:
  
* *$#*: is a special shell variable that indicates how many command line arguments were given. 
* *-eq*: checks if the numbers on either side if it are equal.
* *1*: is a number one

<div class="warning"><b>Spaces in conditions</b>:
Bash is VERY picky about the spaces within the conditions; if you get it wrong it will just behave strangely 
(without warning).  You MUST put a space near the share brackets and between each part of the condition!</div>

So in words our script is saying "if user provided 1 filename, then count the duplicates, otherwise print an error".

<div class="info"><b>Exit-status</b>:
It is a Unix standard that when the user provides incorrect commandline arguments we print a usage message 
and return a *non-zero* exit status.  The *exit status* is a standard way for other programs to know if
our program ran correctly; 0 means everything went as expected, any other number is an error.  If you don't
provide an *exit ..* line then it automatically returns a 0 for you.</div>

\endshowable

---

\showable{5.16) (Advanced) Modify your shell script so that it accepts zero or more FASTQ files on the 
command line argument and outputs the number of duplicated DNA sequences in each file.}{question}\endshowable

\showable{Answer}{answer}

We can add a loop to our script to accept multiple input FASTQ files:

```sh
#!/bin/bash
for file in $@; do
    dups=$(cat $file | paste - - - - | cut -f 2 | sort | uniq -d | wc -l)
    echo "$file $dups"
done
```

There's a lot going on in this script.

The *$@* is a sequence of all command line arguments.

The '*for ...; do*' (a.k.a. for loop) iterates over that sequence one argument at a time, assigning the current argument in 
the sequence to the variable called *file*.

The *$(...)* allow us to capture the output of another command (in-place of the *...*).  In this 
case we capture the output of the pipeline and save it to the variable called *dups*.

If you had multiple FASTQ files available you could run the script like so:

```sh
./fastq_dups.sh sample_1.fastq sample_2.fastq sample_3.fastq
```

And it would produce output like:

```text
sample_1.fastq 56079
sample_2.fastq XXXXX
sample_3.fastq YYYYY
```

\endshowable







## Finished

Well done, you learnt a lot over the last 5 topics and you should be proud of your achievement; it 
was a lot to take in.

From here you should be confortable around the Unix command line and ready to take on the HPC 
Workshop.

You will no-doubt forget a lot of what you learnt here so I encourage you to save a link to this 
Workshop for later reference.

Thank you for your attendance, please don't forget to complete the Melbourne Bioinformatics training survey and give it
back to the Workshop facilitators.










