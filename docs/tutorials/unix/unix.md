<style src="../../includes/media/tute.css" ></style>
<style>em {font-style: normal; font-family: courier new;}</style>

![melbioinf_logo](/img/melbioinf_logo.png){: style="width:350px; padding-right:50px"}       ![unimelb_logo](/img/PRIMARY_A_Vertical_Housed_RGB.png){: style="width:150px"}

# Introduction to Unix

A hands-on workshop covering the basics of the Unix/Linux command line interface.

## Overview

Knowledge of the Unix operating system is fundamental to being
productive on HPC systems. This workshop will introduce you to the
fundamental Unix concepts by way of a series of hands-on exercises.

The workshop is facilitated by experienced Unix users who will be able
to guide you through the exercises and offer assistance where needed.

## Learning Objectives

At the end of the course, you will be able to:

* Log into a Unix machine remotely
* Organise your files into directories
* Change file permissions to improve security and safety
* Create and edit files with a text editor
* Copy files between directories
* Use command line programs to manipulate files
* Automate your workflow using shell scripts

## Requirements

* The workshop is intended for beginners with no prior experience in Unix.
* Attendees are required to bring their own laptop computers.


## Introduction

Before we commence the hands-on part of this workshop we will first give a short 30 minute talk to introduce the Unix concepts.
The [slides](slides.html) are available if you would like.  Additionally the following reference material is available for later
use.

<details>
  <summary>Reference Material</summary>

{!docs/tutorials/unix/intro.md!}

</details>

## Topic 1: Remote log in

In this topic we will learn how to connect to a *Unix* computer via a program called *ssh* and run a few basic commands.

### Connecting to a Unix computer

To begin this workshop you will need to connect to an HPC.  Today we will use *barcoo*.
The computer called *barcoo.vlsci.org.au* is the one that coordinates all the HPC's tasks.

**Server details**:

* **host**: barcoo.vlsci.org.au
* **port**: 22
* **username**: (provided at workshop)
* **password**: (provided at workshop)

{!docs/includes/connecting.md!}


**Note**: for security reasons ssh will not display any characters when you enter your password. This
can be confusing because it appears as if your typing is not recognised by the computer. Don’t be
alarmed; type your password in and press return at the end.

barcoo is a high performance computer for Melbourne Bioinformatics users.  Logging in connects your local computer
(e.g. laptop) to barcoo, and allows you to type commands into the Unix prompt which are run on
the HPC, and have the results displayed on your local screen.

You will be allocated a training account on barcoo for the duration of the workshop. Your
username and password will be supplied at the start of the workshop.

Log out of barcoo, and log back in again (to make sure you can repeat the process).

All the remaining parts assume that you are logged into barcoo over ssh.

### Exercises

#### 1.1) When you've logged into the Unix server, run the following commands and see what they do:

* *who*
* *whoami*
* *date*
* *cal*
* *hostname*
* */vlsci/TRAINING/shared/Intro_to_Unix/hi*

<details>
  <summary>Answer</summary>

* **who**: displays a list of the users who are currently using this Unix computer.
* **whoami**: displays your username (i.e. they person currently logged in).
* **date**: displays the current date and time.
* **cal**: displays a calendar on the terminal.  It can be configured to display more than just
the current month.
* **hostname**: displays the name of the computer we are logged in to.
* **/vlsci/TRAINING/shared/Intro_to_Unix/hi**: displays the text "Hello World"

</details>


## Topic 2: Exploring your home directory

In this topic we will learn how to "look" at the filesystem and further expand our repertoire of Unix commands.

**Duration**: 20 minutes. <!-- See "The shell and the command line" and "The file system" section of the workshop notes. -->

**Relevant commands**: *ls*, *pwd*, *echo*, *man*

Your home directory contains your own private working space.  Your *current working directory* is automatically set
to your *home* directory when you log into a Unix computer.

#### 2.1) Use the *ls* command to list the files in your *home* directory.  How many files are there?

<details>
  <summary>Hint</summary>

Literally, type *ls* and press the *ENTER* key.

</details>

<details>
  <summary>Answer</summary>

```sh
$ ls
exp01  file01  muscle.fq
```

When running the *ls* command with no options it will list files in your current working directory.  The place
where you start when you first login is your *HOME* directory.

**Answer**: 3 (exp01, file01 and muscle.fq)

</details>

---

The above answer is not quite correct.  There are a number of *hidden* files in your home directory as well.

#### 2.2) What *flag* might you use to display *all* files with the *ls* command?  How many files are really there?

<details>
  <summary>Hint</summary>

Take the *all* quite literally.

</details>

<details>
  <summary>Additional Hint</summary>

Type *ls --all* and press the *ENTER* key.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 2.3) What is the full path name of your *home* directory?

<details>
  <summary>Hint</summary>

Remember your *Current Working Directory* starts in your *home* directory.

</details>

<details>
  <summary>Additional Hint</summary>

Try a shortened version of *print working directory*

</details>

<details>
  <summary>Answer</summary>

You can find out the full path name of the current working directory with the *pwd* command. Your home directory
will look something like this:

```sh
$ pwd
/home/trainingXX
```

**Answer**: */vlsci/TRAINING/trainXX*

where *XX* is replaced by some 2 digit sequence.

**Alternate method**:
You can also find out the name of your home directory by printing the value of the *$HOME* shell variable:

```sh
echo $HOME
```

</details>

---

#### 2.4) Run *ls* using the long flag (*-l*), how did the output change?

<details>
  <summary>Hint</summary>

Run *ls -l*

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 2.5) What type of file is *exp01* and *muscle.fq*?

<details>
  <summary>Hint</summary>

Check the output from the *ls -l*.

</details>

<details>
  <summary>Answer</summary>

**Answer**:

* *exp01*: Directory (given the 'd' as the first letter of its permissions)
* *muscle.fq*: Regular File (given the '-')

</details>

---


#### 2.6) Who has permission to *read*, *write* and *execute* your *home* directory?

<details>
  <summary>Hint</summary>

You can also give *ls* a filename as the first option.

</details>

<details>
  <summary>Additional Hint</summary>

*ls -l* will show you the contents of the *CWD*; how might you see the contents of the *parent* directory? (remember
the slides)

</details>

<details>
  <summary>Answer</summary>

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
be the name of your home directory, something like *trainXX*.  Where *XX* is replaced by a two digit string.

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

</details>

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

#### 2.7) Use the *man* command to find out what the *-h* flag does for *ls*

<details>
  <summary>Hint</summary>

Give *ls* as an option to *man* command.

</details>

<details>
  <summary>Additional Hint</summary>

*man ls*

</details>

<details>
  <summary>Answer</summary>

Use the following command to view the *man* page for *ls*:

```sh
$ man ls
```

**Answer**: You should discover that the *-h* option prints file sizes in human readable format

```sh
-h, --human-readable
              with -l, print sizes in human readable format (e.g., 1K 234M 2G)
```

</details>


---

#### 2.8) Use the *-h*, how did the output change of *muscle.fq*?

<details>
  <summary>Hint</summary>

Don't forget the *-l* option too.

</details>

<details>
  <summary>Additional Hint</summary>

Run *ls -lh*

</details>

<details>
  <summary>Answer</summary>

```sh
$ ls -lh
...
-rw-r----- 1 training01 training 2.5K Jun 14 11:28 muscle.fq
```

**Answer**: it changed the output so the *filesize* of *muscle.fq* is now *2.5K* instead of *2461*

</details>









## Topic 3: Exploring the file system

In this topic we will learn how to move around the filesystem and see what is there.

**Duration**: 30 minutes. <!-- See "The file system" section of the workshop notes. -->

**Relevant commands**: *pwd*, *cd*, *ls*, *file*

#### 3.1) Print the value of your current working directory.

<details>
  <summary>Answer</summary>

The *pwd* command prints the value of your current working directory.

```sh
$ pwd
/home/training01
```

</details>

---

#### 3.2) List the contents of the root directory, called '*/*' (forward slash).

<details>
  <summary>Hint</summary>

*ls* expects one or more anonymous options which are the files/directories to list.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 3.3) Use the *cd* command to change your working directory to the root directory.  Did your prompt change?

<details>
  <summary>Hint</summary>

*cd* expects a single option which is the directory to change to

</details>

<details>
  <summary>Answer</summary>

The *cd* command changes the value of your current working directory. To change to the root directory use the
following command:

```sh
$ cd /
```

**Answer**: Yes, it now says the CWD is */* instead of *~*.

Some people imagine that changing the working directory is akin to moving your focus within the file system.
So people often say "move to", "go to" or "charge directory to" when they want to change the working directory.

The root directory is special in Unix. It is the topmost directory in the whole file system.

</details>

---

<div class="info">
<b>Output on ERROR only</b>: Many Unix commands will not produce any output if everything went well; <em>cd</em> is one
such command.  However, it will get grumpy if something went wrong by way of an error message on-screen.
</div>

---

#### 3.4) List the contents of the CWD and verify it matches the list in 3.2

<details>
  <summary>Hint</summary>

*ls*

</details>

<details>
  <summary>Answer</summary>

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



</details>

---

#### 3.5) Change your current working directory back to your home directory. What is the simplest Unix command that will get you back to your home directory from anywhere else in the file system?

<details>
  <summary>Hint</summary>

The answer to exercise 2.6 might give some hints on how to get back to the home directory

</details>

<details>
  <summary>Additional Hint</summary>

*$HOME*, *~*, */vlsci/TRAINING/trainXX* are all methods to name your home directory.  Yet there is a simpler method; the answer
is buried in *man cd* however *cd* doesn't have its own manpage so you will need to search for it.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 3.6) Change your working directory to the following directory:

*/vlsci/TRAINING/shared/Intro_to_Unix*

<details>
  <summary>Answer</summary>

**Answer**: *cd /vlsci/TRAINING/shared/Intro_to_Unix*

</details>

---

#### 3.7) List the contents of that directory. How many files does it contain?

<details>
  <summary>Hint</summary>

*ls*

</details>

<details>
  <summary>Answer</summary>

You can do this with *ls*

```sh
$ ls
expectations.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

**Answer**: 7 files (expectations.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy)

</details>

---

#### 3.8) What kind of *file* is */vlsci/TRAINING/shared/Intro_to_Unix/sleepy*?

<details>
  <summary>Hint</summary>

Take the word *file* quite literally.

</details>

<details>
  <summary>Additional Hint</summary>

*file sleepy*

</details>

<details>
  <summary>Answer</summary>

Use the *file* command to get extra information about the contents of a file:

Assuming your current working directory is */vlsci/TRAINING/shared/Intro_to_Unix*

```sh
$ file sleepy
Bourne-Again shell script text executable
```

Otherwise specify the full path of sleepy:

```sh
$ file /vlsci/TRAINING/shared/Intro_to_Unix/sleepy
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

</details>

---

3.9) What kind of *file* is */vlsci/TRAINING/shared/Intro_to_Unix/hi*?

<details>
  <summary>Hint</summary>

Take the word *file* quite literally.

</details>

<details>
  <summary>Answer</summary>

Use the file command again. If you are in the same directory as *hi* then:

```sh
$ file hi
ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked (uses shared libs), for GNU/Linux
2.6.9, not stripped
```

**Answer**: ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked (uses shared libs), for GNU/Linux

This rather complicated output is roughly saying that the file called *hi* contains a binary executable
program (raw instructions that the computer can execute directly).

</details>

---

#### 3.10) What are the file permissions of the following file and what do they mean?

*/vlsci/TRAINING/shared/Intro_to_Unix/sleepy*

<details>
  <summary>Hint</summary>

Remember the *ls* command, and don't forget the *-l* flag

</details>

<details>
  <summary>Answer</summary>

You can find the permissions of *sleepy* using the *ls* command with the *-l* flag. If you are in the same
directory as *sleepy* then:

```sh
$ ls -l sleepy
-rw-r--r-- 1 arobinson common 183 Feb  9 16:36 sleepy
```

**Answer**: The Answer is dependent on the computer you are connected too however will follow something like above.
We can see that this particular instance of sleepy is owned by the user arobinson, and is part of the common
user group. It is 183 bytes in size, and was last modified on the 9th of February at 4:36pm. The file is
readable to everyone, and write-able only to arobinson.  The digit '1' between the file permission string and
the owner indicates that there is one link to the file. The Unix file system allows files to be referred to
by multiple "links". When you create a file it is referred to by one link, but you may add others later. For
future reference: links are created with the *ln* command.

</details>

---

#### 3.11) Change your working directory back to your home directory ready for the next topic.

<details>
  <summary>Hint</summary>

*cd*

</details>

<details>
  <summary>Answer</summary>

You should know how to do this with the cd command:

```sh
cd
```

</details>









## Topic 4: Working with files and directories

In this topic we will start to read, create, edit and delete files and directories.

**Duration**: 50 minutes.  <!-- See "Working with files" from the workshop notes. -->

**Relevant commands**: *mkdir*, *cp*, *ls*, *diff*, *wc*, *nano*, *mv*, *rm*, *rmdir*, *head*, *tail*, *grep*, *gzip*, *gunzip*

<div class="info">
<b>Hint</b>: Look at the commands above; you will need them roughly in order for this topic.  Use the <em>man</em>
command find out what they do, in particular the NAME, SYNOPSIS and DESCRIPTION sections.
</div>

---


#### 4.1) In your home directory make a sub-directory called test.

<details>
  <summary>Hint</summary>

You are trying to *make a directory*, which of the above commands looks like a shortened version of this?

</details>

<details>
  <summary>Additional Hint</summary>

*mkdir*

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.2) Copy all the files from the following directory into the newly created test directory:

*/vlsci/TRAINING/shared/Intro_to_Unix*

<details>
  <summary>Hint</summary>

You are trying to *copy*, which of the above commands looks like a shortened version of this?

</details>

<details>
  <summary>Additional Hint</summary>

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

</details>

<details>
  <summary>Answer</summary>

Use the *cp* command to copy files.

<div class="info"><b>Wildcards</b>: You could copy them one-by-one, but that would be tedious, so use
the <em>*</em> wildcard to specify that you want to copy all the files.
</div>

There are a number of ways you could do this depending on how you specify the source and destination
paths to *cp*. You only need to perform one of these ways, but we show multiple ones for your reference.

**Answer 1**: From your home directory:

```sh
$ cp /vlsci/TRAINING/shared/Intro_to_Unix/* test
```

**Answer 2**: Change to the test directory and then copy (assuming you started in your home directory):

```sh
$ cd test
$ cp /vlsci/TRAINING/shared/Intro_to_Unix/* .
```

In the example above the '*.*' (dot) character refers to the current working directory. It should be
the test subdirectory of your home directory.

**Answer 3**: Change to the \end{UNIX_TRAINING_FILES_PATH} directory and then copy:

```sh
cd /vlsci/TRAINING/shared/Intro_to_Unix/
cp * ~/test
```

Remember that ~ is a shortcut reference to your home directory.

</details>

---

**Note**: This exercise assumes that the copy command from the previous exercise was successful.

#### 4.3) Check that the file size of *expectations.txt* is the same in both the directory that you copied it from and the directory that you copied it to.

<details>
  <summary>Hint</summary>

Remember *ls* can show you the file size (with one of its flags)

</details>

<details>
  <summary>Additional Hint</summary>

*ls -l*

</details>

<details>
  <summary>Answer</summary>

Use *ls -l* to check the size of files.

You could do this in many ways depending on the value of your working directory. We just show one possible
way for each file:

```sh
$ ls -l /vlsci/TRAINING/shared/Intro_to_Unix/expectations.txt

$ ls -l ~/test/expectations.txt
```

From the output of the above commands you should be able to see the size of each file and check that they
are the same.

**Answer**: They should each be *1033773* bytes

**Alternate**: Sometimes it is useful to get file sizes reported in more "human friendly" units than bytes. If this is
true then try the *-h* option for ls:

```sh
$ ls -lh /vlsci/TRAINING/shared/Intro_to_Unix/expectations.txt
-rw-r--r-- 1 arobinson common 1010K Mar 26  2012 /vlsci/TRAINING/shared/Intro_to_Unix/expectations.txt
```

In this case the size is reported in kilobytes as *1010K*. Larger files are reported in megabytes, gigabytes
etcetera.

</details>

---

**Note**: this exercise assumes your working directory is *~/test*; if not run *cd ~/test*

#### 4.4) Check that the contents of expectations.txt are the same in both the directory that you copied it from and the directory that you copied it to.

<details>
  <summary>Hint</summary>

What is the opposite of *same*?

</details>

<details>
  <summary>Additional Hint</summary>

*diff*erence

</details>

<details>
  <summary>Answer</summary>

Use the *diff* command to compare the contents of two files.

```sh
$ diff /vlsci/TRAINING/shared/Intro_to_Unix/expectations.txt expectations.txt
```

If the two files are identical the *diff* command will NOT produce any output)

**Answer**: Yes, they are the same since no output was given.

</details>

---

#### 4.5) How many lines, words and characters are in expectations.txt?

<details>
  <summary>Hint</summary>

Initialisms are key

</details>

<details>
  <summary>Additional Hint</summary>

*w*ord *c*ount

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.6) Open *~/test/expectations.txt* in the *nano* text editor, delete the first line of text, and save your changes to the file. Exit *nano*.

<details>
  <summary>Hint</summary>

*nano FILENAME*

Once *nano* is open it displays some command hints along the bottom of the screen.

</details>

<details>
  <summary>Additional Hint</summary>

*^O* means hold the *Control* (or CTRL) key while pressing the *o*.  Despite what it displays, you need to type
the lower-case letter that follows the *^* character.

WriteOut is another name for Save.

</details>

<details>
  <summary>Answer</summary>

Take some time to play around with the *nano* text editor.

*Nano* is a very simple text editor which is easy to use but limited in features. More powerful
editors exist such as *vim* and *emacs*, however they take a substantial amount of time to learn.

</details>

---

#### 4.7) Did the changes you made to *~/test/expectations.txt* have any effect on */vlsci/TRAINING/shared/Intro_to_Unix*?

How can you tell if two files are the same or different in their contents?

<details>
  <summary>Hint</summary>

Remember exercise 4.4

</details>

<details>
  <summary>Additional Hint</summary>

Use *diff*

</details>

<details>
  <summary>Answer</summary>

Use *diff* to check that the two files are different after you have made the change to the copy of
*expectations.txt* in your *~/test* directory.

```sh
diff ~/test/expectations.txt \
/vlsci/TRAINING/shared/Intro_to_Unix/expectations.txt
```

You could also use *ls* to check that the files have different sizes.

</details>

---

#### 4.8) In your *test* subdirectory, rename *expectations.txt* to *foo.txt*.

<details>
  <summary>Hint</summary>

Another way to think of it is *moving* it from *expectations.txt* to *foo.txt*

</details>

<details>
  <summary>Additional Hint</summary>

*mv*

Use *man mv* if you need to work out how to use it.

</details>

<details>
  <summary>Answer</summary>

Use the *mv* command to rename the file:

```sh
$ mv expectations.txt foo.txt
$ ls
foo.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

</details>

---

#### 4.9) Rename foo.txt back to expectations.txt.

<details>
  <summary>Answer</summary>

Use the *mv* command to rename the file:

```sh
$ mv foo.txt expectations.txt
$ ls
expectations.txt  hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

Use *ls* to check that the file is in fact renamed.

</details>

---

#### 4.10) Remove the file *expectations.txt* from your *test* directory.

<details>
  <summary>Hint</summary>

We are trying to *remove* a file, check the commands at the top of this topic.

</details>

<details>
  <summary>Additional Hint</summary>

*rm*

</details>

<details>
  <summary>Answer</summary>

Use the *rm* command to remove files (carefully):

```sh
$ rm expectations.txt
$ ls
hello.c  hi  jude.txt  moby.txt  sample_1.fastq  sleepy
```

</details>

---

#### 4.11) Remove the entire *test* directory and all the files within it.

<details>
  <summary>Hint</summary>

We are trying to *remove a directory*.

</details>

<details>
  <summary>Additional Hint</summary>

You could use *rmdir* but there is an easier way using just *rm* and a flag.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.12) Recreate the test directory in your home directory and copy all the files from */vlsci/TRAINING/shared/Intro_to_Unix* back into the test directory.

<details>
  <summary>Hint</summary>

See exercises 4.1 and 4.2

</details>

<details>
  <summary>Answer</summary>

Repeat exercises 4.1 and 4.2.

```sh
$ cd ~
$ mkdir test
$ cp /vlsci/TRAINING/shared/Intro_to_Unix/* test
```

</details>

---

#### 4.13) Change directories to *~/test* and use the *cat* command to display the entire contents of the file *hello.c*

<details>
  <summary>Hint</summary>

Use *man* if you can't guess how it might work.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.14) Use the *head* command to view the first *20* lines of the file *sample_1.fastq*

<details>
  <summary>Hint</summary>

Remember your *best* friend!

</details>

<details>
  <summary>Additional Hint</summary>

Use *man* to find out what option you need to add to display a given number of *lines*.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.15) Use the *tail* command to view the last *8* lines of the file *sample_1.fastq*

<details>
  <summary>Hint</summary>

It's very much like *head*.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.16) Use the *grep* command to find out all the lines in *moby.txt* that contain the word "Ahab"

<details>
  <summary>Hint</summary>

One might say we are 'looking for the *pattern* "Ahab"'

</details>

<details>
  <summary>Additional Hint</summary>

```sh
$ man grep
...
SYNOPSIS
       grep [OPTIONS] PATTERN [FILE...]
...
```

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.17) Use the *grep* command to find out all the lines in *expectations.txt* that contain the word "the" with a case insensitive search (it should count "the" "The" "THE" "tHe" etcetera).

<details>
  <summary>Hint</summary>

One might say we are *ignoring case*.

</details>

<details>
  <summary>Additional Hint</summary>

```sh
$ man grep
...
       -i, --ignore-case
              Ignore case distinctions in both the PATTERN and the input files.  (-i is specified by POSIX.)
...
```

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 4.18) Use the *gzip* command to compress the file *sample_1.fastq*. Use *gunzip* to decompress it back to the original contents.

<details>
  <summary>Hint</summary>

Use the above commands along with *man* and *ls* to see what happens to the file.

</details>

<details>
  <summary>Answer</summary>

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

</details>





## Topic 5: Pipes, output redirection and shell scripts

In this section we will cover a lot of the more advanced Unix concepts; it is here where you will start to see
the power of Unix.  I say *start* because this is only the "tip of the iceberg".

**Duration**: 50 minutes. <!-- See "Processes" from the workshop notes. -->

**Relevant commands**: *wc*, *paste*, *grep*, *sort*, *uniq*, *nano*, *cut*




#### 5.1) How many *reads* are contained in the file *sample_1.fastq*?

<details>
  <summary>Hint</summary>

Examine some of the file to work out how many lines each *read* takes up.

</details>

<details>
  <summary>Additional Hint</summary>

Count the number of lines

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 5.2) How many reads in *sample_1.fastq* contain the sequence *GATTACA*?

<details>
  <summary>Hint</summary>

Check out exercise 4.16

</details>

<details>
  <summary>Answer</summary>

Use *grep* to find all the lines that contain *GATTACA* and "pipe" the output to *wc -l* to count them:

```sh
$ grep GATTACA sample_1.fastq | wc -l
1119
```

**Answer**: *1119*

If you are unsure about the possibility of upper and lower case characters then consider using
the *-i* (ignore case option for grep).

</details>

---

#### 5.3) On what line numbers do the sequences containing *GATTACA* occur?

<details>
  <summary>Hint</summary>

We are looking for the *line numbers*.

</details>

<details>
  <summary>Additional Hint</summary>

Check out the manpage for *grep* and/or *nl*

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 5.4) Use the *nl* command to print each line of *sample_1.fastq* with its corresponding line number at the beginning.

<details>
  <summary>Hint</summary>

Check answer to 5.3.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 5.5) Redirect the output of the previous command to a file called *sample_1.fastq.nl*.

Check the first *20* lines of *sample_1.fastq.nl* with the *head* command. Use the *less* command to
interactively view the contents of *sample_1.fastq.nl* (use the arrow keys to navigate up and down,
*q* to quit and '*/*' to search). Use the search facility in less to find occurrences of
*GATTACA*.

<details>
  <summary>Hint</summary>

Ok that one was tough, *> FILENAME* is how you do it if you didn't break out an internet search for
"redirect the output in Unix"

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 5.6) The four-lines-per-read format of FASTQ is cumbersome to deal with. Often it would be preferable if we could convert it to tab-separated-value (TSV) format, such that each read appears on a single line with each of its fields separated by tabs. Use the following command to convert sample_1.fastq into TSV format:

```sh
$ cat sample_1.fastq | paste - - - - > sample_1.tsv
```

<details>
  <summary>Answer</summary>

The *'-'* (dash) character has a special meaning when used in place of a file; it means use the standard
input instead of a real file.  Note: while it is fairly common in most Unix programs, not all will support it.

The *paste* command is useful for merging multiple files together line-by-line, such that the *Nth*
line from each file is joined together into one line in the output, separated by default with a
*tab* character. In the above example we give paste 4 copies of the contents of *sample_1.fastq*,
which causes it to join consecutive groups of 4 lines from the file into one line of output.

</details>

---

#### 5.7) Do you expect the output of the following command to produce the same output as above? and why?

```sh
$ paste sample_1.fastq sample_1.fastq sample_1.fastq sample_1.fastq > sample_1b.tsv
```

Try it, see what ends up in sample_1b.tsv (maybe use *less*)

<details>
  <summary>Hint</summary>

Use *less* to examine it.

</details>

<details>
  <summary>Answer</summary>

**Answer**: No, in the second instance we get 4 copies of each line.

**Why**: In the first command *paste* will use the input file (standard input) 4 times since the *cat*
command will only give one copy of the file to *paste*, where as, in the second command *paste* will open
the file 4 times.  Note: this is quite confusing and is not necessory to remember; its just an interesting
side point.

</details>

---

#### 5.8) Check that *sample_1.tsv* has the correct number of lines. Use the *head* command to view the first *20* lines of the file.

<details>
  <summary>Hint</summary>

Remember the *wc* command.

</details>

<details>
  <summary>Answer</summary>

We can count the number of lines in *sample_1.tsv* using *wc*:

```sh
$ wc -l sample_1.tsv
```

The output should be *750000* as expected (1/4 of the number of lines in sample_1.fastq).

To view the first *20* lines of *sample_1.tsv* use the *head* command:

```sh
$ head -20 sample_1.tsv
```

</details>

---

#### 5.9) Use the *cut* command to print out the second column of *sample_1.tsv*. Redirect the output to a file called *sample_1.dna.txt*.

<details>
  <summary>Hint</summary>

See exercise 5.3 (for cut) and 5.5 (redirection)

</details>

<details>
  <summary>Answer</summary>

The file sample_1.tsv is in column format. The cut command can be used to select certain columns
from the file. The DNA sequences appear in column 2, we select that column using the -f 2 flag
(the f stands for "field").

```sh
cut -f 2 sample_1.tsv > sample_1.dna.txt
```

Check that the output file looks reasonable using *head* or *less*.

</details>

---

#### 5.10) Use the *sort* command to sort the lines of *sample_1.dna.txt* and redirect the output to *sample_1.dna.sorted.txt*. Use *head* to look at the first few lines of the output file. You should see a lot of repeated sequences of As.

<details>
  <summary>Hint</summary>

Use *man* (sort) and see exercise 5.5 (redirection)

</details>

<details>
  <summary>Answer</summary>

```sh
$ sort sample_1.dna.txt > sample_1.dna.sorted.txt
```

Running *head* on the output file reveals that there are duplicate DNA sequences in the input FASTQ
file.

</details>

---

#### 5.11) Use the *uniq* command to remove duplicate consecutive lines from *sample_1.dna.sorted.txt*, redirect the result to *sample_1.dna.uniq.txt*. Compare the number of lines in sample1_dna.txt to the number of lines in *sample_1.dna.uniq.txt*.

<details>
  <summary>Hint</summary>

I am pretty sure you have already used *man* (or just guessed how to use *uniq*).  You're also a gun at
redirection now.

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 5.12) Can you modify the command from above to produce *only* those sequences of DNA which were duplicated in *sample_1.dna.sorted.txt*?

<details>
  <summary>Hint</summary>

Checkout the *uniq* manpage

</details>

<details>
  <summary>Additional Hint</summary>

Look at the man page for uniq.

</details>

<details>
  <summary>Answer</summary>

Use the *-d* flag to *uniq* to print out only the duplicated lines from the file:

```sh
$ uniq -d sample_1.dna.sorted.txt > sample_1.dna.dup.txt
```

</details>

---

#### 5.13) Write a *shell pipeline* which will print the number of duplicated DNA sequences in sample_1.fastq.

<details>
  <summary>Hint</summary>

That is, *piping* most of the commands you used above instead of redirecting to file

</details>

<details>
  <summary>Additional Hint</summary>

i.e. 6 commands (*cat*, *paste*, *cut*, *sort*, *uniq*, *wc*)

</details>

<details>
  <summary>Answer</summary>

Finally we can 'pipe' all the pieces together into a sophisticated pipeline which starts with a
FASTQ file and ends with a list of duplicated DNA sequences:

**Answer**:
```sh
$ cat sample_1.fastq | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
56079
```

The output file should have *56079* lines.

</details>

---

#### 5.14) (Advanced) Write a shell script which will print the number of duplicated DNA sequences in sample_1.fastq.

<details>
  <summary>Hint</summary>

Check out the *sleepy* file (with *cat* or *nano*); there is a bit of magic on the first line that you will need.

You also need to tell bash that this file can be executed (check out *chmod* command).

</details>

<details>
  <summary>Answer</summary>

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

</details>

---

#### 5.15) (Advanced) Modify your shell script so that it accepts the name of the input FASTQ file as a command line parameter.

<details>
  <summary>Hint</summary>

Shell scripts can refer to command line arguments by their position using special variables called
*$0*, *$1*, *$2* and so on.

</details>

<details>
  <summary>Additional Hint</summary>

*$0* refers to the name of the script as it was called on the command line.
*$1* refers to the first command line argument, and so on.

</details>


<details>
  <summary>Answer</summary>

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
* *-eq*: checks if the numbers on either side of it are equal.
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

</details>

---

#### 5.16) (Advanced) Modify your shell script so that it accepts zero or more FASTQ files on the command line argument and outputs the number of duplicated DNA sequences in each file.

<details>
  <summary>Answer</summary>

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

</details>

## Finished

Well done, you learnt a lot over the last 5 topics and you should be proud of your achievement; it
was a lot to take in.

From here you should be comfortable around the Unix command line and ready to take on the HPC
Workshop.

You will no-doubt forget a lot of what you learnt here so I encourage you to save a link to this
workshop for later reference.

Thank you for your attendance, please don't forget to complete the training survey and give it
back to the workshop facilitators.
