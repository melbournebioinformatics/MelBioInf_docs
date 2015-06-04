<style src="../media/tute.css" ></style>
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


### Connecting to HPC

To begin this workshop you will need to connect to an HPC.  Today we will use the LIMS-HPC.  The computer called *lims-hpc-m* (m is for 
master which is another name for head node) is the one that coordinates all the HPCs tasks.

**Server details**:

* **host**: lims-hpc-m.latrobe.edu.au
* **port**: 6022 
* **username**: trainingXX (where XX is a two digit number, provided at workshop)
* **password**: (provided at workshop) 

{!docs/includes/connecting.md!}



**Note**: for security reasons ssh will not display any characters when you enter your password. This 
can be confusing because it appears as if your typing is not recognised by the computer. Don’t be 
alarmed; type your password in and press return at the end.

LIMS-HPC is a high performance computer for La Trobe Users. Logging in connects your local computer 
(e.g. laptop) to LIMS-HPC, and allows you to type commands into the Unix prompt which are run on 
the HPC, and have the results displayed on your local screen.

You will be allocated a training account on LIMS-HPC for the duration of the workshop. Your 
username and password will be supplied at the start of the workshop.

Log out of LIMS-HPC, and log back in again (to make sure you can repeat the process).

All the remaining parts assume that you are logged into LIMS-HPC over ssh.

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






## Topic 2: exploring your home directory

**Duration**: 20 minutes. See "The shell and the command line" and "The file system" section of the workshop notes.

**Relevant commands**: *ls*, *pwd*, *echo*

Your home directory contains your own private working space. Your current working directory is automatically set 
to your home directory when you log into a Unix computer.

\showable{2.1) How many files are there in your home directory?}{question}\endshowable

\showable{Answer}{answer}

The *ls* command lists files in a directory. When you first log into *LIMS-HPC* you will not have created any files of 
your own, so running the *ls* command (without any arguments) will not show any files. However you do have some files 
in your home directory, they are just not shown by default. By default, the *ls* command does not display files (or 
directories) whose names start with a dot. You can force *ls* to show all files by giving it the *-a* argument:

```sh
ls -a 
```

Now you should see several files in your home directory whose names all begin with a dot. All these files are 
created automatically for your user account. They are mostly configuration options for various programs including 
the shell. It is safe to ignore them for the moment.

\endshowable

---

\showable{2.2) What is the full path name of your home directory?}{question}\endshowable

\showable{Answer}{answer}

You can find out the full path name of the current working directory with the *pwd* command. Your home directory 
will look something like this:

```sh
/home/trainingXY
```

where *XY* is replaced by some 2 digit sequence.

You can also find out the name of your home directory by printing the value of the *$HOME* shell variable:

```sh
echo $HOME
```

\endshowable

---

\showable{2.3) Who has permission to *read*, *write* and *execute* your home directory?}{question}\endshowable

\showable{Answer}{answer}

If you pass the *-l* flag to ls it will display a "long" listing of file information including file permissions.

There are various ways you could find out the permissions on your home directory.

You could list the permissions of all files and directories in the parent directory of your home:

```sh
ls -l $HOME/..
```

In this case we use the shell variable to refer to our home directory. The *..* refers to the parent directory.

You may also refer to your home directory using the *~* (tilde) character:

```sh
ls -l ~/..
```

Both of the solutions above mean the same thing.

You will see a list of files and directories in the parent directory of your home directory. One of them will 
be the name of your home directory, something like *trainingXY*. Its permissions will look something like this 
(where *XY* is replace by a two digit string):

```sh
drwxr-x--- 4 trainingXY training  512 Feb  9 14:18 trainingXY
```

The permission string is *"drwxr-x---"*. The *d* means it is a directory. The *rwx* means that the owner of the directory 
(your user account) can *read*, *write* and *execute* the directory. Execute permissions on a directory means that you 
can *cd* into the directory. The *r-x* means that anyone in the same user group as *training* can read or execute the 
directory. The *---* means that nobody else (other users on the system) can do anything with the directory.

You can combine the *-l* and *-a* arguments to *ls* to get it to list all files in long format. If you do this in your 
home directory then you will see the permissions of all files including the current directory which is referred 
to by the special name '*.*' (a single dot character).

```sh
ls  -la
```

\endshowable

---

\showable{2.4) Can you list the contents of someone else's directory?}{question}\endshowable

\showable{Answer}{answer}

Get the username of another person in the workshop and try to list the contents of their home directory. You can 
refer to the home directory of any user by prefixing a *~* (tilde) character onto the front of their user name.
For example, suppose the username is "*foo*", we can try to list the contents of foo's home directory like so:

```sh
ls ~foo
```

or if we want a long listing of all the files in their home directory then:

```sh
ls -la ~foo
```

Which part of the permissions on a directory lets other users list the contents of the directory?

\endshowable

---

\showable{2.5) Use the man command to find out what the -h option does for ls}{question}\endshowable

\showable{Answer}{answer}

Use the following command to view the *man* page for *ls*:

```sh
man ls
```

You can scroll up and down in the man page using the arrow keys. You can search in the man page using the forward 
slash followed by the search text followed by the *enter* key. You can quit the man page by pressing *q*.

You should discover that the *-h* option does the following thing:

```sh
-h, --human-readable
              with -l, print sizes in human readable format (e.g., 1K 234M 2G)
```

\endshowable








## Topic 3: Exploring the file system

**Duration**: 30 minutes. See "The file system" section of the workshop notes.

**Relevant commands**: *pwd*, *cd*, *ls*, *file*

\showable{3.1) Print the value of your current working directory.}{question}\endshowable

\showable{Answer}{answer}

The *pwd* command prints the value of your current working directory.

\endshowable

---

\showable{3.2) Change your working directory to the root directory, called '*/*' (forward slash).}{question}\endshowable

\showable{Answer}{answer}

The *cd* command changes the value of your current working directory. To change to the root directory use the 
following command:

```sh
cd /
```

Some people imagine that changing the working directory is akin to moving your focus within the file system. 
So people often say "move to", "go to" or "charge directory to" when they want to change the working directory.

The root directory is special in Unix. It is the topmost directory in the whole file system.

\endshowable

---

\showable{3.3) List the contents of the root directory.}{question}\endshowable

\showable{Answer}{answer}

Assuming you have changed to the root directory then this can be achieved with *ls*, or *ls -a* (for all files) or 
*ls -la* for a long listing of all files.

If you are not currently in the root directory then you can list its contents by passing it as an argument to ls:

```sh
ls /
```

Here we see that *ls* can take a filepath as its argument, which allows you to list the contents of directories 
other than your current working directory.

\endshowable

---

\showable{3.4) Change your current working directory back to your home directory. What is the simplest Unix command that 
will get you back to your home directory from anywhere else in the file system?}{question}\endshowable

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

```sh
cd
```

This is a special-case behaviour which is built into *cd* for convenience.

\endshowable

---

\showable{3.5) Change your working directory to */home/group/common/training/Intro_to_Unix/*}{question}\endshowable

\showable{Answer}{answer}

```sh
cd /home/group/common/training/Intro_to_Unix/
```

\endshowable

---

\showable{3.6) List the contents of that directory. How many files does it contain?}{question}\endshowable

\showable{Answer}{answer}

You can do this with *ls*. Remember the *-a* option too.

\endshowable

---

\showable{3.7) What kind of file is */home/group/common/training/Intro_to_Unix/sleepy*?}{question}\endshowable

\showable{Answer}{answer}

Use the *file* command to get extra information about the contents of a file:

Assuming your current working directory is */home/group/common/training/Intro_to_Unix/*

```sh
file sleepy
```

Otherwise specify the full path of sleepy:

```sh
file /home/group/common/training/Intro_to_Unix/sleepy
```

You should see output like this:

```sh
Bourne-Again shell script text executable
```

The "Bource-Again shell" is more commonly known as BASH. The *file* command is telling us that sleepy 
is (probably) a shell script written in the language of BASH.

The file command uses various heuristics to guess the "type" of a file. If you want to know how it works 
then read the Unix manual page like so:

```sh
man file
```

\endshowable

---

\showable{3.8) What kind of file is */home/group/common/training/Intro_to_Unix/hi*?}{question}\endshowable

\showable{Answer}{answer}

Use the file command again. If you are in the same directory as *hi* then:

```sh
file hi
```

or from a different working directory:

```sh
file /home/group/common/training/Intro_to_Unix/hi
```

The output should look like:

```sh
ELF 64-bit LSB executable, x86-64, version 1 (SYSV), dynamically linked (uses shared libs), for GNU/Linux 
2.6.9, not stripped
```

This rather complicated output is roughly saying that the file called *hi* contains a binary executable 
program (raw instructions that the computer can execute directly).

\endshowable

---

\showable{3.9) What are the file permissions of */home/group/common/training/Intro_to_Unix/sleepy*? 
What do they mean?}{question}\endshowable

\showable{Answer}{answer}

You can find the permissions of *sleepy* using the *ls* command with the *-l* argument. If you are in the same 
directory as *sleepy* then:

```sh
ls -l sleepy
```

otherwise from a different directory:

```sh
ls -l /home/group/common/training/Intro_to_Unix/sleepy
```

The output should look like so:

```sh
-rw-r--r-- 1 training01 training 183 Feb  9 16:36 sleepy
```

We can see that this particular instance of sleepy is owned by the user training01, and is part of the training 
user group. It is 183 bytes in size, and was last modified on the 9th of February at 4:36pm. The file is 
readable to everyone, and writeable only to training01.  The digit '1' between the file permission string and 
the owner indicates that there is one link to the file. The Unix file system allows files to be referred to 
by multiple "links". When you create a file it is referred to by one link, but you may add others later. For 
future reference: links are created with the *ln* command.

\endshowable

---

\showable{3.10) Change your working directory back to your home directory.}{question}\endshowable

\showable{Answer}{answer}

You should know how to do this with the cd command:

```sh
cd
```

\endshowable






## Topic 4: Working with files and directories

**Duration**: 50 minutes. See "Working with files" from the workshop notes.

**Relevant commands**: *mkdir*, *cp*, *ls*, *diff*, *wc*, *nano*, *mv*, *head*, *tail*, *grep*, *gzip*, *gunzip*

\showable{4.1) In your home directory make a sub-directory called test.}{question}\endshowable

\showable{Answer}{answer}

Make sure you are in your home directory first. If not *cd* to your home directory.

Use the *mkdir* command to make new directories:

```sh
mkdir test
```

Use the *ls* command to check that the new directory was created.

\endshowable

---

\showable{4.2) Copy all the files from */home/group/common/training/Intro_to_Unix/* into the newly created 
test directory.}{question}\endshowable

\showable{Answer}{answer}

Use the *cp* command to copy files. You could copy them one-by-one, but that would be tedious, so use 
the *\** wildcard to specify that you want to copy all the files.

There are a number of ways you could do this depending on how you specify the source and destination 
paths to *cp*. You only need to perform one of these ways, but we show multiple ones for your reference.

From your home directory:

```sh
cp /home/group/common/training/Intro_to_Unix/* test
```

Change to the test directory and then copy (assuming you started in your home directory):

```sh
cd test
cp /home/group/common/training/Intro_to_Unix/* .
```

In the example above the '*.*' (dot) character refers to the current working directory. It should be 
the test subdirectory of your home directory.

Change to the /home/group/common/training/Intro_to_Unix/ directory and then copy:

```sh
cd /home/group/common/training/Intro_to_Unix/
cp * ~/test
```

Remember that ~ is a shortcut reference to your home directory.

\endshowable

---

\showable{4.3) Check that the file size of *expectations.txt* is the same in both the directory that you copied 
it from and the directory that you copied it to.}{question}\endshowable

\showable{Answer}{answer}

This exercise assumes that the copy command from the previous exercise was successful. 

Use *ls -l* to check the size of files.

You could do this in many ways depending on the value of your working directory. We just show one possible 
way for each file:

```sh
ls -l /home/group/common/training/Intro_to_Unix/expectations.txt

ls -l ~/test/expectations.txt
```

From the output of the above commands you should be able to see the size of each file and check that they 
are the same. They should each be *1033773* bytes.

Sometimes it is useful to get file sizes reported in more "human friendly" units than bytes. If this is 
true then try the *-h* option for ls:

```sh
ls -lh /home/group/common/training/Intro_to_Unix/expectations.txt
```

In this case the size is reported in kilobytes as *1010K*. Larger files are reported in megabytes, gigabytes 
etcetera.

\endshowable

---

\showable{4.4) Check that the contents of expectations.txt are the same in both the directory that you copied 
it from and the directory that you copied it to.}{question}\endshowable

\showable{Answer}{answer}

Use the *diff* command to compare the contents of two files. The example below assumes your working directory 
is *$HOME/test*:

```sh
diff /home/group/common/training/Intro_to_Unix/expectations.txt expectations.txt
```

The above command was too long to fit on one line in this document, so it was split over two lines using 
the backslash character "\" at the end. You don't need to type it on two lines at the unix terminal, but 
you can play with the backslash to see what effect it has if you are interested.

If the two files are identical the *diff* command will not produce any output.

\endshowable

---

\showable{4.5) How many lines, words and characters are in expectations.txt?}{question}\endshowable

\showable{Answer}{answer}

Use the *wc* (for "word count") to count the number of characters, lines and words in a file:

```sh
wc expectations.txt
```

You should see output which looks like this:

```sh
  20415  187465 1033773 expectations.txt
```

which means that there are *20415* lines, *187465* words and *1033773* characters in expectations.txt.

To get just the line count:

```sh
wc -l expectations.txt
```

To get just the word count:

```sh
wc -w expectations.txt
```

To get just the character count:

```sh
wc -c expectations.txt
```

\endshowable

---

\showable{4.6) Open *~/test/expectations.txt* in the *nano* text editor, delete the first line of text, and 
save your changes to the file. Exit *nano*.}{question}\endshowable

\showable{Answer}{answer}

Take some time to play around with the *nano* text editor.

*Nano* is a very simple text editor which is easy to use but limited in features. More powerful 
editors exist such as *vim* and *emacs*, however they take a substantial amount of time to learn.

\endshowable

---

\showable{4.7) Did the changes you made to *~/test/expectations.txt* have any effect on 
*/home/group/common/training/Intro_to_Unix/expectations.txt*? How can you tell if two files are the 
same or different in their contents?}{question}\endshowable

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

\showable{Answer}{answer}

Use the *mv* command to rename the file:

```sh
mv expectations.txt foo.txt
```

Use *ls* to check that the file is in fact renamed.

\endshowable

---

\showable{4.9) Rename foo.txt back to expectations.txt.}{question}\endshowable

\showable{Answer}{answer}

Use the *mv* command to rename the file:

```sh
mv foo.txt expectations.txt
```

Use *ls* to check that the file is in fact renamed.

\endshowable

---

\showable{4.10) Remove the file *expectations.txt* from your *test* directory.}{question}\endshowable

\showable{Answer}{answer}

Use the *rm* command to remove files (carefully):

```sh
rm expectations.txt
```

Use the *ls* command to check that the file has been removed.

\endshowable

---

\showable{4.11) Remove the entire *test* directory and all the files within it.}{question}\endshowable

\showable{Answer}{answer}

You could use the *rm* command to remove each file individually, and then use the *rmdir* command 
to remove the directory. Note that *rmdir* will only remove directories that are empty (do not 
contain files or subdirectories).

A faster way is to pass the *-r* (for recursive) flag to *rm* to remove all the files and the 
directory in one go:

```sh
cd ~
rm -r test
```

<div class="error"><b>Warning</b>: Be very careful with <em>rm -r</em>, it will remove all files and all subdirectories underneath the 
specified directory. This could be catastrophic if you do it in the wrong location! Now is a 
good moment to pause and think about file backup strategies.</div>

\endshowable

---

\showable{4.12) Recreate the test directory in your home directory and copy all the files from 
*/home/group/common/training/Intro_to_Unix/* back into the test directory.}{question}\endshowable

\showable{Answer}{answer}

Repeat exercises 4.1 and 4.2.

```sh
cd ~
mkdir test
cp /home/group/common/training/Intro_to_Unix/* test
```

\endshowable

---

\showable{4.13) Change directories to *~/test* and use the *cat* command to display the entire contents 
of the file *hello.c*}{question}\endshowable

\showable{Answer}{answer}

```sh
cd ~/test
cat hello.c
```

The output should look like this:

```c
#include <stdio.h>
int main(void) {
    printf ("Hello World\n");
    return 0;
}
```

*hello.c* contains the source code of a C program. The compiled executable version of this code 
is in the file called *hi*, which you can run like so:

```sh
./hi
```

\endshowable

---

\showable{4.14) Use the *head* command to view the first *20* lines of the file *sample_1.fastq*}{question}\endshowable

\showable{Answer}{answer}

```sh
head -20 sample_1.fastq
```

The output should look like this:

```text
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

\showable{Answer}{answer}

```sh
tail -8 sample_1.fastq
```

The output should look like this:

```text
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

\showable{Answer}{answer}

```sh
grep Ahab moby.txt
```

If you want to know how many lines are in the output of the above command you can "pipe" it 
into the *wc -l* command:

```sh
grep Ahab moby.txt | wc -l
```

which shows that there are *491* lines in *moby.txt* that contain the word Ahab.

\endshowable

---

\showable{4.17) Use the *grep* command to find out all the lines in *expectations.txt* that contain the 
word "the" with a case insensitive search (it should count "the" "The" "THE" "tHe" etcetera)
.}{question}\endshowable

\showable{Answer}{answer}

Use the *-i* flag to *grep* to make it perform case insensitive search:

```sh
grep -i the expectations.txt
```

Again, "pipe" the output to *wc -l* to count the number of lines:

```sh
grep -i the expectations.txt  | wc -l
```

\endshowable

---

\showable{4.18) Use the *gzip* command to compress the file *sample_1.fastq*. Use *gunzip* to decompress it 
back to the original contents.}{question}\endshowable

\showable{Answer}{answer}

Check the file size of sample_1.fastq before compressing it:

```sh
ls -l sample_1.fastq
```

You will see that it is *90849644* bytes in size.

Compress it with *gzip*:

```sh
gzip sample_1.fastq
```

This will replace *sample_1.fastq* with *sample_1.fastq.gz*. Check the size of the compressed file:

```sh
ls -l sample_1.fastq.gz
```

You will see that is is now *26997595* bytes in size, making it about *0.3* times the size of the 
original file.

Decompress with gunzip:

```sh
gunzip sample_1.fastq.gz
```

And check that the file size has returned to the original size (and name):

```sh
ls -l sample_1.fastq
```

\endshowable





## Topic 5: Pipes, output redirection and shell scripts

**Duration**: 50 minutes. See "Processes" from the workshop notes.

**Relevant commands**: *paste*, *grep*, *sort*, *uniq*, *wc*, *nano*, *cut*

\showable{5.1) How many reads are contained in the file *sample_1.fastq*?}{question}\endshowable

\showable{Hint}{hint}

It might be useful to work out how many lines each read takes

\showable{More|Less}

The word count(*wc*) command has a line option

\endshowable

\endshowable

\showable{Answer}{answer}

We can answer this question by counting the number of lines in the file and dividing by 4:

```sh
wc -l sample_1.fastq
```

There are *3000000* lines in the file representing *750000* reads.

If you want to do simple arithmetic at the command line then you can use the "basic calculator" 
called *bc*:

```sh
echo "3000000 / 4" | bc
```

<div class="info"><b>Note</b>: that the vertical bar character "|" is the unix pipe (and is often 
called the "pipe symbol"). It is used for connecting the output of one command into the input of 
another command. We'll see more examples soon.</div>

bc is suitable for small calculations, but it becomes cumbersome for more complex examples. If 
you want to do more sophisticated calculations then we recommend to use a more general purpose 
programming language (such as Python etcetera).

\endshowable

---

\showable{5.2) How many reads in *sample_1.fastq* contain the sequence *GATTACA*?}{question}\endshowable

\showable{Answer}{answer}

Use *grep* to find all the lines that contain *GATTACA* and "pipe" the output to *wc -l* to count them:

```sh
grep GATTACA sample_1.fastq | wc -l
```

You should see that the answer is *1119*.

If you are unsure about the possibility of upper and lower case characters then consider using 
the *-i* (ignore case option for grep).

\endshowable

---

\showable{5.3) On what line numbers do the sequences containing *GATTACA* occur?}{question}\endshowable

\showable{Answer}{answer}

You can use the *-n* flag to grep to make it prefix each line with a line number:

```sh
grep -n GATTACA sample_1.fastq
```

Or you can use the *nl* command to number each line of sample_1.fastq and then search for *GATTACA* 
in the numbered lines:

```sh
nl sample_1.fastq | grep GATTACA
```

If you just want to see the line numbers then you can "pipe" the output of the above command into 
*cut -f 1*:

```sh
 nl sample_1.fastq | grep GATTACA | cut -f 1
```

*cut* will remove certain columns from the input; in this case it will remove all except column 1
(a.k.a. field 1, hence the *-f 1* option)

\endshowable

---

\showable{5.4) Use the *nl* command to print each line of *sample_1.fastq* with its corresponding line 
number at the beginning.}{question}\endshowable

\showable{Answer}{answer}

```sh
nl sample_1.fastq
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

\showable{Answer}{answer}

```sh
nl sample_1.fastq > sample_1.fastq.nl
```

The greater-than sign "*>*" is the file redirection operator. It causes the standard output of the 
command on the left-hand-side to be written to the file on the right-hand-side.

You should notice that the above command is much faster than printing the output to the screen. 
This is because writing to disk can be performed much more quickly than rendering the output on 
a terminal.

To check that the first 20 lines of the file look reasonable you can use the *head* command like so:

```sh
head -20 sample_1.fastq.nl
```

The *less* command allows you to interactively view a file. The arrow keys move the page up and 
down. You can search using the '*/*' followed by the search term. You can quit by pressing "*q*". Note 
that the *less* command is used by default to display man pages.

```sh
less sample_1.fastq.nl
```

\endshowable

---

\showable+{5.6) The four-lines-per-read format of FASTQ is cumbersome to deal with. Often it would be 
preferable if we could convert it to tab-separated-value (TSV) format, such that each read appears 
on a single line with each of its fields separated by tabs. Use the following command to convert 
sample_1.fastq.nl into TSV format:}{question}

```sh
cat sample_1.fastq | paste - - - - > sample_1.tsv
```

\endshowable

\showable{Answer}{answer}

The *'-'* (dash) character has a special meaning when used in place of a file; it means use the standard
input instead of a real file.  Note: while it is fairly common in most unix programs, not all wil support it.

The *paste* command is useful for merging multiple files together line-by-line, such that the *Nth* 
line from each file is joined together into one line in the output, separated by default with a 
*tab* character. In the above example we give paste 4 copies of the contents of *sample_1.fastq*, 
which causes it to join consecutive groups of 4 lines from the file into one line of output.

\endshowable

---

\showable+{5.7) Do you expect the output of the following command to produce the same output as above? and why?}{question}

```sh
paste sample_1.fastq sample_1.fastq sample_1.fastq sample_1.fastq > sample_1b.tsv
```

Try it, see what ends up in sample_1b.tsv (maybe use *less*)

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

\showable{Answer}{answer}

We can count the number of lines in *sample_1.tsv* using *wc*:

```sh
wc -l sample_1.tsv
```

The output should be *750000* as expected (1/4 of the number of lines in sample_1.fastq).

To view the first *20* lines of *sample_1.tsv* use the *head* command:

```sh
head -20 sample_1.tsv
```

\endshowable

---

\showable{5.9) Use the *cut* command to print out the second column of *sample_1.tsv*. Redirect the 
output to a file called *sample_1.dna.txt*.}{question}\endshowable

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

\showable{Answer}{answer}

```sh
sort sample_1.dna.txt > sample_1.dna.sorted.txt
```

Running *head* on the output file reveals that there are duplicate DNA sequences in the input FASTQ 
file.

\endshowable

---

\showable{5.11) Use the *uniq* command to remove duplicate consecutive lines from *sample_1.dna.sorted.txt*, 
redirect the result to *sample_1.dna.uniq.txt*. Compare the number of lines in sample1_dna.txt to 
the number of lines in *sample_1.dna.uniq.txt*.}{question}\endshowable

\showable{Answer}{answer}

```sh
uniq sample_1.dna.sorted.txt > sample_1.dna.uniq.txt
```

Compare the outputs of:

```sh
wc -l sample_1.dna.sorted.txt
```

and

```sh
wc -l sample_1.dna.uniq.txt
```

View the contents of *sample_1.dna.uniq.txt* to check that the duplicate DNA sequences have been 
removed.

\endshowable

---

\showable{5.12) Can you modify the command from above to produce *only* those sequences of DNA which were 
duplicated in *sample_1.dna.sorted.txt*?}{question}\endshowable

\showable{Hint}{hint}

Look at the man page for uniq.

\endshowable

\showable{Answer}{answer}

Use the *-d* flag to *uniq* to print out only the duplicated lines from the file:

```sh
uniq -d sample_1.dna.sorted.txt > sample_1.dna.dup.txt
```

\endshowable

---

\showable{5.13) Write a *shell pipeline* which will print the number of duplicated DNA sequences in 
sample_1.fastq.}{question}\endshowable

\showable{Answer}{answer}

Finally we can 'pipe' all the pieces together into a sophisticated pipeline which starts with a 
FASTQ file and ends with a list of duplicated DNA sequences:

```sh
cat sample_1.fastq | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
```

The output file should have *56079* lines.

\endshowable

---

\showable{5.14) (Advanced) Write a shell script which will print the number of duplicated DNA sequences 
in sample_1.fastq.}{question}\endshowable

\showable{Hint}{hint}

Check out the *sleepy* file (with *cat* or *nano*); there is a bit of magic on the first line that you will need. 

You also need to tell bash that this file can be executed (check out *chmod*)

\endshowable

\showable{Answer}{answer}

Put the answer to *5.13* into a file called *sample_1_dups.sh* (or whatever you want). Use *nano* to 
create the file. The contents of the file will look like this:

```sh
#!/bin/bash

cat sample_1.fastq | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
```

<div class="info"><b>Note</b>: the first line has special meaning.  If it starts with '<em>#!</em>' (Hash 
then exclamation mark) then it tells bash this file is a script that can be interpreted.  The command 
(including full path) used to intepret the script is placed right after the magic code.</div>

Give everyone execute permissions on the file with chmod:

```sh
chmod +x sample_1_dups.sh 
```

You can run the script like so:

```sh
./sample_1_dups.sh
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
cp sample_1_dups.sh fastq_dups.sh
```

Edit the new shell script file and change it to use the command line parameters:

```sh
#!/bin/bash

cat $1 | paste - - - - | cut -f 2 | sort | uniq -d | wc -l
```

You can run the new script like so:

```sh
./fastq_dups.sh sample_1.fastq
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
It is a unix standard that when the user provides incorrect commandline arguments we print a usage message 
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

Thank you for your attendance, please don't forget to complete the VLSCI training survey and give it
back to the Workshop facilitators.










