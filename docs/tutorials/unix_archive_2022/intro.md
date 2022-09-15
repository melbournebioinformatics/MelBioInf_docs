
### Why use unix

* **Powerful**: Unix computers are typically very powerful in comparision to your desktop/laptop computers.
  Additionally they don't typically use a Graphical User Interface which can free up much resources for actual
  computing.
* **Big data**: Unix programs are designed to handle large data sets
* **Flexible**: Small programs that can be arranged in many ways to solve your problems
* **Automation**: Scripting allows you to do many tasks in one step and repeat steps many times
* **Pipelines**: Unix programs are designed to be 'chained' together to form long multi-step pipelines
* **Science Software**: Lots of Scientific software is designed to run in a Unix environment


### User interface

The Unix interface is a text-based command driven one; often known as a Command Line Interface (CLI).  This means
that you control it by issuing (i.e. typing) commands at the command prompt.  Consequently, the Mouse does not perform
any function in the unix environment.

<img src="../media/unix-screen.png" title="Unix Screen" alt="Unix Screen" width="500px" />


### Command prompt

The command prompt is the first thing you see when you connect to a Unix Computer.  Its purpose is to receive the next
command from you, the user.

<img src="../media/cli-parts.png" title="Command prompt" alt="Command prompt" width="500px" />

There are several parts that make up the command prompt:

* **Time**: the time (when the last command finished)
* **Username**: the username that you are logged in as
* **Hostname**: the name of the computer that your are connected to
* **Current working directory**: the current position within the file system that your are working.  More to follow
* **Prompt**: this is simply a sign to the user that the computer is ready to accept the next command

From this point forward in this document, the command prompt will be simply represented as a '$' rather than
the full command prompt as shown above.  When copying and pasting commands you should NOT copy the '$' sign.


### Command line

Below is an example command with various flags and options.

<img src="../media/cli-command.png" title="Command line" alt="Command line" width="500px" />

There are a number of parts which may be included in a command; each is separated by one of more 'white-space' characters (i.e. space, tab):

* **Command**: this is the name of the program (command) that you want to run
* **Flag**: these turn on (or off) specific features in the program.  They consist of a dash (-) followed by
  a single character.
* **Long flag**: same as flag except they are generally two dashes (--) followed by a word (or two)
* **Option**: set the value of a configurable option.  They are a flag (or long flag) followed by a value
* **Anonymous options**: these are one or more options that are specified in the required order
* **Quoted value**: if you need to specify a space (or tab) in an option then you will need to use double (") or
  single (') quotes on each side of the value.


### File system

The file-system of a unix computer can be thought of as an up-side-down tree.  The topmost directory has a special name
called 'root'; it contains all files and directories that are on the computer system.  It is represented by a single
slash (/).  The figure below shows an example file system with directories (black outline boxes) and files (grey dashed
boxes).

<img src="../media/file-system.png" title="File system" alt="File system" width="400px" />

At the top level we have one file (settings) and one directory (home).  Inside the home directory we have two directories
(user1 and user2) and so on.


#### Absolute file names

Absolute file names are a single unique name for each file and directory within the computer.  They start with the slash (/)
character and follow all the parent directories above the file/directory.

<img src="../media/file-system-settings.png" title="File system" alt="File system" width="400px" />

**Absolute file name**: */settings*

<img src="../media/file-system-file01.png" title="File system" alt="File system" width="400px" />

**Absolute file name**: */home/user1/file01.txt*

<img src="../media/file-system-user2.png" title="File system" alt="File system" width="400px" />

**Absolute file name**: */home/user2*  

**Note**: the final slash is not needed (but generally doesn't hurt if it is present).


#### Current working directory

The *current working directory* is the current location within the file system that you are currently using.
When you first login to a unix computer it will begin with the current working directory set to your home directory, that
is, a place that is unique to you and generally nobody else will have access to it.

<img src="../media/unix-cwd.png" title="Current working directory" alt="Current working directory" width="500px" />

Remember from earlier that the current working directory is shown in the command prompt.


#### Relative file names

Relative file names are a short cut to writing file names that are shorter.  The difference between an absolute file
name is that relative file names do NOT begin with a slash.

<img src="../media/file-system-muscle.png" title="File system" alt="File system" width="400px" />

If your current working directory is set to */home* you can leave this part from the beginning of the filename.

**Relative file name**: *user1/muscle.fq*

(Note: the absence of the leading slash)


**Special file names**:

There are a few further short cuts for typing relative file names:

* *~* (Tilde): is a short cut to your home directory
* *.* (dot): is a short cut for the current directory
* *..* (2x dot): means the parent (or one directory up) from current directory
* *...* (3x dot): does not mean anything (a gotcha for new users).  If you want 2 directories up then chain two
double dot's  e.g. *../..*

**Note**: the special file names above can be used within absolute and relative file name and used multiple times.


<img src="../media/file-system-muscle2.png" title="File system" alt="File system" width="400px" />

Now, if the current working directory is changed to */home/user2* the relative path to muscle.fq is different.

**Relative file name**: *../user1/muscle.fq*
