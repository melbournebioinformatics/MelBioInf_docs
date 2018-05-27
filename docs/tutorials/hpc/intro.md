
### What is an HPC?

An HPC is simply a large collection of server-grade computers working together to solve large problems.

* **Big**: HPCs typically have lots of CPUs and Memory and consequently large jobs.
* **Shared**: There are usually lots of users making use of it at one time.
* **Coordinated**: There is a coordinator program to ensure fair-use between its users.
* **Compute Collection**: HPCs use a number of computers at once to solve lots of large jobs.

<img src="../media/drawing1.png" title="HPC Structure" alt="HPC Structure" width="300px" />

**Figure**: The user (face at top) interacts with their local PC/Laptop through the keyboard and screen.  The PC/Laptop will 
connect to the Head/Login node of the HPC interactively.  The Head/Login node will send the jobs off to the Compute Nodes when
one is available. 

### Why use HPCs?

The main reason we use HPCs is because they are quite big.  Given their size, they are usually very expensive, however through 
sharing the resources the per user/job cost can be kept low.

* **Many CPUs**: HPCs typically have hundreds to tens of thousands of CPUs.  Compare this with the 4 or 8 that your PC/Laptop might have.
* **Large Memory**: Hundreds of GBs to multiple TBs of RAM are typical for each node.
* **Efficient use**: Through sharing the resources each user can have access to a very large computer for a period and hand 
it back for others to use later.


### Software Modules

There are typically hundreds to thousands of software packages installed on an HPC.  Given that each can have its own special 
requirements and multiple versions will be made, software on the HPC will most commonly be packaged and only made available 
to you when you request it.

* **Packaged**: to avoid conflicts between software, each is packaged up into a module and only used on demand.
* **Loadable**: before using a software module you need to load it.
* **Versions**: given not all users want to use the same version of software (and to compare new results with old you might 
   need the same version)    each version is made into its own software module so you have ultimate control.

### Job Submission

Job Submission is the process of instructing the HPC to perform a task for you.  Depending on the HPC software installed on 
your HPC, the process of doing so might be different.

* **SLURM**: this workshop uses an HPC that uses the SLURM HPC software.  Some common alternatives (not covered) are PBS or 
   SGE/OGE.
* **Queues (Partition)**: when a job is submitted it is added to a work queue; in SLURM this is called a Partition.
* **Batch**: HPC jobs are not 'interactive'.  By this we mean, you can't type input into your job's programs and you won't 
   immediately see the output that your program prints on the screen. 

#### Resources

So that SLURM knows how to schedule and fit jobs around each other, you need to specify what resources your job will use.
That is, you need to tell it how many CPUs, RAM, Nodes (servers), and Time you need.

1. **CPUs**: most software is limited using 1 CPU by default but many can use more than one (or you can run multiple copies at once).
   The number of CPUs you specify needs to match how many things your software can do at once.
2. **Memory**: you need to estimate (or guess) how much memory (RAM) your program needs.
3. **Nodes**: most software will only use one of the HPC's Nodes (i.e. one server), but some software can make use of more than
   one to solve the problem sooner.
4. **Time**: like when you are scheduling meetings, SLURM needs to know how long each job will take (maximum) so it can organise
   other jobs afterwards.

#### Job Types

There are two types of jobs that you can submit:

1. **Shared**: a shared job (as the name suggests) is one that shares a node with other jobs.  This is the default and preferred method.
2. **Exclusive**: an exclusive job gets a single (or multiple) nodes to itself.  Given this exclusivity, this type of job must know how 
   to use multiple CPUs as most HPCs have at least 16 CPUs per node.

   
