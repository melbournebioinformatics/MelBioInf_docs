# Workshop test

## Introduction
This document describes a simple series of steps to test a workshop. This test can be used to test for (1) expected behaviour and (2) multi-user performance on the workshop computer (e.g. Nectar virtual computer).

There are two main components to this system.

`run_complete_workshop.sh` contains an ordered series of commands to be run by the participants directly from the workshop. It is expected that the environment and tools are setup appropriately to simulate the workshop environment.

`test_workshop.sh` script performs the test by simulating a workshop environment. Currently, this tool takes in a variable `N` as the expected number of participants to simulate on a single machine (this is the machine where this tool in run). Next, the tool run the `run_complete_workshop.sh` script to simulate the workshop `N` times in parallel, each time creating a unique workshop working directory (prefix `MYTEST_`). At the moment this tool does not have test to check for expected files and outcomes. Requires GNU Parallel.

```bash
!/bin/bash

N=5
seq -w 1 "$N" | parallel --joblog test.log -j 0 --halt 2 'RUNNUMBER={} ./run_complete_workshop.sh'

```


## How to (Variant Calling GATK1)
* Step 1 : Log in to the virtual computer (e.g. Nectar) using the provided credentials
* Step 2: Download the required files to the home directory. Run the `cd` command to make sure you are currently in the home directory.
    * `wget https://github.com/melbournebioinformatics/MelBioInf_docs/raw/master/docs/tutorials/variant_calling_gatk1/files/run_complete_workshop.sh`
    * `wget https://github.com/melbournebioinformatics/MelBioInf_docs/raw/master/docs/tutorials/variant_calling_gatk1/files/test_workshop.sh`
* Run the command below to test:
```bash
# change file permissions
chmod +x run_complete_workshop.sh
chmod +x test_workshop.sh
# run test
./test_workshop.sh
```

* Step 3: After the successful completion of the test, you should see a log file `test_workshop.log`. This file details the time it took to run `run_complete_workshop` script in parallel `N` times.

    ```bash
    Seq	Host	Starttime	JobRuntime	Send	Receive	Exitval	Signal	Command
    3	:	1630852205.612	  1047.621	0	306	0	0	RUNNUMBER=3 ./run_complete_workshop.sh
    1	:	1630852205.608	  1048.805	0	306	0	0	RUNNUMBER=1 ./run_complete_workshop.sh
    4	:	1630852205.614	  1048.803	0	306	0	0	RUNNUMBER=4 ./run_complete_workshop.sh
    2	:	1630852205.610	  1050.313	0	306	0	0	RUNNUMBER=2 ./run_complete_workshop.sh
    5	:	1630852205.617	  1050.391	0	306	0	0	RUNNUMBER=5 ./run_complete_workshop.sh
    ```
    The output from each parallel running of the workshop material will be in its corresponding directories (see blow).

    ```bash
    ls -1

    MYTEST_1
    MYTEST_2
    MYTEST_3
    MYTEST_4
    MYTEST_5
    ```

## Results (Variant Calling GATK1)
The workshop performs as expected for `N=5` users on the current data sets and computer setup.

```bash
16 cpus
32 GB memory
```


Notes: Testing with more users may raise memory issues when loading reference data.
