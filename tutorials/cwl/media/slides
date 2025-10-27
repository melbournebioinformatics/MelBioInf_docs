layout: true
class: content
---
# Common Workflow Language for Bioinformatics

.center[
![](images/cwl_logo.svg)
### Michael Milton
]
---

.container[
.row[
.col-10[
<img src="images/melbioinf_logo.png" style="height: 100px;">
]
.col-2[
<img src="images/melbioinf_unimelb.png" style="height: 100px;">
]
]

.row[
.col-8.offset-2[
Providing bioinformatics support for all researchers and students in Melbourne’s biomedical and biosciences precinct.
]
]
<img src="images/melbioinf_capabilities.png" style="width: 100%">
]

---
class: center, middle

.center[
# Part 1: Introduction
.fa-container[
.fas.fa-book.fa-10x[]
]
]
---
## Motivation

* Most bioinformatics involves running many command-line tools; aligners like `bwa` variant callers like `gatk`, and RNA
Seq tools like `cuffdiff`
* However, once this list of tools reaches a certain quantity and complexity, it becomes hard to  reproduce exactly what you ran and with what parameters
* A bash script may help with this, but proper workflows...

    * Run fully parallel to speed up execution
    * Work automatically with batch systems like SLURM
    * Are written declaratively, allowing the system to work out the optimal order of execution for you
    * Save you having to hard-code input parameters and temporary files
    * Are much more readable than bash
---
## Housekeeping
* Slides
    * These slides are hosted at <https://tinyurl.com/ycenoxbf>
    * I recommend you open them in your browser so you can follow along
* Downloads for Part 1:
    * Rabix Composer:
        * A graphical CWL editor
        * **Absolutely essential for the first part of this workshop**
        * <https://github.com/rabix/composer/releases>
    * Test Data
        * The data we will be using to test our workflows is located [here](data/data.tar.gz).
        * Please download it now so you have a copy for later
    * Docker
        * An engine for running tools inside containers
        * <https://store.docker.com/search?type=edition&offering=community>
---
## CWL Structure
* Tools .fas.fa-wrench[]
    * Are wrappers that describe to the CWL engine how a command-line tool works
    * List its inputs and outputs, and the command to produce one from the other
    * Aren't workflow-specific, so some already exist for commonly used tools
* Workflows .fas.fa-share-alt[]
    * Explain how tools are connected to each other and in what order
    * Are generally project-specific
    * Can be nested inside each other
---

class: center, middle

.center[
# Part 2: Tools
.fa-container[
.fas.fa-wrench.fa-10x[]
]
]

---
## Obtaining Tool Definitions

There are a few useful sources of CWL tool definitions:
* Dockstore
    * <https://dockstore.org>
    * Dockstore - a database of CWL and WDL workflows and tools
    * Once you find a tool definition you like, click:
        * "Files" → "Descriptor Files" → Download
    * Unfortunately many of these definitions are out of date...

.center[
![](images/dockstore_circled.png)
]
---
## Obtaining Tool Definitions

There are a few useful sources of CWL tool definitions:
* Official CWL Workflows Repository
    * <https://github.com/common-workflow-language/workflows>
    * If you clone the entire repo, it will provide a useful library of tool definitions:
    ```bash
    git clone https://github.com/common-workflow-language/workflows
    ```
---
## Obtaining Tool Definitions
.alert.alert-primary[
.alert-heading[
### Exercise
]
.row[
.col-8[
* Clone the CWL Workflows repo
* Open the repo in the Rabix Composer
* Find and open the tool definition for `bwa mem`
* Delete this line:
    ```yaml
    - $import: bwa-docker.yml
    ```
    * Although this is valid CWL, unfortunately Rabix doesn't support import statements
* Open the visual editor
![](images/visual_editor.png)
]
.col-4[
![](images/rabix_open_project.png)
]
]
]
--
.alert.alert-success[
.alert-heading[
### Answer
]
* [`bwa-mem.cwl`](cwl/bwa-mem.cwl)
]


---
## Rabix Tool Overview
.center[
<img src="images/rabix_bwa_mem_overview.png" style="width: 600px">
]
---
## Rabix Tool Overview
* **Docker Image** defines an optional image in which to run this tool
* **The Base Command** defines what command line tool to run
    ![](images/bwa_basecommand.png)
* **Arguments** are confusingly named - they won't be needed for a while. Don't confuse this with **Input Ports**
* The **Command Line** at the bottom builds up a bash command that should correspond to the tool you're wrapping
---
## Rabix Overview
* **The Input Ports** define the input files and flags
    ![](images/bwa_inputs.png)
* **The Output Ports** define the output files
    ![](images/bwa_outputs.png)

---
## Wrapping Samtools
.alert.alert-primary[
.alert-heading[
### Exercise
]
Follow along with the instructions to make a tool wrapper for `samtools sort`
]
---
## Wrapping Samtools
1\. Start by making a new tool definition in Rabix

.center[
![](images/rabix_new_tool.png)
]
---
## Wrapping Samtools
2\. Name it after the tool you're wrapping

.center[
![](images/rabix_samtools_name.png)
]
---
## Wrapping Samtools
3\. Add the "base command" - the fixed part of the command that will never change

* Make sure that each part of the command is on a new line
.center[
![](images/rabix_samtools_base.png)
]
---
## Wrapping Samtools

4\. Define the inputs(s)

.row[
.col-8[
Each input has:
* A name (`ID`)
* A type:
    * `boolean`, `int`, `long`, `float`, `double`, `string`, `File`, or `Directory` for single values
    * An array of the above types, e.g. `File[]`
* Is either optional or required
* A way for this input to be used on the command line. Either:
    * With a numerical `Position` for positional arguments like `gzip some_file`. Here the position is 0 (the first positional argument)
    * With a `Prefix` for named arguments, e.g. `java -jar something.jar`. Here the input prefix is `-jar`

In this case, we want a required, `File`-type, positional input

]
.col-4[
![](images/rabix_samtools_input.png)
]
]
---
## Wrapping Samtools

5\. Define the output(s)

.row[
.col-8[
Each output has:
* A name (`ID`)
* A type (same as an input)
* A glob/filepath which indicates how to find this file. In this case, the file is located wherever we piped stdout

In this case, we want a `File` type output with a fixed filename (`sorted_alignment.bam`)

]
.col-4[
![](images/rabix_samtools_output.png)
]
]
---
## Wrapping Samtools

6\. If the command produces output from stdout, you need to pipe it to a file so that it can be picked up by the output
glob

.center[
![](images/rabix_samtools_other.png)
]
---
## Wrapping Samtools

7\. If everything has been done correctly, you should have the following command line at the bottom of the page:
```bash
samtools sort /path/to/input.ext > sorted_alignment.bam
```

8\. Save your tool definition

![](images/rabix_save.png)

--

.alert.alert-success[
.alert-heading[
### Answer
]
* [`samtools-sort.cwl`](cwl/samtools-sort.cwl)
]

---
## Wrapping Freebayes
.alert.alert-primary[
.alert-heading[
### Exercise
]
* [Freebayes](https://github.com/ekg/freebayes) is a variant caller, which takes a BAM alignment and determines how this alignment differs from the "normal"
reference genome
* Freebayes is called on the command-line as follows:
    ```bash
    freebayes --fasta-reference h.sapiens.fasta NA20504.bam
    ```
* Write a new CWL tool wrapper for Freebayes that supports this command
* If done correctly, your `Command Line` section should show something like:
```bash
freebayes /path/to/input.ext --fasta-reference /path/to/input.ext > variants.vcf
```
]
--
.alert.alert-success[
.alert-heading[
### Answer
]
* [`freebayes.cwl`](cwl/freebayes.cwl)
]

---
## Docker
* We have given CWL instructions on how to *run* these tools, but not how to *get* these tools. For this we can use Docker
* Docker images are tiny virtual machines that have applications pre-installed inside of them
* You can find docker images of bioinformatics tools:
    * In [Biocontainers](https://biocontainers.pro/registry/)
    * On [Bioconda](https://bioconda.github.io/recipes.html)
    * On [Docker Store](https://store.docker.com/)
* These sites will show a `docker pull` commmand in the form `docker pull username/imagename`. Just take the `username/imagename` section
* Once you've found a Docker image, you can plug it into the "Docker Image" section in Rabix:

    ![](images/docker_container_section.png)

---
## Bioconda
* To find Bioconda images
    * Open the [Bioconda search page](https://bioconda.github.io/recipes.html)
    * Search for the tool you want
    * Click any version of that tool
    * On the tool page, follow the link to its container
    * Click the Tag button on the left
    * Pick a version you want, and, click the Download button next to it
    * Choose the Docker Pull (by tag) option

---
## Updating our tool to use Docker
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Find an appropriate Docker image for `bwa`, `samtools sort` and `freebayes`, using Biocontainers
* Once you have found the right images, plug them into the "Docker Image" section
]

--

.alert.alert-success[
.alert-heading[
### Answers
]
* Here are some completed tool definitions:
    * [`bwa.cwl`](cwl/bwa-mem.cwl)
    * [`samtools-sort.cwl`](cwl/samtools-sort.cwl)
    * [`freebayes.cwl`](cwl/freebayes.cwl)
]


Now that we have a way to find the actual tools, we can start actually running CWL...

---
## Running Tools
.alert.alert-primary[
.alert-heading[
### Exercise
]

* Open the Test tab for the BWA tool
* Double click each input bubble and browse for the appropriate file we downloaded at the start of this workshop*
* Double click the main bubble and set the output filename
* Click "Run"
![](images/run_button.png)
* Find your output file by clicking the "Output Directory" button

.center[
<img style="height: 80%" src="images/rabix_bwa_test.png">
]
\\* if this doesn't work you may have to specify secondary files for the reference
]


---
## Secondary Files
* Some files, like indexes, are never considered a main file, but are instead designed to accompany another file,
for example:
    * `.bai` files accompany `bam` alignments
    * `.tbi` indices which accompany `vcf` variant calls.
* Secondary files can accompany both input and output files
* The simplest way to specify secondary files is as a string that will be appended onto the main file
* For example `.bai` means `main_file.bam.bai` is the secondary file name
---
## Secondary Files
* For example, if we open the `reference` Input Port from the `bwa-mem` tool, we will see:

.center[
![](images/bwa_secondary_files.png)
]
* This means that the genome reference file must be accompanied by an `.amb`, `.ann`, `.btw`, `.pac`, and `.sa` index

---
# Adding secondary files
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Freebayes requires an indexed BAM file (it has a `.bai` secondary file)
* Freebayes also requires that the reference genome has a `.fai` index
* Edit your Freebayes tool definition to include these indexes
]
--
.alert.alert-success[
.alert-heading[
### Answers
]
* [`freebayes.cwl`](cwl/freebayes.cwl)
]


---
## Dynamic Expressions
### Motivation
* Sometimes, some of the values in our CWL need to be calculated dynamically
* A common situation where we need dynamic values is when a command creates an output file whose name is based on the input file
* For example, `gzip file.txt` produces `file.txt.gz`
* In order to do this, we can embed an expression in some CWL fields

---
## Dynamic Expressions
### Parameter References
* The most basic expressions you can use are called [Parameter References](https://www.commonwl.org/v1.0/CommandLineTool.html#Parameter_references)
* These consist of an expression in this form: `$(...)`
* For instance, `$(inputs.extractfile)`
* These use a subset of JavaScript syntax that allows for property access and nothing else
* You can use multiple `$(...)` expressions within the one string
---
## Dynamic Expressions
### Variables
* Within the scope of an expression, the following variables are available:
    * `inputs`: a dictionary of all the inputs provided to this tool. If these inputs have the type
        [`File`](https://www.commonwl.org/v1.0/CommandLineTool.html#File) or
        [`Dictionary`](https://www.commonwl.org/v1.0/CommandLineTool.html#Directory), they have special properties
    * `self`: value depends on the specific context of this expression. e.g. when used in `secondaryFiles`, `self`
        is set to the main input or output, so that the function can calculate the secondary files for this file
    * `runtime`:
        * `runtime.outdir`: an absolute path to the designated output directory
        * `runtime.tmpdir`: an absolute path to the designated temporary directory
        * `runtime.cores`: number of CPU cores reserved for the tool process
        * `runtime.ram`: amount of RAM in mebibytes (2**20) reserved for the tool process
        * `runtime.outdirSize`: reserved storage space available in the designated output directory
        * `runtime.tmpdirSize`: reserved storage space available in the designated temporary directory
---
## Dynamic Expressions
### JavaScript Expressions (advanced)
* If you enable JavaScript expressions with the following line, you can use JavaScript in your expressions for even more
power:
    ```yaml
    requirements:
      - class: InlineJavascriptRequirement
    ```
* The `$(...)` syntax is upgraded into a JavaScript expression, which has all the same abilities as a parameter reference,
but can now calculate values:
    * e.g.
        ```javascript
        $(1 + 2)
        ```
* `${...}` denotes a scope in which functions, variables, loops etc may be used, and a `return` statement is needed
    to return the final value
    * e.g.
        ```javascript
        ${
            var x = 1;
            var y = 2;
            return x + y;
        }
        ```
* Refer to the [expressions](https://www.commonwl.org/v1.0/CommandLineTool.html#Expressions) section of the CWL spec
---
## Wrapping Samtools Index

.alert.alert-primary[
.alert-heading[
### Exercise
]
* Part 1
    * Samtools index is invoked as follows:
    ```bash
    samtools index aln.bam
    ```
    * The output from `samtools index` will be a BAM file, with a `.bai` index as a secondary file
    * Use what you have learned from wrapping `bwa` to make a wrapper for the `samtools index` subcommand, including
    base command, inputs, outputs and docker image
    * But what is the output `glob`...?
]
---
## Wrapping Samtools Index

.alert.alert-primary[
.alert-heading[
### Exercise
]
* Part 2
    * The tricky part is that the output BAM file is actually the same as the input file
    * This means your output file Glob will use an expression to use the input file as the output
    * Hint: `$(input.INPUTNAME.basename)` will return the full path to the input named `INPUTNAME`
    * In addition, to grant access to this file, you need to add `$(input.INPUTNAME)` as an **expression** (not a file)
    in the File Requirements section
]

--

.alert.alert-success[
.alert-heading[
### Answer
]
* [`samtools-index.cwl`](cwl/samtools-index.cwl)
]

---
class: center, middle

.center[
# Part 3: Writing Workflows
.fa-container[
.fas.fa-share-alt-square.fa-10x[]
]
]
---
## Workflows in the Rabix Composer

* Refresher: workflows define how your tools connect to each other to form a data flow
* In Rabix, you add tools to your workflow by dragging and dropping from the sidebar
.center[
![](images/drag_tool.png)
]
* To connect the tools, you then drag a line between input ports and output ports
.center[
<video controls autoplay loop>
    <source src="images/connect.mp4" type="video/mp4">
</video>
]
---
* To specify an input that the user must provide, drag it onto empty space
.center[
<video controls autoplay loop>
    <source src="images/rabix_workflow_inputs2.mp4" width="518" height="282" type="video/mp4">
</video>
]
* To specify an output that is not used by another tool, drag it onto empty space
.center[
<video controls autoplay loop>
    <source src="images/rabix_workflow_output2.mp4" width="518" height="282" type="video/mp4">
</video>
]
---
## A Variant Calling Pipeline in Rabix
.alert.alert-primary[
.alert-heading[
### Exercise
]
* In the "Visual Editor" tab, make a basic workflow that connects:
    * `bwa` → `samtools sort` → `samtools index` → `freebayes`
* Set a value for the `output_filename` on BWA
* Now, change to the "Test" tab, and run that workflow with the all of the test data provided
* This should work the same as when we ran BWA
]
---
## A Variant Calling Pipeline in Rabix
.alert.alert-success[
.alert-heading[
### Answer
]
* You should have workflow a bit like this:
<img src="images/germline_workflow.png" style="width: 100%">
* [`germline_workflow.cwl`](cwl/germline_workflow.cwl)
]
---
## Subworkflows
* Often when making workflows you will want to encapsulate a whole pipeline of tools and re-use that together
as a single unit, for example "alignment"
* To do this, simply make a workflow that encapsulates this functionality, and use as a step in another workflow
* In the Rabix composer, you can drag workflows onto a workflow in the same way you add tools, however they have this
logo:
.fas.fa-share-alt.fa-rotate-180[]
.center[
![](images/subworkflow.png)
]

---
# A Tumour-Normal Variant Caller
.alert.alert-primary[
.alert-heading[
### Exercise
]
* First, make an alignment workflow that connects:
    * `bwa` → `samtools sort` → `samtools index`
* Download a Somatic Sniper tool definition [from here](cwl/somatic-sniper.cwl)
* Next, make a tumour-normal workflow that uses this workflow to align both the tumour and the normal reads, and then
feeds the alignments into Somatic Sniper:
    * `Tumour Alignment` → `Somatic Sniper` ← `Normal Alignment`
]
---
# A Tumour-Normal Variant Caller
.alert.alert-success[
.alert-heading[
### Answer
]
* Your alignment workflow should look like this
<img src="images/alignment_workflow.png" style="width: 100%">
* [`alignment_workflow.cwl`](cwl/alignment_workflow.cwl)
]
---
# A Tumour-Normal Variant Caller
.alert.alert-success[
.alert-heading[
### Answer
]
* And your somatic workflow should look like this:
<img src="images/somatic_workflow.png" style="height: 90%">
* [`somatic_workflow.cwl`](cwl/somatic_workflow.cwl)
]
---
class: center, middle

.center[
# Part 4: YAML
.fa-container[
.fas.fa-keyboard.fa-10x[]
]
]
---
## Downloads for Part 4:
* Python:
    * Language runtime used by many CWL executors
    * <https://www.python.org/downloads/>
* cwltool:
    * A simple CWL executor
    * <https://github.com/common-workflow-language/cwltool#install>
* Some text editor
    * Your system editor is probably fine
    * If not, I recommend Atom: <https://atom.io/>
---
## CWLTool
* `cwltool` is the reference CWL runner
* It can do many things with CWL, including validate and run it on the command line
* To validate a tool or workflow definition, open a terminal and type:
    ```bash
    cwltool --validate /path/to/tool.cwl
    ```
* To run a tool or workflow, use:
    ```bash
    cwltool /path/to/tool.cwl --input-name value --another-input value
    ```
* This will be helpful later on!
---
## YAML Format

* CWL tools and workflows are written in a format called `YAML`
* YAML is a lot like JSON, however:
    * Structures like dictionaries (maps) and arrays (lists) don't have delimiters like `{}` and `[]`. Instead, their type is
    implied by their contents (`:` and `-`)
    * Whitespace is used to indicate when an object is nested inside another

.row[
.col-6[
#### YAML
```yaml
a_dictionary:
  key: value
  another_key: another value

an_array:
  - value
  - another value

an_integer: 3

a_string: this is a string
```
]
.col-6[
#### JSON
```json
{
  "a_dictionary": {
    "key": "value",
    "another_key": "another value"
  },
  "an_array": [
    "value",
    "another value"
  ],
  "an_integer": 3,
  "a_string": "this is a string"
}
```
]
]
---
## Tool YAML
```yaml
class: CommandLineTool

cwlVersion: v1.0

baseCommand:
  - samtools
  - sort

inputs:
  - id: alignment
    type: File
    inputBinding:
      position: 0

outputs:
  - id: sorted_alignment
    type: File
    outputBinding:
      glob: sorted_alignment.bam

stdout: sorted_output.bam
```
---
## Tool YAML

* This is a tool (as opposed to a workflow)
    ```yaml
    class: CommandLineTool
    ```
* This follows version 1.0 of the CWL standard:
    ```yaml
    cwlVersion: v1.0
    ```
* The base command is `samtools sort`:
    ```yaml
    baseCommand:
      - samtools
      - sort
      ```
* The docker image
    ```yaml
    requirements:
      - class: DockerRequirement
        dockerPull: biocontainers/samtools
    ```
---
## Tool YAML

.row[
.col-6[
The inputs:
```yaml
inputs:
  - id: alignment
    type: File
    inputBinding:
      position: 0
  - id: sort_by_name
    type: boolean?
    inputBinding:
      prefix: -n
```
* Inputs is a YAML array
* Each entry has a name (`id`), and a type, like in the Rabix Composer
* However, fields to do with command-line binding are put into the `inputBinding` dictionary. This includes `position`
and `prefix`
]
.col-6[
The outputs (including stdout):
```yaml
outputs:
  - id: sorted_alignment
    type: File
    outputBinding:
      glob: sorted_alignment.bam

stdout: sorted_output.bam
```
* Each entry has a name (`id`), and a type, like in the Rabix Composer
* The `glob` key is under `outputBinding`
* `stdout` is its own section, there is no `other` section like in the Rabix Composer
]
]
---
## Dictionary vs Array Syntax
* A number of sections in CWL can be written either as an array, or as a dictionary
* These include:
    * Inputs in a tool or workflow
    * Outputs in a tool or workflow
    * Steps in a workflow
* Rabix uses the array syntax by default, but many tutorials will use the dictionary syntax, as it is generally more
concise

.row[
.col-6[
* I have shown you the array syntax, which looks like this:
    ```yaml
    inputs:
      - id: alignment
        type: File
        inputBinding:
          position: 0
      - id: sort_by_name
        type: boolean?
        inputBinding:
          prefix: -n
    ```
]
.col-6[
* If you turn this into a dictionary, with the `id` fields as the keys, this input list can be written as:
    ```yaml
    inputs:
      alignment:
        type: File
        inputBinding:
          position: 0
      sort_by_name:
        type: boolean?
        inputBinding:
          prefix: -n
    ```
]
]
---
## Stdout Type
* So far when we've wanted to capture the stdout of a tool, we have piped it to a file and then captured that output file
* CWL has a much neater way of doing this - you simply set an output's type to `stdout`
* This simplifies this code:
    ```yaml
    outputs:
      an_output_name:
        type: File
        outputBinding:
          glob: a_stdout_file

    stdout: a_stdout_file
    ```
* To this
    ```yaml
    outputs:
      an_output_name:
        type: stdout
    ```
* In addition, this allows CWL to stream outputs between stages, instead of saving them to the filesystem first
---
## Using the `stdout` Type
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Edit your `samtools sort` and `freebayes` tools to use the `stdout` type
* Edit these files using your text editor only, without the Rabix Composer
* Once you've finished, test your files with `cwltool --validate`
]
---
## Using the `stdout` Type
.alert.alert-success[
.alert-heading[
### Answer
]
* Your answer for `samtools sort` should look a bit like this:
    ```yaml
    class: CommandLineTool
    cwlVersion: v1.0

    baseCommand:
      - samtools
      - sort

    inputs:
      - id: alignment
        type: File
        inputBinding:
          position: 0

    outputs:
      - id: sorted_alignment
        type: stdout
    ```
* [`samtools-sort-stdout.cwl`](cwl/samtools-sort-stdout.cwl)
]
---
## Wrapping Cutadapt (by hand!)
.alert.alert-primary[
.alert-heading[
### Exercise
]
* We have a working germline pipeline, but it doesn't filter out bad quality reads!
* For this we need to add a read trimmer like `cutadapt`
* Cutadapt is invoked like this (`q` is the quality cutoff)
   ```bash
   cutadapt -q 10 input.fastq > output.fastq`
   ```
* Write a tool definition for `cutadapt` by hand
]
---
## Wrapping Cutadapt (by hand!)
.alert.alert-success[
.alert-heading[
### Answer
* [`cutadapt.cwl`](cwl/cutadapt.cwl)
]
.row[
.col-6[
```yaml
class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - cutadapt

requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/cutadapt:1.16--py36_1
```
]
.col-6[
```yaml
inputs:
  quality_cutoff:
    type: int
    inputBinding:
      prefix: -q

  input_fastq:
    type: File
    inputBinding:
      position: 0

outputs:
  trimmed_fastq:
    type: stdout
```
]
]
]
---
## Arguments
* The arguments section does actually have a use!
* It's used for fixed arguments that need more configuration than the `baseCommand` section provides
* `arguments` is a list of `inputBinding` objects, exactly the same as if they were inside an `inputs` entry
* For example, if we wanted `bwa` to always use a fixed output name (since it doesn't matter anyway), we could remove this
    ```yaml
    inputs
      - id: output_filename
        type: string
    ```
* And replace it with this
    ```yaml
    arguments:
      - position: 0
        valueFrom: alignment.sam
    ```
---
## Updating cutadapt
.alert.alert-primary[
.alert-heading[
### Exercise
]
* Cutadapt is invoked like this for paired reads (`q` is a quality cutoff):
    ```bash
    cutadapt -q 10 -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq
    ```
* Using a text editor, edit your tool definition for `cutadapt`
* `-o out.1.fastq` and `-p out.2.fastq` can be `arguments`
* The output needs to be an array of files, and for this you can use a glob such as `out*fastq`
]
---
## Updating cutadapt
.alert.alert-success[
.alert-heading[
### Answer
]

* [`cutadapt-paired.cwl`](cwl/cutadapt-paired.cwl)

.row[
.col-6[
```yaml
class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - cutadapt

arguments:
  - prefix: '-o'
    valueFrom: out.1.fastq
  - prefix: '-p'
    valueFrom: out.2.fastq

requirements:
  - class: DockerRequirement
    dockerPull: quay.io/biocontainers/cutadapt:1.16--py36_1
```
]
.col-6[
```yaml
inputs:
  quality_cutoff:
    type: int
    inputBinding:
      prefix: -q

  input_fastq:
    type: File[]
    inputBinding:
      position: 0

outputs:
  trimmed_reads:
    type: File[]
    outputBinding:
      glob: 'out.*.fastq'

```
]
]
]
---
## Workflows in CWL Files
* Like CWL tools, workflows are represented in a YAML format

.row[
.col-sm[
```yaml
cwlVersion: v1.0
class: Workflow
inputs:
  inp: File
  ex: string

outputs:
  classout:
    type: File
    outputSource: compile/classfile
```
]
.col-sm[
```yaml
steps:
  untar:
    run: tar-param.cwl
    in:
      tarfile: inp
      extractfile: ex
    out: [example_out]

  compile:
    run: arguments.cwl
    in:
      src: untar/example_out
    out: [classfile]
```
]
]
---
* All workflows have the class "Workflow"
    ```yml
    class: Workflow
    ```
* The workflow needs a list of inputs, which will then feed into the tools:
    ```yaml
    inputs:
      inp: File
      ex: string
    ```

* We also pick certain outputs from the tools in this workflow to use as workflow outputs.
    For this, we use the format `outputSource: step_name/output_name`
    ```yaml
    outputs:
      classout:
        type: File
        outputSource: compile/classfile
    ```
---
* The `steps` section is a dictionary of `step_name: step_body` which contains the body of the workflow
* Each step has:
    * A name:
    ```yaml
      compile:
    ```

    * A workflow or tool to run for this stage
    ```yaml
        run: arguments.cwl
    ```

    * A list of inputs with a source for each.
    ```yaml
        in:
          src: untar/example_out
    ```
    In this case, the `src` parameter is provided by the `example_out` parameter of the step named `untar`

    * A list of outputs we intend to use in the wider workflow
    ```yaml
        out: [classfile]
    ```
---
.alert.alert-primary[
.alert-heading[
## Exercise - Adding Cutadapt
]
* Using a text editor, copy or edit your germline workflow, and add the cutadapt tool
    * `cutadapt` will have to be a new stage in the workflow
    * The `bwa` stage will have to be changed to take its inputs from `cutadapt`
* Once you're finished, try to run the workflow with:
    ```bash
    cwltool cwl/germline_workflow_cutadapt.cwl --reference wildtype.fna --input_fastq mutant_R1.fastq --input_fastq mutant_R2.fastq
    ```
* You may have to add the secondary files to the workflow's inputs before it will let you run this
    ```yaml
    inputs:
      - id: reference
        type: File
        secondaryFiles:
          - .amb
          - .ann
          - .bwt
          - .pac
          - .sa
    ```

]
---
.alert.alert-success[
.alert-heading[
### Answer
]
* [`germline_workflow_cutadapt.cwl`](cwl/germline_workflow_cutadapt.cwl)
]
---
.alert.alert-primary[
.alert-heading[
## Final Exercise - Manual Workflow
]
* Rabix really butchered our workflow YAML
* Try to re-write the germline workflow from scratch
* This should now connect:
    * `cutadapt` → `bwa` → `samtools sort` → `samtools index` → `freebayes`
* Test it once you finish using `cwltool`
]
---
class: center, middle

.center[
# That's All!
.fa-container[
.fas.fa-graduation-cap.fa-10x[]
]
]
