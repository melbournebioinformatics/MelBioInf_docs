# Authors:

- Bernie Pope, Melbourne Bioinformatics (formerly VLSCI)

- Catherine de Burgh-Day, Dept. of Physics, The University of Melbourne

# General information

- Python modules are stored in files containing a ".py" suffix (e.g solver.py).

- The main implementation of Python is called CPython (it is written in C). It is byte-code interpreted.

- Python can be used in two modes: interactive and scripted. In interactive mode you enter a program fragment and Python evaluates
  it immediately and then prints the result before prompting for a
  new input. The interactive prompt is usually rendered as the
  *chevron* `>>>`. In scripted mode your program is stored in one
  or more files which are executed as one monolithic entity. Such
  programs behave like ordinary applications.

- Python has automatic memory management (via garbage collection). Memory is allocated automatically as needed and freed
  automatically when no longer used.

## Python 2 versus Python 3

Currently there are two distinct flavours of Python available:

- Python 2 (2.7.10 at the time of writing)

- Python 3 (3.4.3 at the time of writing)

Python 3 is the new and improved version of the language. Python 3 is
not entirely backwards compatible, but the two versions share much in
common. Version 2 is now in maintenance mode; new features will only be
added to version 3. The public transition from 2 to 3 has been slower
than some people would like. You are encouraged to use version 3 where
possible. These notes are generally compatible with both versions, but
we will point out key differences where necessary.

## Indentation for grouping code blocks

- Python uses indentation to group code blocks. Most other languages use some kind of brackets for grouping.

- The recommended style is to use 4 space characters for a single indent (thus 8 spaces for two indents and so forth).

- You are encouraged **not** to use tabs for indentation because there is no standard width for a tab.

- Most good text editors can be configured so that that tab key is rendered as 4 space characters when editing Python code.

## Style Guide

A popular style guide for Python is known as [PEP 0008](https://www.python.org/dev/peps/pep-0008/), there is a
corresponding tool called [pep8](https://pypi.python.org/pypi/pep8) which will check your code against
the guide and report any transgressions.

**Example, Python compared to C:**

Python program for computing factorial:

```python
# Compute factorial of n,
# assuming n >= 0

def factorial(n):
  result = 1
  while n > 0:
      result *= n
      n -= 1
  return result

print(factorial(10))
```

C program for computing factorial:

```c
#include <stdio.h>

/* Compute factorial of n,
  assuming n >= 0 */

int factorial(int n) {
   int result = 1;

   while (n > 0) {
      result *= n;
      n -= 1;
   }
   return result;
}

int main(void) {
   printf("%d\n", factorial(10));
}
```

Things to note:

- The difference in commenting style.

- C programs are statically typed, and you must declare the type of functions and variables. Python is dynamically typed.

-   Code blocks in C are grouped by braces { }; Python uses indentation for grouping.

-   The C program must have a main function. Python does not require a
    main function, it just executes the top-level statements of the
    module.

-   The result returned by the C function is limited to the size of a
    machine integer (say 32 bits). However, the result returned by the
    Python function is unlimited in its size - it can compute
    arbitrarily large factorials (up to the limit of the available
    memory in your computer).

## Comments

Program comments start with a hash character "\#" and continue until the
end of the line. There are no multi-line comment markers, but that can
sometimes be faked with multi-line string literals.

**Examples:**

```python
# This is a comment.
# This is another comment.
x = 5 # This is a comment that follows some code.
'''This is
a multi-line
string literal
which can sometimes act like
a
comment.
'''
```

## Running a Python program

There are many ways to run Python code:

-   You can run the interpreter in interactive mode. On Unix (Linux, OS
    X) you can run the python command at the command line.

-   If you have Python code stores in a file, say example.py, you can
    run it from the command line like so: python example.py

-   You can use one of several integrated programming environments.
    Python ships with a fairly minimal one called [IDLE](https://docs.python.org/2/library/idle.html), though
    many scientists prefer the more comprehensive [IPython](http://ipython.org/).

-   If your Python code was installed as a package (see below), then it
    may be executed like an ordinary application without the user
    being aware of how the program was implemented.

## Objects and types

-   Every value in Python is an *object* (including functions!).

-   Objects can have attributes and methods, which are accessed via the
    dot "." operator.

-   All objects have a type.

-   Types are also objects!

-   Python is dynamically typed: you may get type errors at runtime but
    never at compile time.

-   `type(x)` returns the type of x.

-   Python variables may be assigned to values of different types at
    different points in the program.

**Interactive examples (Python 3):**

```
>>> # Create a list, assign to the variable x
>>> x = [3, 1, 2, 3]
>>> # Ask for the type of the value assigned to x
>>> type(x)
<class 'list'>
>>> # Ask for the type of the first item in the list (an integer)
>>> type(x[0])
<class 'int'>
>>> # Count the number of times 3 appears in the list
>>> # by calling the count method
>>> x.count(3)
2
>>> # Sort the contents of the list in-place.
>>> # Note that this mutates the list object!
>>> # Also note that Python does not print the result in this case.
>>> x.sort()
>>> # Ask Python to show the value of the list
>>> # assigned to the variable x (note it is now sorted)
>>> x
[1, 2, 3, 3]
>>> # Assign x to an object of a different type (a float)
>>> x = 3.142
>>> type(x)
<class 'float'>
```

# Booleans

-   Represent truth values
-   Values: `True`, `False`
-   Type: `bool`
-   Operators: `and`, `or`, `not`
-   `bool(x)` will convert x to a boolean. The heuristic is that empty things and zero-ish things are `False`, everything else is `True` (but the user can override for their own types).
    -   `False` values:
        -   `False`
        -   `0` (zero integer)
        -   `0.0` (zero float)
        -   `{}` (empty dictionary)
        -   `()` (empty tuple)
        -   `[]` (empty list)
        -   `''` (empty string)
        -   `None`
    -   `True` values:
        -   everything else
-   In numerical contexts `True` is considered equal to the integer `1` and
    `False` is considered equal to the integer `0`. However, these
    conversions are a common cause of bugs and should be avoided.
-   Python will automatically test the *truthiness* of a value if it
    appears in a boolean context.

**Interactive examples:**

```python
>>> not True
False
>>> not False
True
>>> not ()
True
>>> not [1,2,3]
False
>>> True and False
False
>>> True and ()
()
```

# Conditional Statements

-   Conditional statements use the keywords: `if`, `elif`, `else`. The syntax
    for a conditional statement is:

```python
if expression:
    statement-block
elif expression:
    statement-block
...
else:
    statement-block
```

-   A conditional statement must have exactly one `if` part. It may have
    zero or more `elif` parts, and a single optional `else` part at the
    end.

-   The `if` and `elif` parts test the value of their boolean expressions.
    If the expression evaluates to something which is `True` or can be
    converted to `True` (see the rules for Booleans above) then the
    statement block immediately beneath that part is executed.
    Otherwise the following condition (if any) is tried. The `else`
    part, if it exists, is always and only executed if no preceding
    condition was `True`.

**Interactive examples:**

```
>>> if []:
...     print("Was considered True")
... else:
...     print("Was considered False")
...
Was considered False
```

# Numbers and basic mathematics

## Integers

-   Represent whole negative and positive numbers (and zero).

-   The range of integer values is unbounded (up to some limit defined
    by how much memory you have on your computer).

-   Python 2 distinguishes between two integer types `int` and `long`, and
    automatically promotes `int` to long where necessary, whereas Python
    3 considers them all one type called `int`.

-   Base ten is the default literal notation: `42` (means `(4 * 10) + 2`)

-   Hexadecimal literals start with `0x`, octal literals start with `0o`,
    binary literals start with `0b`.

-   `int(x)` will try to convert x to an integer, x can be another numeric
    type (including booleans) or a string. You may specify an optional
    base for the conversion.

**Interactive examples (in Python 3):**

```
>>> 2 ** 200
1606938044258990275541962092341162602522202993782792835301376
>>> 0x10
16
>>> 0b10
2
>>> -0 == 0
True
>>> int("123")
123
>>> int("3.142")
Traceback (most recent call last):
File "<stdin>", line 1, in <module>
    ValueError: invalid literal for int() with base 10: '3.142'
```

# Floating Point Numbers

-   Represent a finite approximation to the real numbers.

-   Type: `float`.

-   (On most platforms) Python uses IEEE-754 double precision floating
    point numbers which provide 53 bits of precision.

-   `sys.float_info` contains details about max, min, epsilon etcetera.

-   Literals can be in ordinary notation or in exponential notation:
    -   Ordinary: `3.142`
    -   Exponential: `314.2e-2`

-   Ordinary notation requires a point `.`, but digits following the
    point are optional.

-   Exponential notation does not require a point unless you have a
    fractional component.

-   `float(x)` will try to convert x to a floating point number, x can be
    another numeric type (including booleans) or a string.

-   Numeric operators will automatically convert integer arguments to
    floating point in mixed-type expressions.

-   In Python 3 the division operator `/` computes a floating point result
    for integer arguments. However, in Python 2 it computes an integer
    result for integer arguments.

**Interactive examples:**

```
>>> type(3.142)
<class 'float'>
>>> type(12)
<class 'int'>
>>> 3.142 + 12
15.142
>>> 3.142 == 314.2e-2
True
>>> 3. == 3.0
True
>>> 1/0
Traceback (most recent call last):
File "<stdin>", line 1, in <module>
ZeroDivisionError: division by zero
>>> # Integer divided by integer yields a float in Python 3
>>> 10 / 3
3.3333333333333335
>>> float("123")
123.0
>>> float("3.142")
3.142
```

# Complex Numbers

-   Represent a finite approximation to the complex numbers.

-   Type: `complex`

-   A pair of floating point numbers: real +/- imaginary.

-   The real part is optional (defaults to 0). The imaginary part is
    followed immediately by the character `j`.

**Interactive Examples:**

```
>>> 5j + 3j
8j
>>> 2-5j
(2-5j)
>>> 2-5j + 3j
(2-2j)
```

# Numeric Operators

Name | Operation    | Precedence | Associativity | Notes
-----|--------------|------------| --------------| -------------------------------------------------------------------------------------------------------------------------------
`+`  | add          | low        | left          | Can also be used to concatenate strings together.
`*`  | multiply     | medium     | left          | 
`-`  | subtract     | low        | left          | 
`/`  | divide       | medium     | left          | In Python 3 the result is always a floating point number. In Python 2 the result is an integer if both operands are integers.
`//` | floor-divide | medium     | left          | divide then floor, result is an integer
`**` | exponent     | high       | right         | 
`%`  | modulus      | medium     | left          | remainder after division

**Interactive Examples (Python 3):**

```
>>> 3 + 4 * 5
23
>>> (3 + 4) * 5
35
>>> 10 / 3
3.3333333333333335
>>> 10 // 3
3
>>> 10 % 3
1
>>> 2 ** 3 ** 4
2417851639229258349412352
>>> (2 ** 3) ** 4
4096
```

# Strings

-   Represent text

-   Type: `str`

-   In Python 3, the str type contains Unicode characters.

-   In Python 2, the str type contains ASCII characters (sometimes
    called byte strings). Python 2 has a separate type for unicode
    strings, the type is called unicode; literals of this type are
    prefixed by the letter `u`.

-   String literals must be quoted. There are 3 quoting styles:

    -   single quote characters: `'hello'`

    -   double quote characters: `"hello"`

    -   triple quote characters: `'''hello'''` (three single quotes in a
        row) or `"""hello"""` (three double quote characters in a row)

-   The single quote and double quote versions of strings have the same
    value. The purpose of the different quotation styles is to make it
    convenient to have literal quotation marks inside strings
    (avoiding the need to escape the quote character). For example:

```
>>> "This inverted comma won't be a problem inside quotation marks"
"This inverted comma won't be a problem inside quotation marks"
>>> 'this "quote" will work'
'this "quote" will work'
>>> 'this isn't going to work though'
File "<stdin>", line 1
'this isn't going to work though'
^ SyntaxError: invalid syntax
```

-   Triple quoted strings can be written on multiple lines. The line
    breaks will be preserved within the string. Useful for docstrings
    (see section on functions).

-   The usual set of escape characters are supported:

    -   `\n` newline

    -   `\t` tab

    -   `\\` backslash

    -   `\'` single quote

    -   `\"` double quote

    -   and many more

-   Python does not have a separate type for representing individual
    characters. Instead you use strings of length one.

-   Strings are iterable. If you iterate over a string (using a for
    loop) you process it one character at a time from left to right.

-   Strings can be indexed to obtain individual characters, e.g. `s[5]`

-   Indices are zero-based (but you may also use negative indices to
    access items with respect to the right end of the string).

-   Strings are immutable: you cannot modify a string once it has been
    created.

**Interactive Examples (Python 3):**

```
>>> type("hello")
<class 'str'>
>>> "hello" == 'hello'
True
>>> '''This string
... is on
... multiple
... lines'''
'This string\\nis on\\nmultiple\\nlines'
>>> "bonjour".upper()
'BONJOUR'
>>> len("bonjour")
7
>>> "bonjour".startswith("b")
True
>>> "cat,sat,flat".split(",")
['cat', 'sat', 'flat']
>>> # Print the first 5 Chinese unicode characters
>>> print('\u4E00\u4E01\u4E02\u4E03\u4E04')
一丁丂七丄
>>> x = "floyd"
>>> x[0]
'f'
>>> "hello" + " " + "world"
'hello world'
```

**Example program:**

```python
# Prompt the user to input a string:
input = raw_input("Enter string: ")

# Count the number of vowels in the input string
vowels = 'aeiou'
count = 0

for char in input:
    if char in vowels:
        count += 1

# Print the count to the standard output
print(count)
```

Example usage of the above program from the operating system command
prompt, assuming the program is saved in a file called `vowels.py`:

```no-highlight
python vowels.py
Enter string: abracadabra
5
```

# Lists

-   Represent mutable ordered sequences of values.

-   Type: `list`

-   List literals are written in between square brackets, e.g. `[1, 2, 3]`

-   List elements can be objects of any type (including other lists).

-   Like strings, lists can be indexed like so: `x[3]`

-   Indices are zero-based (but you may also use negative indices to
    access items with respect to the right end of the list).

-   Lists are mutable. You can update items, delete items and add new
    items.

-   Indexing into a list is a constant time (amortised) operation.

**Interactive Examples:**

```
>>> type([1, 2, 3])
<class 'list'>
>>> x = []
>>> len(x)
0
>>> x.append("hello")
>>> x
['hello']
>>> len(x)
1
>>> x[0]
'hello'
>>> x.insert(0, True)
>>> x
[True, 'hello']
>>> del x[1]
>>> x
[True]
>>> x += [42, "Newton", 3.142]
>>> x
[True, 42, 'Newton', 3.142]
```

# Dictionaries

-   Represent finite mappings from keys to values.

-   Are implemented as [hash tables](https://en.wikipedia.org/wiki/Hash_table). The key objects must be hashable
    (which rules out mutable objects, such as lists).

-   Type: `dict`

-   Dictionary literals are written inside curly brackets, with
    key-value pairs separated by colons: e.g. `{12: "XII", 6: "VI"}`

-   Dictionaries can be indexed by keys. If the key exists in the
    dictionary its corresponding value is returned, otherwise a
    `KeyError` exception is raised.

-   The cost of indexing a dictionary is proportional to the time taken
    to hash the key. For many keys this can be considered constant
    time. For variable sized objects, such as strings, this can be
    considered to be proportional to the size of the object.

-   Iterating over a dictionary yields one key at a time. All keys in
    the dictionary are visited exactly once. The order in which the
    keys are visited is arbitrary.

-   You may test if an object is a key of a dictionary using the in
    operator.

**Interactive Examples:**

```
>>> type({12: "XII", 6: "VI"})
<class 'dict'>
>>> friends = {}
>>> friends['Fred'] = ['Barney', 'Dino']
>>> friends
{'Fred': ['Barney', 'Dino']}
>>> friends['Fred']
['Barney', 'Dino']
>>> friends['Barney']
Traceback (most recent call last):
File "\<stdin\>", line 1, in \<module\>
KeyError: 'Barney'
>>> friends['Wilma'] = ['Betty']
>>> friends
{'Fred': ['Barney', 'Dino'], 'Wilma': ['Betty']}
>>> friends.keys()
dict_keys(['Fred', 'Wilma'])
>>> friends.values()
dict_values([['Barney', 'Dino'], ['Betty']])
>>> 'Dino' in friends
False
```

**Example program:**

```python
# Compute and print a histogram of a sequence of integers entered
# on standard input, one number per line

import sys

histogram = {}

# Iterate over each line in the standard input
for line in sys.stdin:
    # Parse the next input as an integer
    next_integer = int(line)
    # Update the histogram accordingly
    if next_integer in histogram:
        # We've seen this integer before
        histogram[next_integer] += 1
    else:
        # First occurrence of this integer in the input
        histogram[next_integer] = 1

# Print each key: value pair in the histogram in ascending
# sorted order of keys
for key in sorted(histogram):
    print("{} {}".format(key, histogram[key]))s
```

Example usage of the above program from the operating system command
prompt, assuming the program is saved in a file called `histo.py`:

```no-highlight
python histo.py
```

User types in a sequence of integers to the program, one per line, and
presses control-d to terminate the input:

```no-highlight
3
43
12
19
3
12
12
43
```

Program prints its output:

```no-highlight
3 2
12 3
19 1
43 2
```

# Tuples

-   Represent *immutable* ordered sequences of values.

-   Very much like lists except they cannot be modified once created.

-   Type: `tuple`

-   Literals are written in between parentheses: `(1, 2, 3)`

-   The can be used as keys in dictionaries (unlike lists).

# Loops

## While loops

-   Iterate until condition is `False`

-   Syntax:

```python
while expression:
    statement_block
```

-   The value of the boolean expression is tested. If it evaluates to
    `True` then the statement block is executed once, before repeating
    the loop. If it evaluates to False then the program continues
    execution immediately after the loop.

**Example:**

```python
def factorial(n):
    result = 1
    while n > 0:
        result *= n
        n -= 1
    return result
```

For loops
---------

-   Iterate over each item in a collection (e.g. list, string, tuple,
    dictionary, file).

-   Syntax:

```python
for variable in expression:
    statement_block
```

-   Each item from the iterator expression is selected and assigned to
    the variable, then the statement block is executed. The loop ends
    when every item in the iterator has been visited.

-   The order of items visited in the iterator depends on the type of
    the iterator. Lists, strings and tuples proceed in a left-to-right
    fashion. Files proceed one line at a time. Dictionaries proceed in
    an arbitrary order.

-   The `range()` function is useful for generating iterators of numbers
    within a range. Note that the lower bound is inclusive and the
    upper bound is exclusive.

**Example:**

```
def factorial(n):
    result = 1
    for item in range(n + 1):
        result *= item
    return result
```

## Break and continue

-   Both types of loops support the `break` and `continue` keywords.

-   `break` terminates the loop immediately.

-   `continue` jumps immediately back to the start of the loop.

-   They can sometimes simplify the conditional logic of a loop, but
    should be used sparingly.

# Functions

-   Allow you to define reusable abstractions. Sometimes called
    *procedures*.

-   Are generally defined at the top level of a module, and can also be
    nested.

-   Type: `function`

-   Named functions are bound to a variable name and may have complex
    bodies.

-   Anonymous functions are used in-line, and may only have expression
    bodies.

-   Named function syntax:

```
def variable(parameter_list):
    statement_block
```

-   Anonymous function syntax:

```
lambda parameter_list: expression
```

**Example:**

```python
def is_leap_year(year):
    if year % 4 == 0 and year % 100 != 0:
        return True
    else:
        return year % 400 == 0

for year in range(2000, 2100 + 1):
    result = is_leap_year(year)
    print("{} {}".format(year, result))
```

**Anonymous function example:**

```python
>>> squared = lambda x: x ** 2
>>> squared(2)
4
>>> list(map(lambda x: x + 1, [1, 2, 3]))
[2, 3, 4]
```

# Input and output

-   The `print` function is useful for displaying text (and other values
    converted to text).

-   In Python 2 `print` was a special keyword. In Python 3 it is a
    function defined in the builtins.

-   Fancy string formatting can be done with the format method on
    strings. Older Python code uses the string interpolation operator
    for the same task `%`, but its use is now discouraged.

-   Files must be opened before than can be manipulated. A file can be
    opened in different modes: read `"r"`, write `"w"`, read-write `"r+"`,
    and append `"a"`. Opening a new file in write or append modes
    creates a new file. Opening an existing file in write mode
    overwrites its contents from the start. Opening an existing file
    in append mode adds new content at the end of the old content.

-   When you are finished processing a file you should close it as soon
    as possible. Closing a file releases limited operating system
    resources, and ensures that any pending buffered writes a flushed
    to the storage system.

-   Certain file types have libraries for convenient processing. One
    example is the CSV (comma separated values) library for processing
    tabular data. It is very handy for working with spreadsheets.

-   The command line arguments of a Python program are contained in a
    list called `sys.argv` (it is a variable exported from the `sys`
    module). For complex program you should consider using a command
    line argument parsing library such as [argparse](https://docs.python.org/3/library/argparse.html).

**Example program:**

```
# Count the number of words and lines in a file

import sys

# Get the input file name from the command line arguments
filename = sys.argv[1]

# Open the file
file = open(filename)

# Count the number of lines in the file
num_lines = 0
num_words = 0

for line in file:
    num_lines += 1
    num_words += len(line.split())

file.close()

print("Number of lines and words in {}: {} {}" \
        .format(filename, num_lines, num_words))
```

# Advanced Topics

## Classes

-   Classes allow you to define your own types.

-   Class definitions may define methods for the type.

-   A class may inherit, and possibly override, some functionality from
    a superclass.

-   Syntax:

```
class variable(superclass_list):
    body
```

-   The name of the class is given by the variable in the definition.
    The superclass list defines the superclasses of the new class
    (very often the base type object is used). The body of the class
    typically defines one or more methods.

-   Instances of classes are created by calling the class name as if it
    were a function.

-   If defined, the special method called `__init__` is used to
    initialise a newly created instance of a class.

-   The first parameter to each method is the object upon which the
    method was called. The convention is to use the variable called
    self, however any variable name will do. Many object oriented
    languages make this variable an implicit parameter called this.

**Example:**

```python
class Vector(object):
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def magnitude(self):
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalise(self):
        magnitude = self.magnitude()
        if magnitude == 0:
            # Somehow we have a degenerate vector.
            return self
        else:
            return self / self.magnitude()

    def angle(self, other):
        dp = self.dot_product(other)
        return acos(dp / self.magnitude() * other.magnitude())

    def dot_product(self, other):
        return self.x * other.x + self.y * other.y + self.z * other.z
```

# Exceptions

-   Exceptions allow Python programs to handle erroneous program
    conditions.

-   An exception is raised (or thrown) at the point of the error and
    handled (or caught) at some other place in the program.

-   Exception handlers have the syntax:

```python
try:
    statement_block
except exception_type as variable:
    statement_block
...
```

-   The statement block after `try` is executed. If no exceptions are
    raised in that block the program continues immediately after the
    exception handler. If an exception is raised in the block then
    program control jumps to the innermost closing `except` clause.
    Except clauses may optionally specify the set of exception types
    that they can handle. If the raised exception is an instance of
    the handled type then the body of the except clause is executed,
    otherwise the next except clause (if any) is tried.

-   If no matching exception handler is found then the program will
    terminate with an unhandled exception error. Python will normally
    print a stack trace at this point for error diagnosis.

-   You may raise your own exceptions using the `raise` keyword.

**Example:**

```python
# alternative version of the histogram code from the section on
# dictionaries
for line in sys.stdin:
    next_integer = int(line)
    try:
        histogram[next_integer] += 1
    except KeyError:
        histogram[next_integer] = 1
```

Modules
-------

-   A module is a file which contains Python code.

-   Any Python file you create is automatically a module.

-   It is considered good programming style to decompose complex
    programs into multiple modules. Each module should collect
    together code with similar purpose.

-   Variables defined at the top level of a module (such as global
    variables, functions and classes) can be imported into other
    modules.

-   Python comes with many standard modules.

-   The `import` keyword is used to import an entire module.

-   You may import a subset of things from a module using the `from ... import ...` syntax.

-   You may import a module with a new name using the `from ... import ... as ...` or `import ... as ...`

-   When a module is first imported in a program, all of its top-level
    statements are executed from top to bottom. Subsequent imports use
    a cached version of its definitions, its statements are not
    re-executed.

-   A special module called `builtins` is imported into every other module
    by default, and it is automatically imported at the interactive
    prompt in the interpreter.

**Interactive Example:**

```python
>>> import math
>>> math.sqrt(100)
10.0
>>> sqrt(100)
Traceback (most recent call last):
File "<stdin>", line 1, in <module>
NameError: name 'sqrt' is not defined
>>> from math import sqrt
>>> sqrt(100)
10.0
>>> import math as m
>>> m.sqrt(100)
10.0
```

# Packages

-   A package is a collection of modules in a hierarchy.

-   Packages are the common way to structure Python libraries.

-   The [Python Package Index (PyPI)](https://pypi.python.org/pypi) is a big collection of open
    source packages contributed by the Python community (PyPI contains more than 64 thousand packages at the time of
    writing).

-   Package installation tools such as [pip](https://en.wikipedia.org/wiki/Pip_(package_manager)), make it easy to install
    packages onto your computer.

-   If you want to make your own Python code easy for others to install
    and use then you should consider making it a package. You can even
    upload it to PyPI.

-   Many people use [virtualenv](http://docs.python-guide.org/en/latest/dev/virtualenvs/) to install packages into a local
    "sandboxed" Python environment. This avoids conflicts with the
    central Python package database on your computer, and allows
    multiple different versions of packages to be used.
