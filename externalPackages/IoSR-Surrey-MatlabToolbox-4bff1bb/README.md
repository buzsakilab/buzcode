# IoSR Matlab Toolbox

A general purpose Matlab toolbox containing functions and classes for: auditory modelling, signal processing, sound source separation, statistics, plotting, etc. See [Contents.m](https://github.com/IoSR-Surrey/Toolbox/blob/master/+iosr/Contents.m) for a full list of functions/classes.

## Installation

Basic installation only requires you to add the install directory to the Matlab search path.

If you wish to perform certain audio / signal processing tasks (especially spatialisation), please navigate to the install directory and type

```
iosr.install
```

This will automatically download the toolbox's dependencies for these tasks, and add the necessary paths to your search path.

## Usage

Use these functions as:

```
iosr.<folderName>.<functionName>(<args>)
```

(Ignoring the '+' in the folder name.) Alternatively, use the `import` directive to add one or more namespaces, e.g.:

```
import iosr.auditory
import iosr.*
```

If using `import`, note that some function names may conflict with built-in Matlab function names (e.g. `quantile`). One method of resolving the conflict and shortening the function call is to create a handle to any functions with conflicting names, e.g.

```
qntl = @iosr.statistics.quantile;
```

Type

```
help iosr
```

for more information.