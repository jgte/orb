# orb

Generic timeseries/satellite geodesy/data processing Matlab software.

**Disclaimer:** much of this software is under development and may not run as expected. There are also many instances of outdated software. Please report any strange behaviour and I'll do my best to get it up and running.

## Directory structure

Each checked-out copy of this repository is intended to be its own self-contained *project*. That means that is should contain all data, metadata, plots and code. It is up to the user to refrain from pushing scripts that are specific to a certain application that isn't really interesting to anyone else (and that may be the case for the scripts listed in the [`Applications`](#Applications) section.)

The following directories make up the repository:

* aux : contains auxiliary data, such as C20 time series, land makes and geopotential coefficients
* data : this is where the data are stored by the `datastorage` class
* metadata : this is where the metadata are stored by the `datastorage` class
* packages : contains auxiliary software (currently `+yaml`)
* plots :  : this is where the plots are stored by the `datastorage` class
* version_patching : contains functions that were not available in old Matlab versions (in case those versions are the only possibility)

## `datastorage`

The following classes make up the *datastorage* infrastructure (in addition to the data, metdata and plots directories):

* datanames : resolves the name of a data product
* dataproduct : handles a data product metadata and relevant directories
* datastorage : main script for loading, processing and saving data products

This infrastructure is intended to make it easier to handle data processing pipelines. The idea is that it is all controlled by the metadata and by splitting the processing into steps, each one corresponding to a certain data product.

The use of this infrastructure has a somewhat steep learning curve and is prone to breaking unless everything lines up perfectly. Fortunately, it is completely de-coupled from the remaining classes. For a good example on how to use it, refer to the `gswarm` class.

## Dummy classes

Each mat-file is a class but many are *dummy* classes, which are only used to group together a bunch of routines that have similar application. Here's an template of those *dummy* classes:

```
classdef dummy
  methods(Static)
    function [out11,...,outi1]=f1(in11,...,ink1)
    [...]
    end
    [...]
    function [out1j,...,outij]=fj(in1j,...,inkj)
    [...]
    end
  end
end
```

This makes it possible to call each of the `f1` ... `fj` members simply as:

```
dummy.fJ(in1J,...,inkJ)
```

The dummy classes may have class constants but do not have non-static members.

The following m-files contain dummy classes:

* cb.m : handles colorbars
* cells.m : utilities relevant to cell arrays
* cluster.m : handles staging and un-staging of data
* file.m : file operations
* num.m : numerical algorithms
* pardecomp.m : parameter decomposition
* plotting.m : wrapper and default-setter for plotting
* str.m : string operations
* structs.m : utilities relevant to structures
* time.m : time-related functions
* url.m : utilities relevant to handling URLs


## Classes

Each of the following classes handle specific physical data types, which are shown according to their class hierarchy:

* attitude : process quaternion or angular data
* orbit : process orbit data
* simpledata : contains primitive methods for data processing, models data as a numeric vector with a common value for the free variable and for the *mask* (valid/invalid data point).
  * simpletimeseries : adds time-related operations
    * simplefreqseries : frequency representation of time series
      * segmentedfreqseries : segmented representation of spectra
    * simplegrid : represents gridded data
    * gravity : handling of Spherical Harmonic coefficients

The `varargs` class is different from all others in the sense that it handles variable-length arguments (TODO: add a better description).


## Interface classes

The following classes make is possible to retrieve data from NRTDM (please ignore if you are unfamiliar with this software package):

* nrtdm
* nrtdm_metadata
* nrtdm_product


## Scripts and functions

* startup.m : handles correct linking of project directories
* test.m : calls all tests in all classes


## Applications

The following classes are examples of using the classes and scripts above.

* csr.m : routines to handle some of CSR's data, specifically those related to accelerometer parameters (highly outdated, possibly to be removed in the future)
* gswarm.m : used for the processing of the Swarm gravity field models


