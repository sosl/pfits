# This is a clone of pfits by George Hobbs and Lawrence Toomey.

For the original, see https://bitbucket.csiro.au/projects/PSRSOFT/repos/pfits/

It was cloned in mid-December 2018 and may diverge henceforth

The code is licensed under CSIRO Open Source Software Licence v1.0, see [OSS_LICENCE.txt](./OSS_LICENCE.txt)

## original README:

The pfits software package is developed to process psrfits data files (both search and fold-mode files).

Software to write:

1. pfits_describe:

summarises the header information in a PSRFITS file.  The output can be written to the screen or to a file.

2. pfits_fv:

code that reproduces some of the "fv" (fits viewer) functionality. This code allows the user to obtain detailed information about the psrfits file.

3. pfits_plot:

This code first determines whether the data is search mode or fold mode.

fold mode:

if 1 channel then it plots the profile (with polarisation information if present)
if multiple channels then it plots of colour scale image for I, Q, U and V (unless just I is requested)

search mode:

the user can select a time range, or nsamp range, or subint range. If > N points (where N is typically 4096) then the data gets averaged and the minimum, mean and maximum value in the average is recorded. The following options are available:

- sum in frequency
- dedisperse in frequency (prior to summing if requested)
- form stokes I
- plot 4-pol data underneath each other
- plot frequency-time plots

4. pfits_dedisperse

read in a psrfits search mode file, dedisperses, frequency-sums (if requested) and a new psrfits file is output.  Have the option to output multiple dispersed data files if requested.

5. pfits_fold

folds a psrfits search mode file with a given pulse period or using a tempo2 predictor.  Outputs a new psrfits fold-mode file

6. pfits_singlepulse

extracts single pulses from a psrfits search mode file.  Can output psrfits fold-mode files for each pulse or calculate various parameters (e.g., flux density) and output those instead.

7. pfits_calibrate


