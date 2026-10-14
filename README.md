# MATLAB reference implementation

This directory contains the standalone MATLAB code accompanying the revised
manuscript described below.

## Manuscript

**Title:** Subband-Selective Sign Subband Adaptive Filtering Algorithms With
Band-Dependent Variable Step Sizes

**Authors:** Bum Yong Park, Jinwoo Yoo, Won Il Lee, and JaeWook Shin

**Journal:** Signal Processing

**Manuscript ID:** SIGPRO-D-26-02018

The code runs the abrupt-system-change experiment for the following six
algorithms presented in the manuscript:

- M-SSAF (SSAF normalization)
- M-SSAF (IWF-SSAF normalization)
- SS-M-SSAF (SSAF normalization)
- SS-M-SSAF (IWF-SSAF normalization)
- VSS-SS-M-SSAF (SSAF normalization)
- VSS-SS-M-SSAF (IWF-SSAF normalization)

`SSAFinit.m` designs the prototype lowpass filter and constructs the analysis
filters directly from the cosine-modulation equation in the manuscript. No
external book code is required.

## Implementation reference

The subband adaptive-filter implementation structure and the cosine-modulated
filter-bank construction were developed with reference to the concepts and
companion MATLAB code associated with:

K. A. Lee, W. S. Gan, and S. M. Kuo, *Subband Adaptive Filtering: Theory and
Implementation*, Wiley, 2009.

The companion MATLAB source files from the book are not included in this
release. `SSAFinit.m` is an independent implementation that generates the
analysis filters from the equation stated in the submitted manuscript; the
synthesis filters are obtained by time-reversing the analysis filters.

## Run

Add this directory to the MATLAB path and run:

```matlab
results = run_demo;
```

The demo uses a colored input, Bernoulli-Gaussian impulsive noise, and an
abrupt system change from `w` to `-w`. It plots the one-run NMSD curves of all
six algorithms. The VSS result structures also contain `mu`, `ctrl`, and
`reset_event`.

MATLAB and Signal Processing Toolbox are required.
