
### Description
This package will read in the scan data, fit and determine the optimal settings for the HWP and QWP so that we have the highest degree of circular polarization in the middle of cavity.

The rootlogon.C file will automatically compile the three classes needed. It should work with rootv6.12 or less.

The plotMacros contain root macros to plot the output.

### Example command for running the scan:
```
Minimize("EFscan_10deg",1,2,"","data/runpoints.txt","/home/compton/gaskelld/entrancefunction/TransferFunction/data/output","ent_scan_09-25-2015_norm_10deg.dat")
```
