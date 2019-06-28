# euler1D
This repository contains python files to solve the 1D euler equations numerically  
with RK2 in time and a 5 point finite volume method in space. Currently, it is  
set up to solve a staple problem (the sod shocktube). The analytical solution of  
this problem is included in two forms. One from sodShock.py and another in the shocktubecalc module. 

## Files
The main files used in this implementation of euler1D is

#### euler1D.py
Executing this file as is will execute the sample problem setup with the sample
domain setup found in the euler1D.txt file. This solution has been validated  
with the analytical solution.  

## Supporting files and modules
nodeFileGenerator.py - this file can be used to create a new domain text file.  
dataPlot.py - this is used to analyze all the results in a directory.  
              update this file appropriately in the file __main__ to change the
              directory and other properties.  
sodShock.py - this is the analytical solution of the sodshock tube problem  
              this solution is validated with the PyPi sodshock calculation

### shocktubecalc module
is a module from PyPi. This is an analytical solution to the  
sodshock tube problem used to validate the numerical solution.  
Note that the functions from sodshock.py are used to save this  
in a consistent form.

### fluid_domain module
domain.py - use this file to create a domain from a given text file.  
decompose.py - use this file to decompose the domain if desired.  
node.py - this file is used by the domain class. Essentially, domain is a
          multidimensional tuple of node objects
### misc
This folder contains former code and other files that are no longer in use
