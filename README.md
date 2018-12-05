# euler1D
This repository contains python files to execute an implementation of the sod shock tube problem.  

## Files
The main files used in this implementation of euler1D is
#### euler1D.py
Executing this file as is will execute the sample problem setup with the sample
domain setup found in the euler1D.txt file.  

## Supporting files and modules
nodeFileGenerator.py - this file can be used to create a new domain text file.  
dataPlot.py - this is used to analyze all the results in a directory.  
              update this file appropriately in the file __main__ to change the
              directory and other properties.  
eulerExact.py - this is the analytical solution of Euler's equations applied
                to the sodshock tube problem. - Currently not properly functioning  
sodshock.py - this is the analytical solution of the sodshock tube problem  with
              constant states - Not validated
### shocktubecalc module
is a module from PyPi. This is a validated analytical solution to the  
sodshock tube problem. Used to compare to the numerical solution.  
Note that the functions from sodshock.py are used to save this  
in a consistent form.

### fluid_domain module
domain.py - use this file to create a domain from a given text file.  
decompose.py - use this file to decompose the domain if desired.  
node.py - this file is used by the domain class. Essentially, domain is a
          multidimensional tuple of node objects
