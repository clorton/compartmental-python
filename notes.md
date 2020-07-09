# Assessment Notes

## GillesPy2

- should work "out-of-the box" on Windows: almost did, required selecting the NumPy solver rather than auto which seemed to choose the C/C++ solver
- has an interface for adding species, parameters, and reactions similar to what I envision for PyCMS (alternative to parsing from EMODL text)
- documentation is "developer" level - lacking in conceptual information, e.g. mass action w/rate vs. custom propensity function
- actively being developed: [latest commits](https://github.com/GillesPy2/GillesPy2/commits/master)

## StochPy

- last update seems to have been in August of 2015
- documentation says to prefer Python 2.7
- uses PySCeS model description language (can also use some SBML)
- events appear to be broken - they assume a fixed value being assigned rather than evaluating their assignment expressions
- events also do not seem to be able to update a parameter, only species