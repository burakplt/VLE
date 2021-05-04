# VLE
Vapor-Liquid Equilibrium calculator with ChemSep chemical database.
***Read Wiki page for further information.***
I created this module originally as a part of Django Web application project. I will modify them in time and update this repository.

Before running any command, you have to create chemical object pickles to increase run-time performance. As long as the project folder 
structure or "Chemical" class definition remains unchanged, the pickle objects are created only once. 

Simply run in "chemsep_operation" file the command "create_chemical", and it will create the required chemical pickle objects in /chemical 
folder. Then you can use all the functionality really fast without parsing XML database every time. 

Django Web application related to this topic is under development, but it can be reached on <s>www.chemicaleasy.com </s> now. 


...Will be updated
