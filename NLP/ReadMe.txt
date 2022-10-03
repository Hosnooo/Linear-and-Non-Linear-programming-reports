Notes: The methods are implemented to support any dimension./ f_as_lamda.m and gradf_as_lamda.m
       are just used to get any function in terms of lamda (1-D).

######To see the report results
1. Put 'Powel.txt','Rosenbrock.txt' and 'runme.m' in the same directory.
2.Open runme.m
3.Run

######To generate new functions
1. Open Generate.m
2. Type down the functions you want to minimize.
3. Make sure you pass it to all methods
4. Make sure that the table function recieves consistent entries;
   the function dimension changes the number of variables (x1,x2,x3,..etc),
   so remember to add them in the vector'vars' and in the table columns.

######To run any stand alone function
1. Just call the function you want in a new file.
