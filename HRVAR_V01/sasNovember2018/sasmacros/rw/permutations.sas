*_________________________________________________________________________________________

	%permutations(dataset,n) generates all the permutations of order n.
	For example, %permutations(test,3) generates a dataset ('test') with three variables (i1, i2, i3)
	that give all of the permutations:
	i1    i2    i3
	1     2     3
 	1     3     2
 	2     1     3
 	2     3     1
 	3     1     2
 	3     2     1
 	%permutations(dataset, n, AllPerms=0) causes only the cyclic permutations to be generated.
 	For example, after %permutations(test,5,AllPerms=0), the dataset 'test' is:
   i1    i2    i3    i4    i5
    1     2     3     4     5
    2     3     4     5     1
    3     4     5     1     2
    4     5     1     2     3
    5     1     2     3     4
__________________________________________________________________________________________;

%macro permutations(dataset,n,AllPerms=1);
%if &AllPerms=1 %then %do;
	data _null_;
		m = fact(&n);
		call symput('m',m);
		run;
	proc plan ordered;
		factors t=&m ordered f=&n perm / noprint;
		output out=perm;
		run;
%end;
%else %do;
	proc plan ordered;
		factors t=&n ordered f=&n cyclic / noprint;
		output out=perm;
		run;
%end;
proc transpose data=perm prefix=i out=&dataset (drop=_name_ t);
	by t;
	run;
%mend Permutations;

/*
*_________________________________________________________________________________________

	Test code
__________________________________________________________________________________________;

options nocenter mprint;
%permutations(Perm0, 3);
proc print data=Perm0;
	title "All permunations";
	run;
%permutations(Perm1, 5, AllPerms=0);
proc print data=Perm1;
	title "Cyclic permunations";
	run;
options nomprint;
*/