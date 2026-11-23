*_________________________________________________________________________________________

	RandomWalkAnalysis
	Version 1.2
	
	Joel Hasbrouck
	Department of Finance, 9-88 KMEC
	Stern School of Business
	New York University
	44 West 4th St.
	New York, NY 10012

	email: jhasbrou@stern.nyu.edu

	Examines all rotations (permutations) of the variables to determine the variance
	of the random walk component and the minimum and maximum information shares.
	
	Parameters

	dsResults	Output dataset that will contain results
	dsCov		Input dataset that contains the covariance matrix of disturbances
	dsCoeff		Input dataset that contains the coefficient matrix
	Keyword Parameters:
	PrintLevel	0	No printed output
				1	Summary printed output (default)
				2	Detail on all permutations
	AllRows		1	(default) Compute the random walk analysis for all rows of
					the coefficient matrix.
				0	Use only the first row of the coefficient matrix

				For the usual price discovery analysis, the rows are identical.
				In this case, set AllRows=0.  (If you do this, you should make sure
				that the rows are indeed identical.  See the RowDifference macro.) 

	A test driver program follows the macro definition.
__________________________________________________________________________________________;

%macro RandomWalkAnalysis(dsResults, dsCov, dsCoeff, PrintLevel=1, AllRows=1, AllPerms=1);
%*	Get number of variables;
data _null_;
	if 0 then set &dsCov nobs=countxx;
	call symput("n",compress(put(countxx,12.)));
	stop;
	run;
%Permutations(dsPerm,&n,AllPerms=&AllPerms);	%*	Build dataset of permutations;
proc iml;

	start main;
	if &PrintLevel>=1 then do;
		print "Random walk analysis";
	end;
	use &dsCoeff;
	read all into bFull [colname=VarNamesb];
	mattrib bFull rowname=VarNamesb label="Coefficient matrix";
	if &PrintLevel>=1 then do;
		print bFull;
	end;
	n = ncol(bFull);
	nRowb = nrow(bFull);
	if n^=nRowb then print "Warning: Coefficient matrix not square.";
	use &dsCov;
	read all into cov [colname=VarNames];
	mattrib cov rowname=VarNames label="Disturbance covariance matrix";
	if &PrintLevel>=1 then do;
		print cov;
	end;
	if n^=nrow(cov) then do;
		print "Covariance matrix not square or not conformable with coefficients.";
		print cov bFull;
		abort;
	end;
	CoeffRowVar="        ";
	Type = "               ";
	Types={"VarRandWalk" "sdRandWalk" "InfoShareMin" "InfoShareMax"};
	create &dsResults var {CoeffRowVar Type};
	create &dsResults.Values var ({"iRow" "VarRandWalk" "sdRandWalk"} || VarNamesb);
	do iRow=1 to nRowb;
		b = bFull[iRow,];
		mattrib b colname=VarNamesb label="Coefficients";
		if &PrintLevel>=1 then do;
			print "Random walk analysis for row " (iRow[1])[format=2.0]
				"Variable:" (VarNamesb[iRow]);
			print b;
		end;
		call RandWalkAnalysis(VarRandWalk, InfoShareMin, InfoShareMax, b, cov);
		sdRandWalk = sqrt(VarRandWalk);
		mattrib sdRandWalk label="sdRandWalk";
		mattrib VarRandWalk label="VarRandWalk";
		mattrib InfoShareMin colname=VarNames label="InfoShareMin";
		mattrib InfoShareMax colname=VarNames label="Info Share Max";
		if &PrintLevel>=1 then do;
			print sdRandWalk;
			print VarRandWalk;
			print InfoShareMin;
			print InfoShareMax;
		end;
		setout &dsResults.Values;
		a = j(3,1,iRow) || 
			((VarRandWalk || sdRandWalk || j(1,n,.)) // (j(1,2,.) || InfoShareMin) // (j(1,2,.) || InfoShareMax));
		append from a;
		*	write row variable names;
		setout &dsResults;
		do i=1 to nrow(a);
			vn = VarNamesb[iRow] || Types[i];
			append from vn;
		end;
		if ^&AllRows then goto FinishUp;
	end;
	FinishUp: 
	finish main;

	start RandWalkAnalysis(VarRandWalk, InfoShareMin, InfoShareMax, b, cov);
	*	Inputs are b (the vector of coefficients) and the covariance matrix;
	*	Check inputs;
	n = ncol(b);
	if n^=nrow(cov) | n^=ncol(cov) then do;
		print b cov "b and cov don''t conform.";
		abort;
	end;

	*	Loop over all permutations to determine min and max;
	use dsPerm;
	read all into Perms;
	nPerm = nrow(Perms);
	InfoShare = j(nPerm,n,0);
	VarRandWalk = j(nPerm,1,0);
	do i=1 to nPerm;
		p = Perms[i,];
		call varcal(VarRandWalk2,InfoShare2,b,cov,p);
		InfoShare[i,] = InfoShare2;
		VarRandWalk[i] = VarRandWalk2;
		if &PrintLevel>=2 then do;
			print "RandWalkAnalysis iteration " i;
			print "Permutation " p;
			print "Information Share " InfoShare2;
		end;
	end;
	InfoShareMin = j(1,n,0);
	InfoShareMax = j(1,n,0);
	do i=1 to n;
		InfoShareMin[i] = min(InfoShare[,i]);
		InfoShareMax[i] = max(InfoShare[,i]);
	end;
	*	Verify that all permutations give the same random walk variance;
	d = log(max(VarRandWalk)/min(VarRandWalk));
	if d>1e-10 then do;
		print "Possible inconsistency in random walk analysis:" d;
	end;
	VarRandWalk = VarRandWalk[1];
	finish RandWalkAnalysis;

	start varcal(VarRandWalk,InfoShare,b,cov,perm);
	*	Performs the variance and info share calculations for a single permutation;
	bP = PermuteMatrix(b,perm,1);
	covP = PermuteMatrix(cov,perm,1);
	VarRandWalkContrib = (bP*t(root(covP)))##2;
	VarRandWalk = sum(VarRandWalkContrib);
	InfoShare = VarRandWalkContrib/VarRandWalk;
	InfoShare = PermuteMatrix(InfoShare,perm,-1);
	finish varcal;

	start PermuteMatrix(a,p,direction);
	*	Reorders the variables in preparation for Cholesky decomposition;
	*	direction=+1 for the forward permutation and -1 to reverse the permutation;
	aP = a;
	if nrow(a)=1 | ncol(a)=1 then do i=1 to max(ncol(a),nrow(a));
		k = p[i];
		if direction=1 then aP[i] = a[k];
		else aP[k] = a[i];
	end;
	else do i=1 to nrow(a);
		do j=1 to ncol(a);
			k = p[i];
			l = p[j];
			if direction=1 then aP[i,j] = a[k,l];
			else aP[k,l] = a[i,j];
		end;
	end;
	return (aP);
	finish PermuteMatrix;

	run;
	quit;
*	Put in character descriptors;
data &dsResults;
	merge &dsResults &dsResults.Values;
	drop iRow;
	*	The following statement is necessary to give mixed upper and lowercase variable names;
	rename CoeffRowVar=CoeffRowVar Type=Type;
	run;
%mend RandomWalkAnalysis;

/*
*_________________________________________________________________________________________

	Test code
__________________________________________________________________________________________;

options nocenter mprint sasautos=(sasautos,".") mautosource;
*	Generate some correlated variates;
title 'Decomposition with simulated data';
data;
	do i=1 to 100;
		p1 = normal(123);
		p2 = p1 + .1*normal(456);
		p3 = p2 + .1*normal(789);
		output;
		drop i;
	end;
proc corr cov outp=acov;
data acov;
	set acov;
	where _type_="COV";
proc print data=acov;
data coeff;
	p1 = .4;
	p2 = .8;
	p3 = .2;
	output;
	p3 = .25;
	output;
	run;
%RandomWalkAnalysis(rwt,acov,coeff);
proc print data=rwt;
	run;
*/
