/*
 *	printds: a utility macro to print out datasets
 *	28 Feb 2018. Title added.
 *
 *	Usage: %printds(data=, obs=20, contents=no, last=)
 *
 *	All parameters are keyword parameters; all have defaults
 *
 *	data			dataset or view (defaults to most recently-created dataset/view)
 *	obs				number of observations to print at start of dataset (default=20)
 *	contents	if contents=yes, then also prints description of dataset (using proc contents)
 *	last			last=xx --> the last xx observations in the dataset are also printed (default: not printed)
 *
 */

%macro printds(data=, obs=20, contents=no, last=, where=);
/*	Save title string	*/
%let TitleText=;
proc sql noprint;
	select Text into :TitleText TRIMMED  /* put the trimmed value into a macro */
  	from Dictionary.Titles
  	where Type="T" & Number=1;
quit;
%*put &=TitleText;
%if &data= %then %let data=&syslast;
%let dsid = %sysfunc(open(&DATA., IS));
%let anobs=%sysfunc(attrn(&dsid,anobs));
%if &anobs=1 %then %do;
	%let nobs = %sysfunc(attrn(&dsid, NLOBS));
	title "&data (&nobs observations)";
%end;
%else title "&data (unknown number of observations)"%str(;);
%if %upcase(&contents)=YES %then proc contents data=&data%str(; run;);
%if &obs>0 %then %do;
	proc print data=&data (obs=&obs);
		%if %length(&where)>0 %then where &where;;
		run;
%end;
%if &last>0 AND &anobs=1 %then %do;
	%let firstObs=%eval(&nobs-&last+1);
	%if &firstObs<1 %then %let firstObs=1;
	title "&data (last &last observations)";
	proc print data=&data (firstobs=&firstobs); run;
%end;
%if %length(&TitleText)>0 %then title "&TitleText";
%else title;;
%let rc=%sysfunc(close(&dsid));
%mend printds;

/*title "blah";*/
/*title;*/
/*data x;*/
/*	do i=1 to 10;*/
/*		output;*/
/*	end;*/
/*	run;*/
/*%printds;*/
/*proc print datra=x;*/
/*run;*/
