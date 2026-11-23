%macro getnobs(dsIn,n);
%global &n;
%let &n=-1;
%let dsid = %sysfunc(open(&dsIn., IS));
%let anobs=%sysfunc(attrn(&dsid,anobs));
%if &anobs=1 %then %let &n = %sysfunc(attrn(&dsid, NLOBS));
%let rc=%sysfunc(close(&dsid));
%mend;
