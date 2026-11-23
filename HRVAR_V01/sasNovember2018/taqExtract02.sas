*________________________________________________________________________________________________________________________________
	taqExtract02
	Build csv files from wrds datasets.
	THIS RUNS LOCALLY
;;
dm log 'clear';
x 'cd c:\active\pdl\sas';
options nocenter number mautosource ls=200 label noxwait;
libname temp "c:\temp";
libname wrds "c:\active\pdl\wrdsSasDatasets";
%let wrdsDir=c:\active\pdl\wrdsSasDatasets;
options append=(sasautos=("c:\active\CLS\CLSsas\sasmacros",'.'));

goptions reset;
symbol1 i=join v=none;
%macro plot(dsIn,plotVar,symbol=);
proc tabulate data=&dsIn (where=(
	%if %length(&symbol)>0 %then sym_root="&symbol" and ;
	timeMid>'9:35:00't and timeMid<'15:55't)) missing format=comma11.;
	var &plotVar;
	class date sym_root;
	table sym_root*date,&plotvar*(n nmiss);
	table sym_root*date, n*f=comma7. (min max mean std)*&plotvar*f=dollar8.2;
	run;
proc gplot data=&dsIn (where=(
	%if %length(&symbol)>0 %then sym_root="&symbol" and ;
	timeMid>'9:35:00't and timeMid<'15:55't));
	plot &plotVar*timeMid=1;
	by date sym_root;
	format timeMid time.;
	run;
	quit;
%mend;
%*plot(wrds.part_nbb,nbb);
%*plot(wrds.part_nbo,nbo);
%plot(wrds.lex_nbb,nbb,symbol=NVDA)
%runquit;


/*proc export data=wrds.part_nbb outfile="&wrdsDir.\part_nbb.csv" dbms=csv replace;*/
/*proc export data=wrds.part_nbo outfile="&wrdsDir.\part_nbo.csv" dbms=csv replace;*/
/*proc export data=wrds.sip_nbb outfile="&wrdsDir.\sip_nbb.csv" dbms=csv replace;*/
/*proc export data=wrds.sip_nbo outfile="&wrdsDir.\sip_nbo.csv" dbms=csv replace;*/

proc export data=wrds.lex_nbb outfile="&wrdsDir.\lex_nbb.csv" dbms=csv replace;
proc export data=wrds.lex_nbo outfile="&wrdsDir.\lex_nbo.csv" dbms=csv replace;
proc export data=wrds.exl_nbb outfile="&wrdsDir.\exl_nbb.csv" dbms=csv replace;
proc export data=wrds.exl_nbo outfile="&wrdsDir.\exl_nbo.csv" dbms=csv replace;

/*proc export data=wrds.ctDark outfile="&wrdsDir.\ctDark.csv" dbms=csv replace;*/
/*proc export data=wrds.ctLit outfile="&wrdsDir.\ctLit.csv" dbms=csv replace;*/
run;
