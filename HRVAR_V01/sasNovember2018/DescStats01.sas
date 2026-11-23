*________________________________________________________________________________________________________________________________
	DescStats01
	THIS RUNS LOCALLY
;;
dm log 'clear';
x 'cd c:\active\pdl\sasNovember2018';
options nocenter number mautosource ls=200 label noxwait;
libname temp "c:\temp";
libname wrds "..\wrdsSasDatasets";
libname agora "..\agora";
options append=(sasautos=(".\sasmacros",'.'));

/*data temp.cqOct3;*/
/*	set agora.taqq;*/
/*	where date='03oct2016'd;*/
/*	run;*/
/*data temp.ctOct3;*/
/*	set agora.taqt;*/
/*	where date='03oct2016'd;*/
/*	run;*/

proc format;
	value $exch (default=11) 'A'='NYSE MKT' 'B'='BX (NASDAQ)' 'C'='NSX' 'D'='FINRA ADF' 'I'='ISE'
		'J'='EDGA' 'K'='EDGX' 'M'='CHX' 'N'='NYSE' 'P'='NYSE ARCA' 'T'='NASDAQ' 'Q'='NASDAQ'
		'V'='IEX' 'W'='CBOE' 'X'='NASDAQ PSX' 'Y'='Bats BYX' 'Z'='Bats BZX';
	value znz (notsorted) 0='Zero' low-<0='Neg' 0<-high='Pos';
	run;
title 'trades';
data ct;
	set temp.ctOct3;
	where date='03Oct2016'd and (not (ex in('C','M'))) and part_time>='9:45't and part_time<='16:00't;
	by date sym_root;
	tValue = price*size;	*	million $;
	d = (time_m - part_time)*1000;
	run;
options orientation=portrait;
ods rtf file='./odsDescStats/tables234.rtf';
proc tabulate data=ct noseps missing order=formatted;
	class sym_root ex;
	var tValue d price size;
	table sym_root,n*f=comma10. mean*price*f=dollar8.2 mean*size*f=comma6.1 mean*tValue*f=dollar10.2 / nocellmerge misstext=' ' box='Table 2.';
	table sym_root*(ex all),n*f=comma10. (min p1 p10 p50 p90 p99 max)*d*f=comma10.2 / nocellmerge misstext=' ' box='Table 3A.';
	table sym_root*(ex all),(n*f=comma10. pctn<ex all>*f=6.1 sum*tValue*f=comma15. pctsum<ex all>*f=6.1*tValue) / nocellmerge misstext=' ' box='Table 4A.';
	format ex $exch.;
	run;
title 'quotes';
data cq;
	set temp.cqOct3;
	where date='03Oct2016'd and (not (ex in('C','M'))) and part_time>='9:45't and part_time<='16:00't;
	d = (time_m - part_time)*1000;
	run;
proc tabulate data=cq noseps missing order=formatted;
	class sym_root ex;
	var d;
	table sym_root,n*f=comma10.;
	table sym_root*(ex all),n*f=comma10. (min p1 p10 p50 p90 p99 max)*d*f=comma10.2 / nocellmerge misstext=' ' box='Table 3B.';
	table sym_root*(ex all),(n*f=comma10. pctn<ex all>*f=6.1) / nocellmerge misstext=' ' box='Table 4B.';
	format ex $exch.;
	run;
ods rtf close;
