/*
	Multiscale and bridged VECMs
*/
dm log 'clear';
x 'cd c:\active\pdl\sasNovember2018';
options nocenter number mautosource append=(sasautos=('.',"./sasmacros","./sasmacros/rw")) 
	orientation=portrait ls=200 label noxwait;
libname this '.';

options nomprint nosymbolgen;
title;
footnote;
%let i0=1;

symbol1 v=none i=join color=blue w=2 l=1;
symbol2 v=none i=join color=blue w=2 l=2;
symbol3 v=dot  i=join color=blue w=2 l=3;
symbol4 v=squareFilled  i=join color=blue w=2 l=20;
symbol5 v=none i=join color=blue w=4 l=2;
symbol6 v=none i=join color=red  w=2 l=2;

goptions reset;
goptions ftext="albany amt" htext=11pt htitle=11pt hsize=6.5in vsize=5in device='emf';
axis1 logbase=10 label=('Time ("seconds")');
axis2 label=none;
legend1 label=none;

ods rtf close;
ods rtf file="./odsVarsim04/varsim04_&i0..rtf" startpage=no;
*_____________________________________________________________________________________________________________________________________________________________
	Simulation of two-price VECM with linear coefficients
;;
%macro vecmSim1(dsOut, nStartup, nT, gam, pMax, phi11, phi22);
data _null_;
	put "vecSim1. dsOut=&dsOut";
	file print;
	put 'varsim04: vecSim1';
	put "dsOut=&dsOut, nStartup=&nStartup, nT=&nT, gam=&gam pMax=&pMax";
	put "phi11=&phi11";
	put "phi22=&phi22";
	run;
data &dsOut;
	array dpL(2,&pMax) _temporary_;
	array dp(2);
	array p(2);
	array u(2) _temporary_;
	array gam(2) _temporary_ &gam;
	array phi1(&pMax) _temporary_ &phi11;
	array phi2(&pMax) _temporary_ &phi22;
	do i=1 to 2;
		dp(i) = 0;
		p(i) = 0;
		do j=1 to &pMax;
			dpl(i,j)=0;
		end;
	end;
	do t=-&nStartup to &nT;
		u(1) = normal(456);
		u(2) = normal(567);
		cv = p1 - p2;
		dp1 = gam(1)*cv;
		dp2 = gam(2)*cv;
		do j=1 to &pMax;
			dp1 = dp1 + phi1(j)*dpL(1,j);
			dp2 = dp2 + phi2(j)*dpL(2,j);
		end;
		do i=1 to 2;
			dp(i) = dp(i) + u(i);
			p(i) = p(i) + dp(i);
		end;
		if t>0 then output;
		drop i j;
		do i=1 to 2;
			do j=&pMax to 2 by -1;
				dpl(i,j) = dpl(i,j-1);
			end;
			dpl(i,1) = dp(i);
		end;
	end;
	keep p1 p2 t;
	run;
/*goptions reset;*/
/*symbol1 i=join v=none color=blue;*/
/*symbol2 i=join v=none color=red;*/
/*proc gplot data=&dsOut;*/
/*	plot p1*t=1 p2*t=2 / overlay;*/
/*	run;*/
/*	quit;*/
%mend vecmSim1;
%vecmSim1(temp.dsvecm1, 1000, 1000000, (-.01 .02), 100, (0.1 9*.02 90*0.005), (0.1 9*.02 90*0.005));

%let irfMax=1000;
*_____________________________________________________________________________________________________________________________________________________________
	Estimation at time scale 1
;;
%macro ts1;
title "VECM estimation at time scale 1";
footnote "Time scale 1";
proc model data=temp.dsvecm1 plots=(none) outmodel=temp.model1;
	endogenous p1 p2;
	cv = lag(p1 - p2);
	p1 = alpha1*cv + zlag(p1) + phi111*zlag(dif(p1)) + phi119*( lag2(p1)-lag11(p1) ) + phi1190*( lag11(p1)-lag101(p1) );
	p2 = alpha2*cv + zlag(p2) + phi221*zlag(dif(p2)) + phi229*( lag2(p2)-lag11(p2) ) + phi2290*( lag11(p2)-lag101(p2) );
	ods output CovResiduals=temp.cov1;
	fit p1 p2 / outest=est1 covs;
	run;
%*printds(data=temp.cov1);
title "IRF computation at time scale 1";
data temp.estMultiple1;	*	Make copies of estimates to use in irf calculations;
	set est1 (in=in1) est1 (in=in2);
	if in1 then iVar=1;
	if in2 then iVar=2;
	run;
data irf;
	array p(2);
	do iVar=1 to 2;
		do t=-100 to &irfMax;	*	include slots for startup lags!;
			drop i;
			if t<=0 then do i=1 to 2;
				p(i) = 0;
			end;
			if t=0 then p(iVar)=&i0;
			if t>0 then do i=1 to 2;
				p(i) = .;
			end;
			output;
		end;
	end;
	run;
proc model model=temp.model1 noprint;
	solve / estdata=temp.estMultiple1 data=irf out=temp.irf1 (drop=_lag_ _errors_ _type_ _mode_ where=(t>=0)) outlags forecast time=t;
	by iVar;
	run;
%printds(data=temp.irf1, where=t<=20);
proc gplot data=temp.irf1 (where=(t>0)) uniform;
	plot p1*t=1 p2*t=2 / overlay legend=legend1 haxis=axis1 vaxis=axis2 autovref;
	by iVar;
	run;
	quit;
title "Random walk analysis at time scale 1";
data maSum1;
	set temp.irf1;
	by iVar;
	if last.iVar then output;
	run;
%*printds(data=maSum1);
proc transpose data=maSum1 out=temp.maSumT1 (drop=_:) prefix=p;
	id iVar;
	var p1 p2;
	run;
%*printds(data=temp.maSumT1);
%randomWalkAnalysis(rwOut,temp.cov1,temp.maSumT1);
%mend ts1;
%ts1;

*_____________________________________________________________________________________________________________________________________________________________
	Truncated estimation at time scale 1
;;
%macro ts1t;
title "Truncated estimation at time scale 1";
footnote "Truncated estimation, Time scale 1";
proc model data=temp.dsvecm1 plots=(none) outmodel=temp.model1Trunc;
	endogenous p1 p2;
	cv = lag(p1 - p2);
	p1 = alpha1*cv + zlag(p1) + phi111*zlag(dif(p1)) + phi119*( lag2(p1)-lag11(p1) );
	p2 = alpha2*cv + zlag(p2) + phi221*zlag(dif(p2)) + phi229*( lag2(p2)-lag11(p2) );
	ods output CovResiduals=temp.cov1trunc;
	fit p1 p2 / outest=est1trunc covs;
	run;
%*printds(data=cov1t);
title "(truncted) IRF computation at time scale 1";
data temp.estMultiple1Trunc;	*	Make copies of estimates to use in irf calculations;
	set est1trunc (in=in1) est1trunc (in=in2);
	if in1 then iVar=1;
	if in2 then iVar=2;
	run;
data irf;
	array p(2);
	do iVar=1 to 2;
		do t=-100 to &irfMax;	*	include slots for startup lags!;
			drop i;
			if t<=0 then do i=1 to 2;
				p(i) = 0;
			end;
			if t=0 then p(iVar)=&i0;
			if t>0 then do i=1 to 2;
				p(i) = .;
			end;
			output;
		end;
	end;
	run;
proc model model=temp.model1Trunc noprint;
	solve / estdata=temp.estMultiple1Trunc data=irf out=temp.irf1Trunc (drop=_lag_ _errors_ _type_ _mode_ where=(t>=0)) outlags forecast time=t;
	by iVar;
	run;
proc gplot data=temp.irf1Trunc (where=(t>0)) uniform;
	plot p1*t=1 p2*t=2 / overlay legend=legend1 haxis=axis1 vaxis=axis2 autovref;
	by iVar;
	run;
	quit;
title "(truncated) Random walk analysis at time scale 1";
data maSum1;
	set temp.irf1Trunc;
	by iVar;
	if last.iVar then output;
	run;
%*printds(data=maSum1);
proc transpose data=maSum1 out=temp.maSumT1t (drop=_:) prefix=p;
	id iVar;
	var p1 p2;
	run;
%*printds(data=temp.maSumT1);
%randomWalkAnalysis(rwOut,temp.cov1trunc,temp.maSumT1t);
%mend ts1t;
%ts1t;

*_____________________________________________________________________________________________________________________________________________________________
	Estimation at time scale 10
;;
%macro ts10;
title "VECM estimation at time scale 10";
footnote "Time scale 10";
data dsvecm10;
	set temp.dsvecm1;
	where mod(t,10)=0;
	run;
proc model data=dsvecm10 plots=(none) outmodel=temp.model10;
	endogenous p1 p2;
	cv = lag(p1 - p2);
	p1 = alpha1*cv + zlag(p1) + phi111*zlag(dif(p1)) + phi119*( lag2(p1)-lag11(p1) );
	p2 = alpha2*cv + zlag(p2) + phi221*zlag(dif(p2)) + phi229*( lag2(p2)-lag11(p2) );
	ods output CovResiduals=temp.cov10;
	fit p1 p2 / outest=est10 covs;
	run;
data temp.estMultiple10;	*	Make copies of estimates to use in irf calculations;
	set est10 (in=in1) est10 (in=in2);
	if in1 then iVar=1;
	if in2 then iVar=2;
	run;
data irf;
	array p(2);
	do iVar=1 to 2;
		drop i;
		do i=1 to 2; p(i)=0; end;
		do t=-100 to -1;	*	include slots for startup lags!;
			output;
		end;
		t = 0;
		p(iVar)=&i0;
		output;
		do i=1 to 2; p(i)=.; end;
		do t=10 to &irfMax by 10;
			output;
		end;
	end;
	run;
proc model model=temp.model10 noprint;
	solve / estdata=temp.estMultiple10 data=irf out=temp.irf10 (drop=_lag_ _errors_ _type_ _mode_ where=(t>=0)) outlags forecast time=t;
	title "IRF calculation at time scale 10";
	by iVar;
	run;
proc gplot data=temp.irf10 (where=(t>0)) uniform;
	plot p1*t=1 p2*t=2 / overlay legend=legend1 haxis=axis1 vaxis=axis2 autovref;
	by iVar;
	run;
	quit;
title "Random walk analysis at time scale 10";
data maSum10;
	set temp.irf10;
	by iVar;
	if last.iVar then output;
	run;
%*printds(data=maSum10);
proc transpose data=maSum10 out=temp.maSumT10 (drop=_:) prefix=p;
	id iVar;
	var p1 p2;
	run;
%*printds(data=maSumT);
%randomWalkAnalysis(rwOut,temp.cov10,temp.maSumT10);
data _null_;
	file print;
	put "To rescale the rw volatility per 10 seconds to rw volatility per second, divide by sqrt(10)."
	put "8.5795 / sqrt(10) = 2.713 [the value reported in the table";
	run;

%mend ts10;
%ts10;

*_____________________________________________________________________________________________________________________________________________________________
	Merge all three irfs
;;
data temp.irf3;
	merge temp.irf1 temp.irf1Trunc (rename=(p1=p1t p2=p2t)) temp.irf10 (rename=(p1=p1_10 p2=p2_10));
	by ivar t;
	if t>0 then output;
	run;
*_____________________________________________________________________________________________________________________________________________________________
	Graph all three [Not used in paper]
;;
%macro graph3;
%annomac;
data dsAnnotate;
	xsys='1';
	hsys='D';
	%label(1,1,"testing",black,0,0,11,'NONE',1);
	run;
data dsAnnotate;	*	duplicate to cover both shocks;
	set dsAnnotate (in=in1) dsAnnotate (in=in2);
	if in1 then iVar=1;
	if in2 then iVar=2;
run;
%printds(data=dsAnnotate);
proc gplot data=temp.irf3 (where=(t>0)) uniform;
	plot p1*t=1 p1t*t=3 p1_10*t=5 / overlay legend=legend1 haxis=axis1 vaxis=axis2 autovref annotate=dsAnnotate;
	by iVar;
	run;
	quit;
%mend;

*_____________________________________________________________________________________________________________________________________________________________
	Bridged IRF (based on full specifications, not used in paper)
;;
%macro bridge1;
footnote "Bridged IRF calculation";
title "Bridged IRF calculation";
data b1;
	set temp.irf1;
	by iVar;
	if t=10 then output;
	run;
%printds(data=b1);
data irf;
	set b1;
	by iVar;
	array p(2);
	if first.iVar then do;
		p1Save = p1;
		p2Save = p2;
		drop p1Save p2Save;
		do t=-100 to -1;
			p1 = 0.;
			p2 = 0.;
			output;
		end;
		t = 0;
		p(iVar) = &i0;
		output;
		t = 10;
		p1 = p1Save;
		p2 = p2Save;
		output;
		do t=20 to 1000 by 10;
			p1=.;
			p2=.;
			output;
		end;
	end;
	run;
%printds(data=irf, where=t>=-1 and t<=100);
proc model model=temp.model10 noprint;
	solve / estdata=temp.estMultiple10 data=irf out=irfBridge (drop=_lag_ _errors_ _type_ _mode_ where=(t>=0)) outlags forecast time=t;
	by iVar;
	run;
%printds(data=irfBridge, obs=50);
data irf1b;
	merge temp.irf1 irfBridge (rename=(p1=p1Bridge p2=p2Bridge));
	by iVar t;
	run;
%printds(data=irf1b, where=t>=0 and t<=20);
proc gplot data=irf1b (where=(t>0)) uniform;
	plot p1*t=1 p2*t=2 p1Bridge*t=3 p2Bridge*t=4 / overlay legend=legend1 haxis=axis1 vaxis=axis2 autovref;
	by iVar;
	run;
	quit;
data maSumBridge;
	set irfBridge;
	by iVar;
	if last.iVar then output;
	run;
%printds(data=maSumBridge);
proc transpose data=maSumBridge out=maSumTBridge (drop=_:) prefix=p;
	id iVar;
	var p1 p2;
	run;
%printds(data=maSumTBridge);
title "Bridged random walk analysis";
%randomWalkAnalysis(rwOut,temp.cov10,maSumTBridge);
%mend;
%*bridge1;

*_____________________________________________________________________________________________________________________________________________________________
	Bridged IRF (with truncated estimation at time scale 1, used in paper)
;;
%macro bridge2;
footnote "Truncated bridged IRF calculation";
title "Truncated bridged IRF calculation";
data b1;
	set temp.irf1Trunc;
	by iVar;
	if t=10 then output;
	run;
%*printds(data=b1);
data irf;
	set b1;
	by iVar;
	array p(2);
	if first.iVar then do;
		p1Save = p1;
		p2Save = p2;
		drop p1Save p2Save;
		do t=-100 to -1;
			p1 = 0.;
			p2 = 0.;
			output;
		end;
		t = 0;
		p(iVar) = &i0;
		output;
		t = 10;
		p1 = p1Save;
		p2 = p2Save;
		output;
		do t=20 to &irfMax by 10;
			p1=.;
			p2=.;
			output;
		end;
	end;
	run;
%*printds(data=irf, where=t>=-1 and t<=100);
proc model model=temp.model10 noprint;
	solve / estdata=temp.estMultiple10 data=irf out=temp.irfBrTr (drop=_lag_ _errors_ _type_ _mode_ where=(t>=0)) outlags forecast time=t;
	by iVar;
	run;
%printds(data=temp.irfBrTr, obs=50);
/*data irf1b;*/
/*	merge temp.irf1 temp.irf1Trunc (where=(t>0 and t<=10) rename=(p1=p1Tr p2=p2Tr)) irfBrTr (rename=(p1=p1BrTr p2=p2BrTr));*/
/*	by iVar t;*/
/*	run;*/
/*%printds(data=irf1b, where=t>=0 and t<=20);*/
/*proc gplot data=irf1b (where=(t>0)) uniform;*/
/*	plot p1*t=1 p2*t=2 p1Tr*t=3 p2Tr*t=4 p1BrTr*t=5 p2BrTr*t=6 / overlay legend=legend1 haxis=axis1 vaxis=axis2 autovref;*/
/*	by iVar;*/
/*	run;*/
/*	quit;*/
data maSumBridge;
	set temp.irfBrTr;
	by iVar;
	if last.iVar then output;
	run;
%printds(data=maSumBridge);
proc transpose data=maSumBridge out=temp.maSumTBridge (drop=_:) prefix=p;
	id iVar;
	var p1 p2;
	run;
%printds(data=temp.maSumTBridge);
title "Truncated bridged random walk analysis";
%randomWalkAnalysis(rwOut,temp.cov1trunc,temp.maSumTBridge);
%mend;
%bridge2;
data temp.irf4;
	merge temp.irf3 temp.irfBrTr (rename=(p1=p1tb p2=p2tb));
	by iVar t;
	run;
%annomac;
data dsAnnotate;
	xsys='2';
	ysys='2';
	hsys='D';
	%*label(1,2,"testing",black,0,0,11,'NONE',7);
	%label(1, 2.1,"Estim",black,0,0,11,'NONE',>);
/*	%label(100,1.1,"10 est & irf",black,0,0,11,'NONE',7);*/
/*	%label(90,1.2,"100 est & irf",black,0,0,11,'NONE',7);*/
/*	%label(20,2.5,"10 est, 100 est, bridged irf",black,0,0,11,'NONE',7);*/
	run;
/*data dsAnnotate;	*	duplicate to cover both shocks;*/
/*	set dsAnnotate (in=in1) dsAnnotate (in=in2);*/
/*	if in1 then iVar=1;*/
/*	if in2 then iVar=2;*/
/*run;*/
%printds(data=dsAnnotate);
symbol;
symbol1 v=none i=join color=blue w=3 l=1;
symbol2 v=none i=join color=blue w=1 l=1;
symbol3 v=none  i=join color=blue w=3 l=2;
symbol4 v=none 	i=join color=blue w=1 l=3;
proc gplot data=temp.irf4 (where=(t>0 and iVar=1));
	plot p1*t=1 p1t*t=2 p1_10*t=3 p1tb*t=4 / overlay haxis=axis1 vaxis=axis2 autovref legend=legend1;
	label p1='Short and long' p1t='Short' p1_10='Long' p1tb='Bridged';
	by iVar;
	run;
	quit;

title "Random walk analysis: long and short";
%randomWalkAnalysis(rwOut1,temp.cov1,temp.maSumT1,AllRows=0);
%printds(data=rwOut1);
title "Random walk analysis: short";
%randomWalkAnalysis(rwOut1t,temp.cov1trunc,temp.maSumT1t,AllRows=0);
title "Random walk analysis: long";
%randomWalkAnalysis(rwOut10,temp.cov10,temp.maSumT10,AllRows=0);
title "Random walk analysis: bridged";
%randomWalkAnalysis(rwOutBridge,temp.cov1trunc,temp.maSumTBridge,AllRows=0);
data rwAll;
	set rwOut1 (in=in1) rwOut1t (in=in2) rwOut10 (in=in3) rwOutBridge (in=in4);
	if in1 then model="Short and long";
	if in2 then model="Short";
	if in3 then model="Long";
	if in3 then varRandwalk=varRandwalk/10;
	if in4 then model="Bridged";
	run;
proc tabulate data=rwAll order=data format=comma6.3;
	var varRandwalk p1 p2;
	class model type;
	table model='',mean=''*(varRandWalk p1*type p2*type) / nocellmerge;
	run;
ods rtf close;
