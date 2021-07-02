/*********************************************************************************************/
title1 'Exploring AD Incidence Definition';

* Author: PF;
* Purpose: 	Dementia Incidence Rates for all methods;

options compress=yes nocenter ls=150 ps=200 errors=5 mprint merror
	mergenoby=warn varlenchk=error dkricond=error dkrocond=error msglevel=i;
/*********************************************************************************************/

%include "../header";

%let minyear=2001;
%let maxyear=2018;

* merge to sample;
data dem_samp;
	merge &outlib..adrdinc_&minyear._&maxyear. (in=a) &outlib..samp_3yrffsptd_0618 (in=b);
	by bene_id;
	keep bene_id birth_date death_date race_bg sex scen: drop: age: first: drop: insamp:;
run;

%macro inc_dropv(byear,eyear,samp=);
	
%do year=&byear %to &eyear;
	%let nextyear=%eval(&year+1);
	%let lastyear=%eval(&year+1);

	data inc_dropv&year.;
		set dem_samp;
		
		if insamp&year.;
		
		death_yr=year(death_date);
		
		year=&year;
	
		if dropdx ne 1 then do;
			if .<year(scen_dx_inc)<&year then do;
				dx_atrisk=0;
				dx_priorv=1;
			end;
			else do;
				dx_inc=0;
				dx_atrisk=1;
				if death_yr=&year then dx_deathexit=1;
				else if insamp&nextyear ne 1 then dx_sampexit=1;
				if year(scen_dx_inc)=&year then do;
					dx_inc=1;
					if find(scen_dx_vtype,"1") then dx_scen1=1;
					if find(scen_dx_vtype,"4") then do;
						dx_scen4=1; 
						death_dx1=1;
					end;
				end;
			end;
		end;
		
		if dropdxrx ne 1 then do;
			if .<year(scen_dxrx_inc)<&year then do;
				dxrx_atrisk=0;
				dxrx_priorv=1;
			end;
			else do;
				dxrx_inc=0;
				dxrx_atrisk=1;
				if death_yr=&year then dxrx_deathexit=1;
				else if insamp&nextyear ne 1 then dxrx_sampexit=1;
				if year(scen_dxrx_inc)=&year then do;
					dxrx_inc=1;
					if find(scen_dxrx_vtype,"1") and find(scen_dxrx_inctype,"1") then dxrx_scen1=1;
					if find(scen_dxrx_vtype,"2") or find(scen_dxrx_inctype,"2") then dxrx_scen2=1;
					if find(scen_dxrx_vtype,"4") then do;
						dxrx_scen4=1;
						if find(scen_dxrx_inctype,"1") then death_dxrx1=1;
						if find(scen_dxrx_inctype,"2") then death_dxrx2=1;
					end;
				end;
			end;
		end;

		if dropdxsymp ne 1 then do;
			if .<year(scen_dxsymp_inc)<&year then do;
				dxsymp_atrisk=0;
				dxsymp_priorv=1;
			end;
			else do;
				dxsymp_atrisk=1;
				dxsymp_inc=0;
				if death_yr=&year then dxsymp_deathexit=1;
				else if insamp&nextyear ne 1 then dxsymp_sampexit=1;
				if year(scen_dxsymp_inc)=&year then do;
					dxsymp_inc=1;
					* counting only when contributes to final;
					if find(scen_dxsymp_vtype,"1") and find(scen_dxsymp_inctype,"1") then dxsymp_scen1=1;
					if find(scen_dxsymp_vtype,"3") or find(scen_dxsymp_inctype,"3") then dxsymp_scen3=1;
					if find(scen_dxsymp_vtype,"4") then do;
						dxsymp_scen4=1;
						if find(scen_dxsymp_inctype,"1") then death_dxsymp1=1;
						if find(scen_dxsymp_inctype,"3") then death_dxsymp3=1;
					end;
				end;
			end;
		end;

		if dropdxrxsymp ne 1 then do;
			if .<year(scen_dxrxsymp_inc)<&year then do;
				dxrxsymp_atrisk=0;
				dxrxsymp_priorv=1;
			end;
			else do;
				dxrxsymp_inc=0;
				dxrxsymp_atrisk=1;
				if death_yr=&year then dxrxsymp_deathexit=1;
				else if insamp&nextyear ne 1 then dxrxsymp_sampexit=1;
				if year(scen_dxrxsymp_inc)=&year then do;
					dxrxsymp_inc=1;
					if find(scen_dxrxsymp_vtype,"1") and find(scen_dxrxsymp_inctype,"1") then dxrxsymp_scen1=1;
					if find(scen_dxrxsymp_vtype,"2") or find(scen_dxrxsymp_inctype,"2") then dxrxsymp_scen2=1;
					if find(scen_dxrxsymp_vtype,"3") or find(scen_dxrxsymp_inctype,"3") then dxrxsymp_scen3=1;
					if find(scen_dxrxsymp_vtype,"4") then do;
						dxrxsymp_scen4=1;
						if find(scen_dxrxsymp_inctype,"1") then death_dxrxsymp1=1;
						if find(scen_dxrxsymp_inctype,"2") then death_dxrxsymp2=1;
						if find(scen_dxrxsymp_inctype,"3") then death_dxrxsymp3=1;
					end;
				end;
			end;
		end;

	run;
%end;

%macro stats(subgroup=,output=);
%do year=&byear %to &eyear;

	* Incident rate statistics;
	proc means data=inc_dropv&year noprint nway missing;
		class &subgroup year;
		var dx_inc dxrx_inc dxsymp_inc dxrxsymp_inc;
		output out=incstats_dropv&output.&year. (drop=_type_) mean(dx_inc dxrx_inc dxsymp_inc dxrxsymp_inc)=dx_inc_rate dxrx_inc_rate dxsymp_inc_rate dxrxsymp_inc_rate
		lclm(dx_inc dxrx_inc dxsymp_inc dxrxsymp_inc)=dx_inc_lclm dxrx_inc_lclm dxsymp_inc_lclm dxrxsymp_inc_lclm
		uclm(dx_inc dxrx_inc dxsymp_inc dxrxsymp_inc)=dx_inc_uclm dxrx_inc_uclm dxsymp_inc_uclm dxrxsymp_inc_uclm
		sum(dx_inc dxrx_inc dxsymp_inc dxrxsymp_inc)=dx_inc_total dxrx_inc_total dxsymp_inc_total dxrxsymp_inc_total
		sum(dx_atrisk dxrx_atrisk dxsymp_atrisk dxrxsymp_atrisk)=;
	run;
	
	* Age, death & sample exit statistics;
	proc means data=inc_dropv&year noprint nway missing;
		class &subgroup year;
		output out=othstats_dropv&output.&year (drop=_type_ _freq_) mean(age_beg&year)=age_beg_avg
		p25(age_beg&year)=age_beg_p25 median(age_beg&year)=age_beg_med p75(age_beg&year)=age_beg_p75
		sum(dx_deathexit dxrx_deathexit dxsymp_deathexit dxrxsymp_deathexit dx_sampexit dxrx_sampexit dxsymp_sampexit dxrxsymp_sampexit)=;
	run;
	                                                            
	* Break down of denominator drops;
	proc means data=inc_dropv&year noprint nway missing;
		class &subgroup year;
		output out=drops_dropv&output.&year (drop=_type_ _freq_)
		sum(dx_priorv dxrx_priorv dxsymp_priorv dxrxsymp_priorv)=;
	run;
	
	* Breakdown of numerator - not mutually exclusive, except for scenario 4;
	proc means data=inc_dropv&year noprint nway missing;
		class &subgroup year;
		output out=numerator_dropv&output.&year (drop=_type_ _freq_)
		sum(dx_scen1 dx_scen4 dxrx_scen2 dxrx_scen4 dxrx_scen1 dxsymp_scen1 dxsymp_scen3 dxsymp_scen4 dxrxsymp_scen1-dxrxsymp_scen4
		death_dx1 death_dxrx1-death_dxrx2 death_dxsymp1 death_dxsymp3 death_dxrxsymp1-death_dxrxsymp3)=;
	run;
%end;

***** Merging all statistics together;
data ADRDinc_stats_dropv&samp.&output.;
	format dx_inc_lclm dxrx_inc_lclm dxsymp_inc_lclm dxrxsymp_inc_lclm
		   dx_inc_uclm dxrx_inc_uclm dxsymp_inc_uclm dxrxsymp_inc_uclm;
	merge incstats_dropv&output.&byear-incstats_dropv&output.&eyear
		  drops_dropv&output.&byear-drops_dropv&output.&eyear
		  numerator_dropv&output.&byear-numerator_dropv&output.&eyear
		  othstats_dropv&output.&byear-othstats_dropv&output.&eyear;
	by &subgroup year;
run;

ods excel file="./output/adrdinc_ma_stats&output..xlsx";
proc print data=ADRDinc_stats_dropv&samp.&output.; run;
ods excel close;
%mend;

%stats(subgroup=,output=);
%stats(subgroup=race_bg,output=race);
%stats(subgroup=sex,output=sex);

%mend;

%inc_dropv(&minyear.,&maxyear.,samp=);