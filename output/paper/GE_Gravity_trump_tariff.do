****************************************************************
*This file contains the Stata code for the GE Gravity Analysis.*
****************************************************************
clear all
set mem 1g
set more off
clear matrix
set memory 500m
set matsize 8000
set maxvar 30000
cd "/home/bepp/hchulkim/Documents/gravity/Gravity_data"
cap log close
log using ge_gravity, text replace

******************
*I. Prepapre Data*
******************
use gravity_data, clear
***************************************************************
*1. Keep only latest year. Add and create additional variables*
***************************************************************
keep if year==2006
gen BRDR=1
replace BRDR=0 if exporter==importer
gen LN_DIST=ln(DIST)
*******************************
*2. Create aggregate variables*
*******************************
bysort exporter: egen output=sum(trade)
bysort importer: egen expndr=sum(trade)
****************************************
*3. Chose a country for reference group*
****************************************
replace exporter="ZZZ" if exporter=="DEU"
replace importer="ZZZ" if importer=="DEU"
gen expndr_deu0=expndr if importer=="ZZZ"
egen expndr_deu=mean(expndr_deu0)
*************************
*4. Create Fixed Effects*
*************************
qui tab exporter, gen(exp_fe_)
qui tab importer, gen(imp_fe_)
******************************
*5. Set additional parameters*
****************************** 
gen double phi=expndr/output if exporter==importer
scalar sigma=5
save ge_ppml_data, replace

************************
*II. GE Border Analysis*
************************
*************************************
*Step 1: Estimate `Baseline' Gravity*
*************************************
*********************************************************************************
*NOTE: This is the first place, where you may make changes, e.g., to account for*
*the border between Canada (CAN) and the US (USA) at the estimation stage.      *
*********************************************************************************
ppml trade LN_DIST CNTG BRDR exp_fe_* imp_fe_1-imp_fe_68, iter(30) noconst
predict trade_0_pred, mu
******************************************
*Step 1.a: Conscrut `Baseline' GE Indexes*
******************************************
forvalues i=1(1)68{
qui replace exp_fe_`i'=exp_fe_`i'*exp(_b[exp_fe_`i'])
qui replace imp_fe_`i'=imp_fe_`i'*exp(_b[imp_fe_`i'])
}
qui replace exp_fe_69=exp_fe_69*exp(_b[exp_fe_69])
egen all_exp_fes_0=rowtotal(exp_fe_1-exp_fe_69)
egen all_imp_fes_0=rowtotal(imp_fe_1-imp_fe_69) 
gen output_BLN=output
gen expndr_BLN=expndr
gen omr_BLN=output_BLN*expndr_deu/(all_exp_fes_0)
gen imr_BLN=expndr_BLN/(all_imp_fes_0*expndr_deu)
******************************************************************************
*NOTE: This is the second place, where you should consider making changes, so* 
*that the baseline trade costs are consistent with the estimating equation.  *
******************************************************************************
gen t_ij_BLN=exp(_b[LN_DIST]*LN_DIST+_b[CNTG]*CNTG+_b[BRDR]*BRDR)
gen trade_BLN=(output_BLN*expndr_BLN*t_ij_BLN)/(imr_BLN*omr_BLN)
gen double exp_BLN=trade_BLN if exporter!=importer
bysort exporter: egen tot_exp_BLN=sum(exp_BLN)
****************************************
*Step 2: Define Counterfactual Scenario*
****************************************
*********************************************************************************
*NOTE: This is the third and final place, where you should make changes. This is*
*where you will define the counterfactual vector of direct bilateral trade costs*
*********************************************************************************
gen t_ij_base = exp(_b[LN_DIST]*LN_DIST + _b[CNTG]*CNTG + _b[BRDR]*BRDR)
gen TARIFF=1

replace TARIFF=1.6 if importer == "USA" && exporter == "CHN"
replace TARIFF=1.6 if importer == "CHN" && exporter == "USA"
replace TARIFF=1.1 if importer == "USA" & exporter != "CHN" & exporter != "USA"
replace TARIFF=1.1 if exporter == "USA" & importer != "CHN" & importer != "USA"
*replace TARIFF=1.1 if importer == "USA" && exporter == "USA"

gen t_ij_CFL=t_ij_base *(TARIFF^(1-sigma))
gen t_ij_CFL_1=log(t_ij_CFL)
****************************************
*Step 3: Estimate `Conditional' Gravity*
****************************************
drop exp_fe_* imp_fe_*
qui tab exporter, gen(exp_fe_)
qui tab importer, gen(imp_fe_)
ppml trade exp_fe_* imp_fe_1-imp_fe_68, iter(30) noconst offset(t_ij_CFL_1)
predict trade_1_pred, mu
**********************************************
*Step 3.b: Construct `Conditional' GE Indexes*
**********************************************
forvalues i=1(1)68{
qui replace exp_fe_`i'=exp_fe_`i'*exp(_b[exp_fe_`i'])
qui replace imp_fe_`i'=imp_fe_`i'*exp(_b[imp_fe_`i'])
}
qui replace exp_fe_69=exp_fe_69*exp(_b[exp_fe_69])
egen all_exp_fes_1=rowtotal(exp_fe_1-exp_fe_69)
egen all_imp_fes_1=rowtotal(imp_fe_1-imp_fe_69) 
gen omr_CDL=output*expndr_deu/(all_exp_fes_1)
gen imr_CDL=expndr/(all_imp_fes_1*expndr_deu)
gen output_CDL=output
gen expndr_CDL=expndr
gen trade_CDL=(output_CDL*expndr_CDL*t_ij_CFL)/(imr_CDL*omr_CDL)
gen double exp_CDL=trade_CDL if exporter!=importer
bysort exporter: egen tot_exp_CDL=sum(exp_CDL)
*******************************************
*Step 4: Estimate `Full Endowment' Gravity*
*******************************************
*************************************
*Step 4.a: Prepare Initial Variables*
*************************************
gen expndr_1=expndr
gen double expndr_deu_1=expndr_deu
gen double temp1=all_exp_fes_0 if exporter==importer
bysort importer: egen double all_exp_fes_0_imp=mean(temp1)
gen double temp2=all_exp_fes_1 if exporter==importer
bysort importer: egen double all_exp_fes_1_imp=mean(temp2)
gen double p_FLL_exp_1=(all_exp_fes_1/all_exp_fes_0)^(1/(1-sigma))
gen double p_FLL_imp_1=(all_exp_fes_1_imp/all_exp_fes_0_imp)^(1/(1-sigma))
gen double p_FLL_exp_0=0
gen double imr_FLL_1=expndr_1/(all_imp_fes_1*expndr_deu_1)
gen double imr_FLL_ch_1=1
gen double omr_FLL_1=output*expndr_deu_1/(all_exp_fes_1)
gen double omr_FLL_ch_1=1
drop *temp*
******************************************************
*Step 4.b: Iterate to obtain `Full Endowment' indexes*
******************************************************
local i=3
local diff_all_exp_fes_sd=1
local diff_all_exp_fes_max=1
while (`diff_all_exp_fes_sd'>0.01) | (`diff_all_exp_fes_max'>0.01) {
		*******************************************************
		*i. Create new dependent variable and estimate gravity*
		*******************************************************
		gen double trade_`=`i'-1'=trade_`=`i'-2'_pred*p_FLL_exp_`=`i'-2'*p_FLL_imp_`=`i'-2'/(omr_FLL_ch_`=`i'-2'*imr_FLL_ch_`=`i'-2')
		drop exp_fe_* imp_fe_*
		qui tab exporter, gen (exp_fe_)
		qui tab importer, gen (imp_fe_)
		capture ppml trade_`=`i'-1' exp_fe_* imp_fe_*, offset(t_ij_CFL_1) noconst iter(30) 
		predict trade_`=`i'-1'_pred, mu
		*********************************
		*ii. Update output & expenditure*
		*********************************
		bysort exporter: egen double output_`=`i'-1'=total(trade_`=`i'-1'_pred)
		gen double expndr_temp_`=`i'-1'=phi*output_`=`i'-1' if exporter==importer
		bysort importer: egen double expndr_`=`i'-1'=mean(expndr_temp_`=`i'-1')
		gen double expndr_deu0_`=`i'-1'=expndr_`=`i'-1' if importer=="ZZZ"
		egen double expndr_deu_`=`i'-1'=mean(expndr_deu0_`=`i'-1')
		***********************************************
		*iii. Update factory-gate prices & resistances*
		***********************************************
		forvalues j=1(1)68{
		qui replace exp_fe_`j'=exp_fe_`j'*exp(_b[exp_fe_`j'])
		qui replace imp_fe_`j'=imp_fe_`j'*exp(_b[imp_fe_`j'])
		}
		qui replace exp_fe_69=exp_fe_69*exp(_b[exp_fe_69])
		egen all_exp_fes_`=`i'-1'=rowtotal(exp_fe_1-exp_fe_69)
		egen all_imp_fes_`=`i'-1'=rowtotal(imp_fe_1-imp_fe_69) 
		gen double temp=all_exp_fes_`=`i'-1' if exporter==importer
		bysort importer: egen double all_exp_fes_`=`i'-1'_imp=mean(temp)
		gen double p_FLL_exp_`=`i'-1'=((all_exp_fes_`=`i'-1'/all_exp_fes_`=`i'-2')/(expndr_deu_`=`i'-1'/expndr_deu_`=`i'-2'))^(1/(1-sigma))
		gen double p_FLL_imp_`=`i'-1'=((all_exp_fes_`=`i'-1'_imp/all_exp_fes_`=`i'-2'_imp)/(expndr_deu_`=`i'-1'/expndr_deu_`=`i'-2'))^(1/(1-sigma))
		gen double omr_FLL_`=`i'-1'=(output_`=`i'-1'*expndr_deu_`=`i'-1')/all_exp_fes_`=`i'-1'
		gen double omr_FLL_ch_`=`i'-1'=omr_FLL_`=`i'-1'/omr_FLL_`=`i'-2'
		gen double imr_FLL_`=`i'-1'=expndr_`=`i'-1'/(all_imp_fes_`=`i'-1'*expndr_deu_`=`i'-1')
		gen double imr_FLL_ch_`=`i'-1'=imr_FLL_`=`i'-1'/imr_FLL_`=`i'-2'
		****************************************************************************
		*iv. Iterate Until Convergence: Until change in factory-gate prices is zero*
		****************************************************************************
		gen double diff_p_FLL_exp_`=`i'-1'=p_FLL_exp_`=`i'-2'-p_FLL_exp_`=`i'-3'
		sum diff_p_FLL_exp_`=`i'-1'
		local diff_all_exp_fes_sd=r(sd)
		local diff_all_exp_fes_max=abs(r(max))
		local i=`i'+1
		drop temp* 
}
*************************************************
*Step 4.c: Construct `Full Endowment' GE Indexes*
*************************************************
forvalues j=1(1)68{
qui replace imp_fe_`j'=imp_fe_`j'*exp(_b[imp_fe_`j'])
}
gen double p_FLL=((all_exp_fes_`=`i'-2'/all_exp_fes_0)/(expndr_deu_`=`i'-2'/expndr_deu))^(1/(1-sigma))
gen output_FLL=p_FLL*output
gen double expndr_FLL_temp=phi*output_FLL if exporter==importer
bysort importer: egen double expndr_FLL=mean(expndr_FLL_temp)
gen double omr_FLL=output_FLL*expndr_deu_`=`i'-2'/(all_exp_fes_`=`i'-2')
gen double imr_FLL=expndr_`=`i'-2'/(all_imp_fes_`=`i'-2'*expndr_deu_`=`i'-2')
gen trade_FLL=(output_FLL*expndr_FLL*t_ij_CFL)/(imr_FLL*omr_FLL)
gen double exp_FLL=trade_FLL if exporter!=importer
bysort exporter: egen tot_exp_FLL=sum(exp_FLL)
save full_static_all, replace
**************************************
*Step 5: Construct Percentage Changes*
**************************************
******************************
*a. On Export/Production Side*
******************************
use full_static_all, clear
collapse p_FLL output* expndr* tot_exp*, by(exporter)
gen p_FLL_ch=(p_FLL-1)*100
gen tot_exp_FLL_ch=(tot_exp_FLL-tot_exp_BLN)/tot_exp_BLN*100
rename exporter country
save prod_all, replace
*******************************
*b. On Import/Consumption Side*
*******************************
use full_static_all, clear
collapse imr_FLL imr_BLN imr_CDL, by(importer)
rename importer country
gen imr_FLL_ch=(imr_FLL^(1/(1-sigma))-imr_BLN^(1/(1-sigma)))/imr_BLN^(1/(1-sigma))*100
save cons_all, replace
************************************
*c. Combine and Export to Tex/Excel*
************************************
joinby country using prod_all
gen rGDP_BLN=output_BLN/(imr_BLN^(1/(1-sigma)))
gen rGDP_FLL=output_FLL/(imr_FLL^(1/(1-sigma)))
gen rGDP_FLL_ch=(rGDP_FLL-rGDP_BLN)/rGDP_BLN*100
keep country tot_exp_FLL_ch p_FLL_ch imr_FLL_ch rGDP_FLL_ch
order country tot_exp_FLL_ch p_FLL_ch imr_FLL_ch rGDP_FLL_ch
*dataout, save(border_indexes_tariff_a) word replace dec(2)
*dataout, save(border_indexes_tariff_b) word replace dec(2)
*dataout, save(border_indexes_tariff_c) word replace dec(2)
dataout, save(border_indexes_tariff_d) word replace dec(2)
