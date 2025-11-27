**********************************************************************
*Estimatig Gravity: Some Applications and Practice/Homework Quastions*
**********************************************************************
*****************************************************************************
*NOTE: The following is a standard preamble of a do file, which sets various*
*dimensions, including the working directory and the log file, which will   *
*record everything that you do. The "text" option is useful for sharing.    *
*****************************************************************************
clear all
set more off
clear matrix
set memory 500m
set matsize 8000
set maxvar 30000
cd "/home/bepp/hchulkim/Documents/gravity/Gravity_data"
cap log close
log using gravity_applications, text replace

******************************************************************************************
*A. APPLICATIONS. The following are 4 questions that gave me the opportunity to introduce* 
*you to various approaches and practical tips for Stata coding, including some nice, new *
*and flexible Stata commands, and some important approaches and aspects for estimating   * *gravity, including the recommendations covered in class. Since we will not have time to *
*go over the file in class, I have added detailed notes, which should help. I will find*
*time (e.g., office hours) to address any questions that you may have on the applications*
******************************************************************************************
*************************************************
*I. APPLICATION 1: Traditional Gravity Estimates*
*************************************************
***************************************************************************
*NOTE: This is a panel dataset for aggregate manufacturing, which includes*
*domestic trade, for 69 countries over the period 1986-2006.              *
***************************************************************************
use gravity_data, clear
describe
summarize, detail
*********************
*1. Prepare the data*
*********************
**************************
*a. Prepare interval data*
**************************
keep if year == 1986 | year == 1990 | year == 1994 | year == 1998 | year == 2002 | year == 2006 
*****************************************************************************************
*NOTE: Using interval data is consistent with the older gravity literature. You should  *
*use data for consequtive years, i.e., all the data. For now, we will stick to intervals*
*for speed. At the end, the whole data will be used.                                    *
*****************************************************************************************
**************************************************************
*b. Create and transform some variables that are needed later*
**************************************************************
********************************************************
*i. Transform some varibales needed for the estimations*
********************************************************
********************************************************
*NOTE: Create the log of trade for the OLS regressions.*
********************************************************
gen ln_trade=ln(trade) 
gen ln_DIST=ln(DIST) 
label variable ln_DIST "Log of distance"
**********************************************************************
*NOTE: Create total output and expenditure. This should remind you of*
*the market clearing condition. Note that we have domestic trade!    *
**********************************************************************
bysort exporter year: egen output=sum(trade)  
gen ln_OUTPUT=ln(output)
bysort importer year: egen expenditure=sum(trade)
gen ln_EXPEND=ln(expenditure)
****************************************************
*ii. Create pair ID and symmetric pair ID variables*
****************************************************
************************************************************
*NOTE: The following creates an asymmetric country-pair ID.*
************************************************************
egen pair_id=group(exporter importer) 
********************************************************************
*NOTE: The following lines create a symmetric country-pair ID. I am* 
*curious if you can come up with a more efficient code to do this? *
********************************************************************
gen first = cond(exporter < importer, exporter, importer) 
gen second = cond(exporter < importer, importer, exporter)
gen symm = first+"_"+second
egen symm_pair_id = group(symm)
*******************************
*iii. Save data for future use*
********************************
save applications_data_1, replace 
*************************************
*2. Obtain traditional OLS Estimates*
*************************************
*********************
*a. Estimate Gravity*
*********************
reg ln_trade ln_DIST CNTG LANG CLNY ln_OUTPUT ln_EXPEND if exporter!=importer, cluster(pair_id) 
****************************************************
*b. Store the so they can be used for a nice table.*
****************************************************
estimates store ols
*************************************************
*c. Test if the output elasticity is really one.*
*************************************************
test ln_OUTPUT=1 
*************************************************************
*d. Perform the RESET/Ramsey test for model misspecification*
*************************************************************
************************************************************************
*NOTE: This is the implementation from Santos-Silva and Tenreyro (2006)*
************************************************************************
predict fit, xb
gen fit2=fit^2
reg ln_trade ln_DIST CNTG LANG CLNY ln_OUTPUT ln_EXPEND fit2 if exporter!=importer, cluster(pair_id)
test fit2=0
drop fit*
******************************************************************************
*NOTE: OLS fails the Ramsey test, just as in Santos-Silva and Tenreyro (2006)*
****************************************************************************** 
**********************************************
*3. Use Remoteness Indexes to control for MRs*
**********************************************
******************************
*a. Create Remoteness indexes*
******************************
*******************************************************************************
*NOTE: These are size-weighted distances on the importer and the exporter side*
*******************************************************************************
use applications_data_1, clear
collapse DIST [aw=expenditure], by(exporter year)
gen REM_EXP=ln(DIST)
save rmtns_exp, replace
use applications_data_1, clear
collapse DIST [aw=output], by(importer year)
gen REM_IMP=ln(DIST)
save rmtns_imp, replace
**************************************************************
*b. Combine data and estimate gravity with remoteness indexes*
**************************************************************
use applications_data_1, clear
joinby exporter year using rmtns_exp
joinby importer year using rmtns_imp
reg ln_trade ln_DIST CNTG LANG CLNY ln_OUTPUT ln_EXPEND REM_EXP REM_IMP if exporter!=importer, cluster(pair_id) 
estimates store rmtns
***********************
*c. Perform RESET Test*
***********************
predict fit, xb
gen fit2=fit^2
reg ln_trade ln_DIST CNTG LANG CLNY ln_OUTPUT ln_EXPEND REM_EXP REM_IMP fit2 if exporter!=importer, cluster(pair_id) 
test fit2=0
drop fit*
*****************************************
*4. Use Fixed Effects to control for MRs*
*****************************************
*************************
*a. Create Fixed Effects*
*************************
*************************************************************************************
*NOTE: The fixed effects are created explicitly, but later we will just use the IDs.*
*************************************************************************************
egen exp_time=group(exporter year)
tab exp_time, gen(exp_time_fe)
egen imp_time=group(importer year)
qui tab imp_time, gen(imp_time_fe)
*************************************************************************
*b. Estimate gravity with exporter-time and importer-time fixed effects.*
*************************************************************************
*************************************************************************************
*NOTE: We can no longer include any variable that varies acorss countries (and time)*
*************************************************************************************
reg ln_trade ln_DIST CNTG LANG CLNY exp_time_fe* imp_time_fe* if exporter!=importer, cluster(pair_id) 
******************************************************************************
*NOTE: Look at the dropped fixed effects. In cross-section, when the model is* 
*estimated with a constant, two of them will be dropped (one exporter and one* 
*importer). How many will be dropped in a panel setting?                     *
******************************************************************************
**********************************************************************************
*c. Estimate OLS gravity with the fast "reghdfe" command, and store the estimates*
**********************************************************************************
reghdfe ln_trade ln_DIST CNTG LANG CLNY if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id)
estimates store fes
******************************************************************************
*NOTE: Compare the speed of convergence between reg and reghdfe. Pretty cool.*
*It gets much better for non-linear estimations, e.g., PPML. Will show below.*
******************************************************************************
***********************
*c. Perform RESET Test*
***********************
predict fit, xb
gen fit2=fit^2
reghdfe ln_trade ln_DIST CNTG LANG CLNY fit2 if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id)
test fit2=0
drop fit*
***************************
*5. Use the PPML estimator*
***************************
*************************************
*a. Estimate Gravity & Store Results*
*************************************
***************************************************************************************
*NOTE: "ppml" is the original command from SST (2006). We will see alternatives below.*
***************************************************************************************
ppml trade ln_DIST CNTG LANG CLNY exp_time_fe* imp_time_fe* if exporter!=importer, cluster(pair_id)
estimates store ppml
********************
*b. Construct an R2*
********************
*************************************************************************
*NOTE: Since PPML is a pseudo-maximum likelihood estimator, the R^2 that*
*it reports is not valid. The simple code below delivers proper R^2.    *
*************************************************************************
predict trade_hat
qui cor trade_hat trade 
display "R-squared 0" (`r(rho)')^2 
drop trade_hat
***********************
*c. Perform RESET Test*
***********************
predict fit, xb
gen fit2=fit^2
ppml trade ln_DIST CNTG LANG CLNY fit2 exp_time_fe* imp_time_fe* if exporter!=importer, cluster(pair_id)
test fit2=0
drop fit*
***********************************************************
*NOTE: PPML passes the Ramsey test, just as in SST (2006).*
***********************************************************
********************************
*d. Alternative PPML estimators*
********************************
***********************************************************************
*NOTE: I have commented some of thos out, but good to be aware of them*
***********************************************************************
*glm trade ln_DIST CNTG LANG CLNY exp_time_fe* imp_time_fe* if exporter!=importer, cluster(pair_id) family(poisson) diff iter(30)
*ppml_panel_sg trade ln_DIST CNTG LANG CLNY if exporter!=importer, ex(exporter) im(importer) y(year) nopair cluster(pair_id)
********************************************************
*NOTE: I recomment using "ppmlhdfe". It is really nice!*
********************************************************
ppmlhdfe trade ln_DIST CNTG LANG CLNY if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id) 
********************************************************
*NOTE: You can store the estimates of the fixed effects*
********************************************************
ppmlhdfe trade ln_DIST CNTG LANG CLNY if exporter!=importer, absorb(exp_time imp_time, savefe) cluster(pair_id) 
********************************************************
*NOTE: You can store the estimates of the fixed effects* 
*as different variables, which will not be replaced.   *
********************************************************
ppmlhdfe trade ln_DIST CNTG LANG CLNY if exporter!=importer, absorb(expfe=exp_time impfe=imp_time) cluster(pair_id) 
********************************************
*NOTE: You can use alterantive clusterings.*
********************************************
ppmlhdfe trade ln_DIST CNTG LANG CLNY if exporter!=importer, absorb(exp_time imp_time, savefe) cluster(exporter importer year)
**************************************************
*6. Export estimates and save data for future use*
**************************************************
*********************************************************************
*NOTE: May use extension '.rtf' for word and '.csv' for excel tables*
*********************************************************************
esttab ols rmtns fes ppml using gravity_applications.rtf, append title(Traditional Gravity Estimates) mtitles(OLS RMTNS FES PPML) b(3) se(3) scalars(N r2) star(+ 0.10 * .05 ** .01) drop(exp_time_fe* imp_time_fe*) staraux nogaps 
save applications_data_2, replace

*************************************************************
*II. APPLICATION 2: A Simple Solution to the Distance Puzzle*
*************************************************************
**************************************************************
*1. Call data from previous application (applications_data_2)*
**************************************************************
use applications_data_2, clear
*****************************************************************
*2. Create time-varying distance variables (i.e., for each year)*
*****************************************************************
******************************************************************
*NOTE: We are looping over 4 years due to the intervals, however,*
*when all the data are used, then the loop should reflect it.    *
******************************************************************
forvalues i=1986(4)2006{
gen ln_DIST_`i'=ln_DIST if year==`i'
replace ln_DIST_`i'=0 if ln_DIST_`i'==.
}
******************************************************************************
*3. Estimate gravity with OLS, exp_time & importer-time FES, the time-varying* 
*distance variables and the standard gravity covariates CNTG LANG CLNY. Use  *
*ONLY international trade flows, as has been standard in the literature.     *
******************************************************************************
reghdfe ln_trade ln_DIST_1986 ln_DIST_1990 ln_DIST_1994 ln_DIST_1998 ln_DIST_2002 ln_DIST_2006 CNTG LANG CLNY if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id) 
estimates store dist_puzzl_ols
*****************************************************************
*NOTE: Seems like the effects of distance have not changed much.*
*****************************************************************
******************************************************************
*4. Test whether the impact of distance has changed significantly* 
*between 1986 and 2006. Hint: Use lincom/nlcom commands.         *
******************************************************************
***************************************************************
*NOTE: "nlcom" and "lincom" are nice, as they delives standard* 
*errors too, which are based on the Delta method.             *
***************************************************************
nlcom 100*(_b[ln_DIST_2006]-_b[ln_DIST_1986])/_b[ln_DIST_1986]
lincom _b[ln_DIST_2006]-_b[ln_DIST_1986]
********************************************************************
*5. Reproduce the results from the previous specification with PPML*
********************************************************************
ppmlhdfe trade ln_DIST_1986 ln_DIST_1990 ln_DIST_1994 ln_DIST_1998 ln_DIST_2002 ln_DIST_2006 CNTG LANG CLNY if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id) 
estimates store dist_puzzl_ppml
nlcom 100*(_b[ln_DIST_2006]-_b[ln_DIST_1986])/_b[ln_DIST_1986]
lincom _b[ln_DIST_2006]-_b[ln_DIST_1986]
***************************************
*NOTE: PPML does not solve the puzzle.*
***************************************
*****************************************************************************
*6. Reproduce the previous specification but with internal trade flows added* 
*to the sample. Use internal distance to proxy for internal trade costs.    * 
*****************************************************************************
******************************
*a. Generate/modify variables*
******************************
*********************************************************************************
*NOTE: If domestic trade flows are used then, at a minimum, you should introduce* 
*a border/same country dummy for domestic vs. international trade. Of course, by*
*construction, the estimates on SMCTRY and BRDR will be the same in magnitude & *
*opposite in sign.                                                              *
*********************************************************************************
gen SMCTRY=1 if exporter==importer
replace SMCTRY=0 if SMCTRY==.
gen BRDR=1 if exporter!=importer
replace BRDR=0 if BRDR==.
***********************************************************************************
*NOTE: The following code created time-varying variables for domestic distance too*
***********************************************************************************
gen ln_DIST_INTRA=ln_DIST*SMCTRY
forvalues i=1986(4)2006{
replace ln_DIST_`i'=0 if SMCTRY==1
}
******************************
*b. Estimate distance effects*
******************************
*******************************************************************************
*NOTE: Thus far, we control/proxy for domestic trade costs only with distance.*
*******************************************************************************
ppmlhdfe trade ln_DIST_1986 ln_DIST_1990 ln_DIST_1994 ln_DIST_1998 ln_DIST_2002 ln_DIST_2006 CNTG LANG CLNY ln_DIST_INTRA, absorb(exp_time imp_time) cluster(pair_id) 
estimates store dist_puzzl1
nlcom 100*(_b[ln_DIST_2006]-_b[ln_DIST_1986])/_b[ln_DIST_1986]
lincom _b[ln_DIST_2006]-_b[ln_DIST_1986]
*****************************
*NOTE: The puzzle is solved.*
*****************************
*****************************************************************************
*7. Reproduce previsous specification but, in addition to internal distance,*
*control for Home bias with a dummy variable SMTRY=0 for internal trade.    *
*****************************************************************************
ppmlhdfe trade ln_DIST_1986 ln_DIST_1990 ln_DIST_1994 ln_DIST_1998 ln_DIST_2002 ln_DIST_2006 CNTG LANG CLNY ln_DIST_INTRA SMCTRY, absorb(exp_time imp_time) cluster(pair_id) 
estimates store dist_puzzl2
nlcom 100*(_b[ln_DIST_2006]-_b[ln_DIST_1986])/_b[ln_DIST_1986]
lincom _b[ln_DIST_2006]-_b[ln_DIST_1986]
**********************************************************************************
*NOTE: Puzzle remains solved, but note the estimates on SMCTRY and ln_DIST_INTRA.*
**********************************************************************************
*********************************************************************
*8. Now use Country-specific Fixed Effects to capture time-invariant* 
*internal trade costs comprehensively.                              *
*********************************************************************
*********************************************************
*a. Create country-specific SMCTRY variables, which will* 
*capture all time-invbariant domestic trade costs.      *
*********************************************************
encode exporter, gen(exp_id)
gen SMCTRY_CNTRY=exp_id*SMCTRY
****************************************************************************************
*NOTE: We can create and add the fixed effects explicitly, but better to "absorb" them.*
****************************************************************************************
qui tab SMCTRY_CNTRY, gen(SMCTRY_CNTRY_FE)
ppmlhdfe trade ln_DIST_1986 ln_DIST_1990 ln_DIST_1994 ln_DIST_1998 ln_DIST_2002 ln_DIST_2006 CNTG LANG CLNY, absorb(exp_time imp_time SMCTRY_CNTRY) cluster(pair_id) 
estimates store dist_puzzl3
nlcom 100*(_b[ln_DIST_2006]-_b[ln_DIST_1986])/_b[ln_DIST_1986]
lincom _b[ln_DIST_2006]-_b[ln_DIST_1986]
**************************************************
*9. Export estimates and save data for future use*
**************************************************
esttab dist_puzzl_ols dist_puzzl_ppml dist_puzzl1 dist_puzzl2 dist_puzzl3 using gravity_applications.tex, append title(A Simple Solution to the Distance Puzzle In Trade) mtitles(OLS PPML INTRA ENDG LEAD PHSNG) b(3) se(3) scalars(N r2) star(+ 0.10 * .05 ** .01) drop() staraux nogaps 
save applications_data_3, replace

****************************************************
*III. APPLICATION 3: Estimating the Effects of RTAs*
****************************************************
****************************************************************************************
*1. Start with the data from the previous quation and add RTA data (rta_membership.dta)*
****************************************************************************************
use applications_data_3, clear
joinby exporter importer year using rta_membership, unmatched(master)
tab _merge
drop _merge
**************************************************************************
*2. Estimate gravity with OLS, exp_time & importer-time FES, the standard* 
*gravity covariates ln_DIST CNTG LANG CLNY RTA. Use ONLY international   *
*trade flows, as has been standard in the literature. Use reghdfe command*
**************************************************************************
reghdfe ln_trade ln_DIST CNTG LANG CLNY RTA if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id)
estimates store ols
****************************************************************************************
*NOTE: most gravity estimates look Ok, but the key RTA estimate is zero, economically &*
*statistically. This is a problem! Can you solve it? (PPML, pair FEs, domestic flows)? *
****************************************************************************************
********************************************************************
*3. Reproduce the previous results with PPML. Use ppmlhdfe command.*
********************************************************************
ppmlhdfe trade ln_DIST CNTG LANG CLNY RTA if exporter!=importer, absorb(exp_time imp_time) cluster(pair_id) 
*****************************************************************
*NOTE: The following generates the trade-volume effects of RTAs.*
*****************************************************************
nlcom (exp(_b[RTA])-1)*100
estimates store ppml
**********************************************************************************
*NOTE: PPML "helps" -- the RTA estimate is positive and significant. Still small.*
**********************************************************************************
*************************************************************************************
*4. Reproduce previous results with intra-national trade and a single SMCTRY control*
*************************************************************************************
ppmlhdfe trade ln_DIST CNTG LANG CLNY RTA SMCTRY, absorb(exp_time imp_time) cluster(pair_id) 
estimates store ppml_intra
****************************************************
*NOTE: No big difference, but we will revisit this.*
****************************************************
***********************************************************************************
*5. Produce RTA estimates with Avergare Treatment Effects, i.e. pair fixed effects*
*Hint: Experiment with symmetric and asymmetric pair fixed effects.               * 
***********************************************************************************
ppmlhdfe trade RTA, absorb(symm_pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_symm_pair
ppmlhdfe trade RTA, absorb(pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_asymm_pair
parmest, saving(pair_fes, replace)
*****************************************************************************************
*NOTE: It seems that the RTA effects with symmatric and asymmetric FEs are very similar.*
*This will change when we allow for asymmeties in the trade policy effects below. The   *
*safe bet is to always use asymmetric pair FEs.                                         *
*****************************************************************************************
***************************************************************************
*6. Reproduce the last two specifications but without domestic trade flows*
***************************************************************************
ppmlhdfe trade RTA if exporter!=importer, absorb(symm_pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_symm_pair_nodom
ppmlhdfe trade RTA if exporter!=importer, absorb(pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_asymm_pair_nodom
************************************
*NOTE: Domestic trade flows matter.*
************************************
************************************************
*7. Test for "strict exogeneity" with RTAs LEAD*
************************************************
ppmlhdfe trade RTA RTA_LEAD4, absorb(symm_pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_lead
********************************************************************************
*8. Allow for Phasing-in Effects of RTAs. Hint: Add RTA_LAG4 RTA_LAG8 RTA_LAG12*
********************************************************************************
ppmlhdfe trade RTA RTA_LAG4 RTA_LAG8 RTA_LAG12, absorb(symm_pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_phasing
******************************************************************************
*9. Construct the total RTA effect using command lincom and all RTA estimates*
******************************************************************************
lincom _b[RTA]+_b[RTA_LAG4]+_b[RTA_LAG8]+_b[RTA_LAG12]
********************************************************************************
*10. Reproduce the specification from part 8 but also control for globalization*
*using time-varying border variables for each year in the sample.              *
********************************************************************************
****************************
*a. Create Border Variables*
****************************
gen INTL_BRDR=1 if exporter!=importer
replace INTL_BRDR=0 if INTL_BRDR==.
forvalues i=1986(4)2006{
gen INTL_BRDR_`i'=1 if INTL_BRDR==1 & year==`i'
replace INTL_BRDR_`i'=0 if INTL_BRDR_`i'==.
}
*****************************************************
*b. Estimate gravity and construct total RTA effects*
*****************************************************
ppmlhdfe trade RTA RTA_LAG4 RTA_LAG8 RTA_LAG12 INTL_BRDR_1986 INTL_BRDR_1990 INTL_BRDR_1994 INTL_BRDR_1998 INTL_BRDR_2002 INTL_BRDR_2006, absorb(symm_pair_id exp_time imp_time) cluster(pair_id) 
estimates store ppml_glbzn
*******************************************************************************
*NOTE: One of the borders is dropped. Why? What will happen to the RTA and the*
*border estimates if you drop the border for 1986, i.e., INTL_BRDR_1986?      *
*******************************************************************************
******************************************************************************
*c. Construct the total RTA effect using command lincom and all RTA estimates*
******************************************************************************
lincom _b[RTA]+_b[RTA_LAG4]+_b[RTA_LAG8]+_b[RTA_LAG12]
********************************************************************
*11. Reproduce the last specification but, instead of including the* 
*globalization  dummies explicitly, add them to the absorb option. *
********************************************************************
egen brdr_time=group(INTL_BRDR year)
ppmlhdfe trade RTA RTA_LAG4 RTA_LAG8 RTA_LAG12, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
*************************************************************************************
*NOTE: The pracitcal benefit of this specification is that it is faster. Note that  *
*you can include any dummy variable in the absorb option, which will help with speed*
*************************************************************************************
****************************************************
*12. Save data for furture use and export estimates*
****************************************************
esttab ols ppml ppml_intra ppml_symm_pair ppml_asymm_pair ppml_symm_pair_nodom ppml_asymm_pair_nodom ppml_lead ppml_phasing ppml_glbzn using gravity_applications.tex, append title(Estimating The Effects of Regional Trade Agreements) mtitles(OLS PPML INTRA SYMM ASYYM LEAD PHSNG GLBZN) b(3) se(3) scalars(N r2) star(+ 0.10 * .05 ** .01) drop() staraux nogaps 
save applications_data_4, replace

********************************************************
*IV. APPLICATION 4: Gravity and RTAs with Staggered DiD*
********************************************************
***********************************************************************************
*NOTE: The staggered DiD methods are just being adopted in the gravity literature.*
*Thus, there are many opportunities for further improvements in this area.        *
*********************************************************************************** 
*******************************************************
*1. Start with the original dataset "gravity_data.dta"*
*******************************************************
use gravity_data, clear
**************************************************************
*2. Add the data on complete trade sanctions (sanct_data.dta)*
**************************************************************
joinby exporter importer year using rta_membership
**************************************************************
*3. Create IDs for fixed effects and any additional variables*
***************************************************************
egen exp_time=group(exporter year)
egen imp_time=group(importer year)
egen pair_id=group(exporter importer)
gen INTL_BRDR=1 if exporter!=importer
replace INTL_BRDR=0 if INTL_BRDR==.
egen brdr_time=group(INTL_BRDR year)
**************************************************
*3. Obtain "TWFE" estimates with standard command*
**************************************************
ppmlhdfe trade RTA, absorb(pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store twfe_ppml
*************************************************
*4. Obtain "TWFE" estimates with "jwdid" command*
*************************************************
xtset pair_id year, yearly
egen FT_RTA=csgvar(RTA), tvar(year) ivar(pair_id)
jwdid trade, ivar(pair_id) tvar(year) gvar(FT_RTA) method(ppmlhdfe) fevar(exp_time imp_time brdr_time)  hettype(twfe)
estimates store twfe_jwdid
**********************************************************************
*NOTE: As expected "ppmlhdfe" and "jwdid" with option "hettype(twfe)"* 
*deliver identical results, i.e., "jwdid" nests "ppmlhdfe".          *
**********************************************************************
****************************************************
*5. Obtain "ETWFE" estimates (with "jwdid" command)*
****************************************************
*************************
*a. Obtain the estimates*
*************************
jwdid trade, ivar(pair_id) tvar(year) gvar(FT_RTA) method(ppmlhdfe) fevar(exp_time imp_time brdr_time)
******************************
*b. Report the average effect*
******************************
estat simple, predict(xb) 
**************************************************************************************
*NOTE: Consistent with Nagengst and Yotov (AEJ: Applied, 2025), the ETWFE estimate of*
*the impact of RTAs is significantly larger than the corresponding TWFE estimate.    *
************************************************************************************** 
*******************************
*c. Create an event-type graph*
*******************************
estat event, predict(xb)
estat plot, ytitle("Treatment effect") xtitle("Years from treatment onset")
***************************
*Make graph a bit prettier*
***************************
estat plot, ytitle("Treatment effect") xtitle("Years from treatment onset") ///
xlabel(0/17, valuelabels angle(10) labsize(small)) graphregion(color(white))
**************************************************************************************
*NOTE: I added this so I can refer you to some graph options that I have found useful*
**************************************************************************************
graph export etwfe_post.pdf, replace
*****************************************************************
*NOTE: We see some interesting patterns at the end of the graph.*
*What is your explanation?                                      *
*****************************************************************
********************************************************
*d. Amend the analysis to allow for possible pre-trends*
********************************************************
jwdid trade, ivar(pair_id) tvar(year) gvar(FT_RTA) method(ppmlhdfe) fevar(exp_time imp_time brdr_time) never
estat event, predict(xb)
estat plot, ytitle("Treatment effect") xtitle("Years from treatment onset")
graph export etwfe_all.pdf, replace
******************************************************************
*NOTE: This graph, or part of it, can be used to study pre-trends*
******************************************************************

********************************************************************************
*B. PRACTICE/PROBLEM SET QUESTIONS. The following are a few questions that will*
*reninforce and extend on the previous applications.                           *
********************************************************************************
**************************************************
*Question 1: Evaluating the impact of WTO & NAFTA*
**************************************************
*******************************************************************************
*In this question, you are asked to evaluate the impact of NAFTA. Obtain the  *
*results from each specification below. Then combine them into a single table.* 
*Interpret each set of results.                                               *
*******************************************************************************   
***********************************************************************
*1. Start with the dataset that was just saved in Application 3, i.e.,* applications_data_4.dta, and reproduce the latest specification (11)  * 
*from application 3 but with a single RTA variable, i.e., do not allow*
*for phasing-in effects. (This is just for simplicity.)               *
***********************************************************************
use applications_data_4, clear
ppmlhdfe trade RTA, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store rta_wto1
*************************************************************
*2. Merge the "wto_data.dta" and add the WTO variable to the*
*previous specification. Then interpret the WTO estimates.  *                  
*************************************************************
joinby exporter importer year using wto_membership
ppmlhdfe trade RTA WTO, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store rta_wto2
**********************************************************************
*3. Now estimate the same model but without domestic trade flows and *
*compare the WTO estimates from this and the previous specifications.*
**********************************************************************
ppmlhdfe trade RTA WTO if exporter != importer, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store rta_wto3

* save data
esttab rta_wto1 rta_wto2 rta_wto3 ///
    using "/home/bepp/hchulkim/Documents/gravity/output/tables/tableQ2.tex", ///
    replace label nodepvar nogaps compress booktabs ///
    b(2) se(2) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(RTA WTO) ///
    varlabels(RTA "RTA member" ///
              WTO "WTO member") ///
    stats(N, labels("Observations"))
  

********************************************************************************
*4. Create a dummy variable for NAFTA and add it to the specification from part*
*(2). What is the impact of NAFTA? Interpret your results. Hint: Define NAFTA as*
*a dummy variable that takes a value of one for trade between Canada, USA, and *
*Mexico for the years after 1993, i.e., from 1994 onwards until the end of the * 
*sample. There are different ways to do this. I give you an example.           *
********************************************************************************
gen NAFTA=1 if exporter=="CAN" & importer=="USA" & year>=1994
replace NAFTA=1 if importer=="CAN" & exporter=="USA" & year>=1994
replace NAFTA=1 if exporter=="CAN" & importer=="MEX" & year>=1994
replace NAFTA=1 if importer=="CAN" & exporter=="MEX" & year>=1994
replace NAFTA=1 if exporter=="USA" & importer=="MEX" & year>=1994
replace NAFTA=1 if importer=="USA" & exporter=="MEX" & year>=1994
replace NAFTA=0 if NAFTA==.
*********************************************************************************
*5. You should've obtained and interpreted the estimates from the previous part *
*as deviations from the average FTA estimate. Now obtain the effects of NAFTA in*
*levels. Hint: subtract the NAFTA dummy from the RTA dummy and re-estimate.     *
*********************************************************************************
ppmlhdfe trade RTA WTO NAFTA, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store nafta1

gen RTA_ex=RTA-NAFTA
ppmlhdfe trade WTO NAFTA RTA_ex, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store nafta2
*************************************************************************
*4. Obtain estimates of the effects of NAFTA on trade between each pair,*
*e.g. MEX-CAN, CAN-USA, USA-MEX. Interpret/discuss your results.        * 
*************************************************************************
gen NAFTA_can_us=1 if exporter=="CAN" & importer=="USA" & year>=1994
replace NAFTA_can_us=1 if importer=="CAN" & exporter=="USA" & year>=1994
replace NAFTA_can_us=0 if NAFTA_can_us==.

gen NAFTA_can_mex=1 if exporter=="CAN" & importer=="MEX" & year>=1994
replace NAFTA_can_mex=1 if importer=="CAN" & exporter=="MEX" & year>=1994
replace NAFTA_can_mex=0 if NAFTA_can_mex==.

gen NAFTA_us_mex=1 if exporter=="USA" & importer=="MEX" & year>=1994
replace NAFTA_us_mex=1 if importer=="USA" & exporter=="MEX" & year>=1994
replace NAFTA_us_mex=0 if NAFTA_us_mex==.

ppmlhdfe trade WTO RTA_ex NAFTA_can_us NAFTA_can_mex NAFTA_us_mex, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store nafta3
*****************************************************************************
*5. Obtain estimates of the effects of NAFTA on trade between each pair,    *
*and also depending on the direction of trade flows (e.g., Mexico to Canada)*
*vs. Canada to Mexico. Interpret/discuss your results.                      *
*****************************************************************************
gen NAFTA_CAN_USA=1 if exporter=="CAN" & importer=="USA" & year>=1994
replace NAFTA_CAN_USA=0 if NAFTA_CAN_USA==.

gen NAFTA_USA_CAN=1 if exporter=="USA" & importer=="CAN" & year>=1994
replace NAFTA_USA_CAN=0 if NAFTA_USA_CAN==.

gen NAFTA_CAN_MEX=1 if exporter=="CAN" & importer=="MEX" & year>=1994
replace NAFTA_CAN_MEX=0 if NAFTA_CAN_MEX==.

gen NAFTA_MEX_CAN=1 if exporter=="MEX" & importer=="CAN" & year>=1994
replace NAFTA_MEX_CAN=0 if NAFTA_MEX_CAN==.

gen NAFTA_USA_MEX=1 if exporter=="USA" & importer=="MEX" & year>=1994
replace NAFTA_USA_MEX=0 if NAFTA_USA_MEX==.

gen NAFTA_MEX_USA=1 if exporter=="MEX" & importer=="USA" & year>=1994
replace NAFTA_MEX_USA=0 if NAFTA_MEX_USA==.

ppmlhdfe trade WTO RTA_ex NAFTA_CAN_USA NAFTA_USA_CAN NAFTA_CAN_MEX NAFTA_MEX_CAN NAFTA_USA_MEX NAFTA_MEX_USA, absorb(symm_pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store nafta4
************************************************************************************
*6. Instead of symmetric pair fixed effects, now use asymmetric pair fixed effects.*
*Did your estimates change? How? Why? Discuss.                                     *
************************************************************************************
ppmlhdfe trade WTO RTA_ex NAFTA_CAN_USA NAFTA_USA_CAN NAFTA_CAN_MEX NAFTA_MEX_CAN NAFTA_USA_MEX NAFTA_MEX_USA, absorb(pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store nafta5

* save data
esttab nafta1 nafta2 nafta3 nafta4 nafta5 ///
    using "/home/bepp/hchulkim/Documents/gravity/output/tables/tableQ2_2.tex", ///
    replace label nodepvar nogaps compress booktabs ///
    b(2) se(2) star(* 0.10 ** 0.05 *** 0.01) ///
    keep(RTA WTO RTA_ex NAFTA_can_us NAFTA_can_mex NAFTA_us_mex NAFTA_CAN_USA NAFTA_USA_CAN NAFTA_CAN_MEX NAFTA_MEX_CAN NAFTA_USA_MEX NAFTA_MEX_USA) ///
    varlabels(RTA "RTA member" ///
              WTO "WTO member" ///
	      RTA_ex "RTA exclude NAFTA" ///
	      NAFTA_can_us "NAFTA (CAN-USA)" ///
	      NAFTA_can_mex "NAFTA (CAN-MEX)" ///
	      NAFTA_us_mex "NAFTA (USA-MEX)" ///
	      NAFTA_CAN_USA "NAFTA (CAN-USA)" ///
	      NAFTA_USA_CAN "NAFTA (USA-CAN)" ///
	      NAFTA_CAN_MEX "NAFTA (CAN-MEX)" ///
	      NAFTA_MEX_CAN "NAFTA (MEX-CAN)" ///
	      NAFTA_USA_MEX "NAFTA (USA-MEX)" ///
	      NAFTA_MEX_USA "NAFTA (MEX-USA)") ///
    stats(N, labels("Observations"))
  
*******************************************************************
*Question 2: Treatments with entry and exit: The case of sanctions*
*******************************************************************
*******************************************************
*1. Start with the original dataset "gravity_data.dta"*
*******************************************************
use gravity_data, clear
**************************************************************
*2. Add the data on complete trade sanctions (sanct_data.dta)*
**************************************************************
joinby exporter importer year using Sanctions_data
**************************************************************
*3. Create IDs for fixed effects and any additional variables*
***************************************************************
egen exp_time=group(exporter year)
egen imp_time=group(importer year)
egen pair_id=group(exporter importer)
gen INTL_BRDR=1 if exporter!=importer
replace INTL_BRDR=0 if INTL_BRDR==.
egen brdr_time=group(INTL_BRDR year)
*************************************************************************
*3. Obtain a single "ETWFE" estimate of the impact of complete sanctions* 
*on trade from  a PPML model with three-way fixed effects.              *
*************************************************************************
ppmlhdfe trade SANCT, absorb(pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store sanct1

xtset pair_id year, yearly
egen FT_SANC=csgvar(SANCTION), tvar(year) ivar(pair_id)
jwdid trade, ivar(pair_id) tvar(year) gvar(FT_SANC) method(ppmlhdfe) fevar(exp_time imp_time brdr_time)
estat simple, predict(xb) 
* estimates
*-0.816 (0.114), p=0.000
***************************************************************************
*4. Extend the previous specification to allow for the impact of sanctions*
*to vary before their imposition, during their imposition, and after their*
*lifting, i.e., obtain "event-type" estimates. Allow for 5 leads, 5 pahsing
*in effects, and 5 effects after the sanctions are lifted.                *
***************************************************************************
**************************
*a. Create leads and lags*
**************************
bysort pair_id (year): egen start_s = min(cond(SANCT==1, year, .))
bysort pair_id (year): egen end_s   = max(cond(SANCT==1, year, .))

gen rel_start = year - start_s if start_s < .
gen rel_end   = year - end_s   if end_s   < .

forvalues k = 1/5 {
    gen lead_m`k' = rel_start == -`k'
}

forvalues k = 0/4 {
    gen phase`k' = rel_start == `k' & SANCT == 1
}

forvalues k = 1/5 {
    gen post`k' = rel_end == `k' & SANCT == 0
}


*********************
*b. Obtain estimates*
*********************
ppmlhdfe trade lead_m1 lead_m2 lead_m3 lead_m4 lead_m5 phase0 phase1 phase2 phase3 phase4 post1 post2 post3 post4 post5, absorb(pair_id exp_time imp_time brdr_time) cluster(pair_id) 
estimates store sanct2

*******************
*c. Create a graph*
*******************
estimates restore sanct2
parmest, norestore

* Keep only the event-study dummies
keep if parm=="lead_m1" | parm=="lead_m2" | parm=="lead_m3" | parm=="lead_m4" | parm=="lead_m5" ///
    | parm=="phase0" | parm=="phase1" | parm=="phase2" | parm=="phase3" | parm=="phase4" ///
    | parm=="post1" | parm=="post2" | parm=="post3" | parm=="post4" | parm=="post5"


* Relative year (x-axis): you can tweak this mapping if your coding differs
gen rel_year = .

* Leads: -1, -2, ..., -5
replace rel_year = -1 if parm=="lead_m1"
replace rel_year = -2 if parm=="lead_m2"
replace rel_year = -3 if parm=="lead_m3"
replace rel_year = -4 if parm=="lead_m4"
replace rel_year = -5 if parm=="lead_m5"

* Phase-in: 0,1,2,3,4
replace rel_year = 0 if parm=="phase0"
replace rel_year = 1 if parm=="phase1"
replace rel_year = 2 if parm=="phase2"
replace rel_year = 3 if parm=="phase3"
replace rel_year = 4 if parm=="phase4"

* Listed / post: 5,6,7,8,9  (you can choose other spacing if you like)
replace rel_year = 5 if parm=="post1"
replace rel_year = 6 if parm=="post2"
replace rel_year = 7 if parm=="post3"
replace rel_year = 8 if parm=="post4"
replace rel_year = 9 if parm=="post5"

* Group indicator for colors
gen group = .
replace group = 1 if strpos(parm,"lead_m")
replace group = 2 if strpos(parm,"phase")
replace group = 3 if strpos(parm,"post")

label define grouplab 1 "Leads" 2 "Phase-in" 3 "Listed"
label values group grouplab

twoway ///
    /// confidence intervals for all points
    (rcap min95 max95 rel_year, sort) ///
    /// leads
    (scatter estimate rel_year if group==1, ///
        mcolor(navy) msymbol(O)) ///
    /// phase-in
    (scatter estimate rel_year if group==2, ///
        mcolor(orange) msymbol(D)) ///
    /// listed / post
    (scatter estimate rel_year if group==3, ///
        mcolor(forest_green) msymbol(S)), ///
    ///
    xline(-1, lpattern(dash) lcolor(gs8)) ///
    yline(0, lpattern(dash) lcolor(gs8)) ///
    xtitle("Years relative to sanctions onset/lifting") ///
    ytitle("Treatment effect (PPML estimate)") ///
    legend(order(2 "Leads" 3 "Phase-in" 4 "Lifted") ///
           pos(6) ring(0) cols(1)) ///
    xlabel(-5(1)9)
graph save "sanctions_eventstudy.gph", replace
graph export "/home/bepp/hchulkim/Documents/gravity/output/figures/sanctions_eventstudy.png", replace width(2000)


************************************
*d. Report and discuss your results*
************************************
** In the report.


