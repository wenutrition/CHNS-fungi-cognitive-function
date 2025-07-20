
*-------describe the variable of the CHNS cross-sectional cohort--

use "D:\CNHS\data\data_pre\phenotype18_raw",clear
duplicates drop SampleID,force
merge 1:1 SampleID using "D:\真菌分析\data_pre\2018_fungi_clr_trans.dta"
keep if _merge==3
drop rice-others vegs-diet_index _merge
rename Idind idind
drop if MMSE_score==.
save "D:\CNHS\data\data_pre\CHNS_MMSE",replace

merge 1:1 idind using "D:\真菌分析\data_pre\2018_ffqdata.dta"
keep if _merge==3 | _merge==1 // 7411 participants
table1, vars(age contn \ sex cat \ BMI contn\ diet_index contn\ whole_grain contn \ greenveg contn \otherveg contn \fruit contn \nut contn\fish contn \bean contn \pourtry contn \tea contn\ redmeat contn \fryfood contn \sweet contn \ d3_energy contn \ MMSE_score contn)format(%8.2f) onecol test pdp(2) saving ("D:\真菌aging\results_paper_final_final\Supplemental files\CHNS_cross_des_overall.xlsx", replace)

preserve
keep MMSE_score d3_energy diet_index whole_grain greenveg otherveg fruit nut fish bean pourtry tea redmeat fryfood sweet
order  MMSE_score d3_energy diet_index whole_grain greenveg otherveg fruit nut fish bean pourtry tea redmeat fryfood sweet
drop if MMSE_score==. | diet_index==.
su diet_index,de
su MMSE_score,de
save  "D:\真菌分析\data_pre\cross_data_des.dta",replace
restore

replace mmse=2 if mmse==0.5 | mmse==1
replace mmse=1 if mmse==0.25
tab mmse
table1, by(mmse) vars(age contn \ sex cat \ BMI contn\ diet_index contn\ whole_grain contn \ greenveg contn \otherveg contn \fruit contn \nut contn\fish contn \bean contn \pourtry contn \tea contn\ redmeat contn \fryfood contn \sweet contn \ d3_energy contn \ MMSE_score contn)format(%8.2f) onecol test pdp(2) saving ("D:\真菌aging\results_paper_final_final\Supplemental files\CHNS_cross_des.xlsx", replace)

*-------describe the variable of the CHNS prospective cohort--

use "D:\真菌分析\data_pre\2015_ffqdata.dta",clear
destring idind,replace
rename idind Idind
merge 1:1 Idind using "D:\真菌分析\data_pre\2015_fungi_mmse_count_pros.dta" 
drop if global_mmse1<7 // 
drop if mmse==.
replace mmse=2 if mmse==0.5 | mmse==1
replace mmse=1 if mmse==0.25
tab mmse

table1,  vars(age contn \ sex cat \ BMI contn\ diet_index contn\ whole_grain contn \ greenveg contn \otherveg contn \fruit contn \nut contn\fish contn \bean contn \pourtry contn \tea contn\ redmeat contn \fryfood contn \sweet contn \ d3_energy contn \ MMSE_score contn)format(%8.2f) onecol test pdp(2) saving ("D:\真菌aging\results_paper_final_final\Supplemental files\CHNS_pros_overall_des.xlsx", replace)
table1, by(mmse) vars(age contn \ sex cat \ BMI contn\ diet_index contn\ whole_grain contn \ greenveg contn \otherveg contn \fruit contn \nut contn\fish contn \bean contn \pourtry contn \tea contn\ redmeat contn \fryfood contn \sweet contn \ d3_energy contn \ MMSE_score contn)format(%8.2f) onecol test pdp(2) saving ("D:\真菌aging\results_paper_final_final\Supplemental files\CHNS_pros_des.xlsx", replace)


*-------------------------Dataset for Permonova analyses------------------

*---Mind diet and components and MMSE with fungi--

use  "D:\CNHS\data\data_pre\18_phenotype_fungi.dta",clear // 11456
drop if MMSE_score==. // 7411 keep 
merge 1:1 idind using  "D:\真菌分析\data_pre\CHNS_2018ffq.dta" // 11456 
keep if _merge==3 // 6962
su age // 40-96

preserve
keep g2-g761
save   "D:\CNHS\data\data_pre\18micro.dta",replace
restore

preserve
keep age sex MMSE_score whole_grain greenveg otherveg fruit  nut tea bean fish pourtry sweet redmeat fryfood diet_index
save   "D:\CNHS\data\data_pre\18phenotype.dta",replace
restore

*----Fungi-mmse PCOA analysis

*--cross sectional--

use  "D:\CNHS\data\data_pre\FIanalysis_plot.dta",clear
preserve
replace mmse=1 if mmse==0.5
drop if mmse==.
keep  mmse shannon_entropy-pielou_evenness
save   "D:\CNHS\data\data_pre\age_mmse_diversity.dta",replace
restore
preserve
tab mmse,missing
drop if mmse==.
keep SampleID mmse g2-g761
save   "D:\CNHS\data\data_pre\age_pcoa_mmse.dta",replace // 7411
restore

*--prospective---

use "D:\真菌分析\data_pre\2015_fungi_mmse_count_pros.dta",clear // 3519 participants with fungi at 15, and mmse at 18
su global_mmse1,de // 5分为下1%，7分为5%，9分为下10%
count if global_mmse1<=5 // 47
count if global_mmse1<7 // 113
replace mmse=1 if mmse==0.5
drop if mmse==.
drop if global_mmse1<7
tab mmse,missing 


*------------------------MMSE related fungal taxa------------

*--CHNS-2018--

use "D:\CNHS\data\data_pre\phenotype18_raw",clear
duplicates drop SampleID,force
merge 1:1 SampleID using "D:\真菌分析\data_pre\2018_agepredict_fillter_clr.dta"
keep if _merge==3
keep if age>=18

foreach v of varlist age shannon_entropy-g761 FI MMSE_score pss10 disability disease1{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace

foreach outcome of varlist  FI MMSE_score pss10 disability disease1  {
foreach var of varlist  shannon_entropy-pielou_evenness g2-g761 { 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final\statistic.xlsx",  firstrow(variables) sheet("CHNS_funFI_CLR1") sheetreplace  
restore

*-----Mapping the results---

import excel "D:\真菌aging\results_paper_final\statistic.xlsx", sheet("CHNS_funFI_CLR1") firstrow clear 
rename micro outcome
rename var micro
save  "D:\真菌分析\data_pre\fungi_FI_18_clr.dta" ,replace

use  "D:\真菌分析\data_pre\fungi_FI_18_clr.dta",clear
merge m:m micro using  "D:\真菌分析\data_pre\18fungi_mapping.dta"
keep if _merge==1 | _merge==3
sort p
order outcome Genus p
keep if _merge==3
save  "D:\真菌分析\data_pre\18_fung_FI_association_clr.dta" ,replace

*--CHNS-2015 prospective analysis--

use "D:\真菌分析\data_pre\2015_fungi_mmse_count_pros.dta",clear // 3519 participants with fungi at 15, and mmse at 18
drop g1-g751 
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
keep if _merge==3
drop  if global_mmse1<=5 // 47
drop if MMSE_score==.
foreach v of varlist  g1-g751 MMSE_score{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace

foreach outcome of varlist   MMSE_score   {
foreach var of varlist   g1-g751 { 
    mixed `outcome' `var' age i.sex BMI ||region:, covariance(unstructured) 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final\statistic.xlsx",  firstrow(variables) sheet("CHNS15_mmse_pros") sheetreplace  
restore

*-----Mapping the results---

import excel "D:\真菌aging\results_paper_final\statistic.xlsx", sheet("CHNS15_mmse_pros") firstrow clear 
drop micro
rename var micro
merge m:m micro using  "D:\真菌分析\data_pre\18fungi_mapping.dta"
keep if _merge==3
sort p
order Family Genus rr p
rename rr rr_15
rename p p_15
rename lul lul_15
rename uul uul_15
drop _merge
save "D:\真菌分析\data_pre\15_fungi_mmse_prosig.dta",replace

*------GNHS prospective-----

use "D:\真菌aging\final_results\GNHS_F3phe.dta",clear
gen ID=substr(SampleID,3,6)
merge 1:1 ID using "D:\CNHS\data_pre\GNHS_phe_fungi.dta"
keep if _merge==3 // 1147
drop if 帕金森病==1 | 老年痴呆==1 // drop 9
spearman mmse g80
regress mmse g80 age i.sex i.edu3 BMI 

*--some participants with error

foreach v of varlist 空间定向力得分4 空间定向力得分5 复述能力得分 {
		   replace `v'=1 if `v'>1
         }
		 
egen mmse_=rowtotal(时间定向力得分1-结构能力得分)
regress mmse_ g80 age i.sex i.edu3 BMI 

drop outcome
rename mmse_ outcome
 foreach micro of varlist g1-g204{
	count if `micro'==0 
	if r(N)>1118{	
	drop  `micro' 
	}
}

foreach v of varlist age pielou_e_its-g203 outcome{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
tempname coef
tempfile res
postfile `coef' str200( micro var) float(n rr lul uul p)   using "`res'", replace

foreach outcome of varlist  outcome  {
foreach var of varlist  pielou_e_its-g196 { 
    regress `outcome' `var' age i.sex BMI 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final\statistic.xlsx",  firstrow(variables) sheet("GNHS_fungi_mmse") sheetreplace  
restore

*---mapping the results---

import excel "D:\真菌aging\results_paper_final\statistic.xlsx", sheet("GNHS_fungi_mmse") firstrow clear 
rename var g_GNHS
merge m:m g_GNHS using "D:\真菌分析\data_pre\GNHS_codemap.dta"
keep if _merge==3
sort p 
rename rr rr_GNHS
rename lul lul_GNHS
rename uul uul_GNHS
rename p p_GNHS
drop _merge
save "D:\真菌分析\data_pre\GNHS_fungi_mmsesig.dta",replace
	
*-------overlaped signals----

*----15 vs 18------

import delimited "D:\真菌aging\results_paper_final\fungi_mmse.csv",  clear 
keep if p_adjust<0.05 
merge 1:1 micro using "D:\真菌分析\data_pre\15_fungi_mmse_prosig.dta"
keep if _merge==3 
keep if p_15<0.05
order Family Genus p_15 p rr_15 rr 
 
*---combine the data for plot---

import delimited "D:\真菌aging\results_paper_final\fungi_mmse.csv",  clear 
keep if p_adjust<0.05 // 27
keep micro rr lul uul p_adjust genus
gen Cohort="Discovery"
rename p_adjust p
save "D:\真菌分析\data_pre\18_fungi_mmse_siginput.dta",replace

preserve
keep micro genus
merge 1:1 micro using "D:\真菌分析\data_pre\15_fungi_mmse_prosig.dta"
keep if _merge==1 | _merge==3 
rename rr_15 rr
rename lul_15 lul
rename uul_15 uul
rename p_15 p
gen Cohort="Validation1"
keep micro rr lul uul p genus Cohort
save "D:\真菌分析\data_pre\15_fungi_mmse_siginput.dta",replace
restore

preserve
keep micro genus
rename micro g_CHNS
merge m:m g_CHNS using "D:\真菌分析\data_pre\GNHS_fungi_mmsesig.dta"
keep if _merge==1 | _merge==3 
rename rr_GNHS rr
rename lul_GNHS lul
rename uul_GNHS uul
rename p_GNHS p
gen Cohort="Validation2"
keep micro rr lul uul p genus Cohort
save "D:\真菌分析\data_pre\GNHS_fungi_mmse_siginput.dta",replace
 

*---------------------------Fungi diet associations------------------

*--2018 cross-sectional analysis

use "D:\CNHS\data\data_pre\phenotype18_raw",clear
duplicates drop SampleID,force
merge 1:1 SampleID using "D:\真菌分析\data_pre\2018_fungi_clr_trans.dta"
keep if _merge==3
drop rice-others vegs-diet_index _merge
rename Idind idind
merge 1:1 idind using "D:\真菌分析\data_pre\2018_ffqdata.dta"
keep if _merge==3 // 8736 participants
drop _merge 
merge 1:1 idind using   "D:\CNHS\data\data_pre\18_fungidiversity.dta" 
keep if _merge==3 // 8736 participants
drop _merge
foreach v of varlist  shannon_entropy-pielou_evenness g4-g761 greenveg-fryfood diet_index{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
sort region
egen code_id=group(region)

tempname coef
tempfile res
postfile `coef' str200(outcome exposure ) float(n rr lul uul p)  using "`res'", replace

foreach outcome of varlist shannon_entropy-pielou_evenness g4-g761 { 
foreach var of varlist greenveg-fryfood diet_index{ 
    mixed `outcome' `var' age i.sex BMI d3_energy  ||code_id:, covariance(unstructured) 
	test `var'
	post `coef'   ("`outcome'") ("`var'") (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
    }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final_final\statistic.xlsx",  firstrow(variables) sheet("ffq_diet_fungi_2018") sheetreplace  
restore 

*---------Interaction analysis----------

use "D:\真菌分析\data_pre\2015_ffqdata.dta",clear
destring idind,replace
rename idind Idind
merge 1:1 Idind using  "D:\真菌分析\data_pre\2015_fungi_mmse_count_pros.dta"
keep if _merge==3
keep Idind SampleID MMSE_score global_mmse diet_index d3_energy age sex BMI 
rename  Idind idind
merge 1:1 idind using "D:\真菌分析\data_pre\2015_all_dataset_fill.dta"
keep if _merge==3
drop g1-g751 _merge
merge 1:1 SampleID using "D:\真菌分析\data_pre\fungi_clr.dta"
keep if _merge==3
drop _merge // 
drop if global_mmse<7
drop if diet_index==. | MMSE_score==. // 3300
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_fungi_mmse_prosid.dta"
keep if _merge==3 // 3205
xtile cat_bac = g139, nq(2)
regress MMSE_score  c.diet_index##c.g139   age i.sex BMI i.education  
sort cat_bac
tab cat_bac
by cat_bac:regress MMSE_score  diet_index   age i.sex BMI i.education 
