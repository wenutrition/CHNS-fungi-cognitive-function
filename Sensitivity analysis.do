
*----------------------------------------NA revise--------------------------------------------

*---------------Species level analysis------------

clear
clear matrix
clear mata
set max_memory 32g
set maxvar 20000

import excel "D:\Zhenglab  data\Microbiota\fungi\CHNS_phase4\fungi\species.xlsx",sheet("Sheet1") firstrow clear 
save  "D:\CNHS\data_pre\CHNS18_fungi_species.dta",replace

use "D:\CNHS\data\data_pre\phenotype18_raw",clear
duplicates drop SampleID, force
merge 1:1 SampleID using "D:\CNHS\data_pre\CHNS18_fungi_species.dta"
keep if _merge==3
drop if MMSE_score==.

foreach v of varlist   s378 s379 s165-s180 MMSE_score insomnia_score HbA1C tg chol LDL_CH HDL_CH glu{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
sort region
egen code_id=group(region)

tempname coef
tempfile res
postfile `coef' str200(outcome exposure ) float(n rr lul uul p)  using "`res'", replace

foreach outcome of varlist MMSE_score insomnia_score HbA1C tg chol LDL_CH HDL_CH glu{ 
foreach var of varlist s378 s379 s165-s180 { 
    mixed `outcome' `var' age i.sex BMI   ||code_id:, covariance(unstructured) 
	test `var'
	post `coef'   ("`outcome'") ("`var'") (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
    }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final\statistic_update.xlsx",  firstrow(variables) sheet("species_mmse") sheetreplace  
restore 



*------Furter adjustments---------

use "D:\真菌分析\data_pre\mmse_fungi_GNHS_raw.dta",clear
regress mmse g80 age i.sex BMI 
regress mmse g80 age i.sex BMI smoke alc heart_disease ganzang_disease weichang_disease huxi_disease miniao_disease neifenmi_disease shenjing_disease guge_disease
 
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
    regress `outcome' `var' age i.sex BMI smoke alc heart_disease ganzang_disease weichang_disease huxi_disease miniao_disease neifenmi_disease shenjing_disease guge_disease
 
test `var'
 post `coef'   ("`outcome'") ("`var'")  (e(N)) (_b[`var'])  (_b[`var']-1.96*_se[`var'])  (_b[`var']+1.96*_se[`var'])  (chi2tail(1,(_b[`var']/_se[`var'])^2))  
 }
 }

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final\statistic.xlsx",  firstrow(variables) sheet("GNHS_fungi_mmse_furtheradjust") sheetreplace  
restore

*-------MIND diet and components intra associations--------

use "D:\CNHS\data\data_pre\phenotype18_raw",clear
drop if SampleID==""
duplicates drop SampleID,force 
merge 1:1 SampleID using  "D:\真菌分析\data_pre\18allfungi_glevel.dta" 
keep if _merge==3 
drop rice-others _merge diet_index
rename Idind idind
merge 1:1 idind using  "D:\真菌分析\data_pre\CHNS_2018ffq.dta" 
keep if _merge==3 //
su age 
su diet_index 
foreach v of varlist age diet_index whole_grain fruit greenveg otherveg nut tea bean fish pourtry  sweet redmeat fryfood{
           egen s_`v' = std(`v')
		   drop `v'
		   rename s_`v' `v'
         }
mixed diet_index whole_grain fruit greenveg otherveg nut tea bean fish pourtry  sweet redmeat fryfood  age i.sex BMI   ||region:, covariance(unstructured) 



*----------------keep the age>40 participants for analysis----

use  "D:\真菌分析\data_pre\18_fung_statistic.dta" ,clear
keep if age>40
mixed  g280  age i.sex BMI ||region:, covariance(unstructured) 

egen std_age=std(age)
drop age
rename std_age age

tempname coef
tempfile res
postfile `coef' str200( micro ) float(n rr lul uul p)  using "`res'", replace

foreach var of varlist shannon_entropy-g761 { 
    mixed  `var'  age i.sex BMI ||region:, covariance(unstructured) 
test age
 post `coef'   ("`var'") (e(N)) (_b[age])  (_b[age]-1.96*_se[age])  (_b[age]+1.96*_se[age])  (chi2tail(1,(_b[age]/_se[age])^2))  
 }
 

  postclose `coef'
  preserve
use "`res'", clear
export excel using "D:\真菌aging\results_paper_final\statistic.xlsx",  firstrow(variables) sheet("CHNS_funage_clr_overall_age40") sheetreplace  
restore


import excel "D:\真菌aging\results_paper_final\statistic.xlsx", sheet("CHNS_funage_clr_overall_age40") firstrow clear 
rename (n-p) =_overall
save  "D:\真菌分析\data_pre\18_fung_age_clr_overall40.dta" ,replace

use "D:\真菌分析\data_pre\18_fung_age_clr_overall40.dta",clear
merge 1:1 micro using  "D:\真菌分析\data_pre\18fungi_mapping.dta"
count if p_overall<0.05 // 29
drop _merge
save "D:\真菌aging\plot data\age_fungi_clr40.dta",replace

*--comebine the age related fungi across difference cohorts----

use "D:\真菌aging\plot data\age_fungi_clr40.dta",clear
rename micro g_CHNS
merge m:m g_CHNS using "D:\真菌aging\plot data\age_fungi_GNHS.dta"
drop _merge
merge m:m g_CHNS using "D:\真菌aging\plot data\age_fungi_zousan.dta"
drop _merge
merge m:m g_CHNS using "D:\真菌aging\plot data\age_fungi_15.dta"
drop if n_overall==.
drop _merge
save  "D:\真菌aging\plot data\age_fungi_association.dta",replace
keep if g_CHNS=="g280" | g_CHNS=="g139" |g_CHNS=="g137" |g_CHNS=="g186"	|g_CHNS=="g332"	|g_CHNS=="g717"	|g_CHNS=="g702"	|g_CHNS=="g564"	|g_CHNS=="g700"	|g_CHNS=="g305"	|g_CHNS=="g149"	|g_CHNS=="g347"	|g_CHNS=="g78"										
count if p_GNHS_overall<0.05 & p_zousan_overall<0.05 
save  "D:\真菌aging\plot data\age_fungi_association_plotupdate.dta",replace

*---------------Sacc age associations update------------------

import excel "D:\真菌aging\results_paper_final_final\statistic.xlsx", sheet("sacc_ageupdate") firstrow clear

metan b CI_lower CI_upper, label(namevar=Cohorts) nooverall nowt nobox ///
    fixed effect(b) ///
    boxopt(mcolor(navy*0.8) msymbol(circle)) ///  // 修改为圆形
    pointopt(msymbol(circle) mcolor(navy*0.8) msize(medium)) /// // 统一设定圆形和固定大小
    ciopt(lcolor(navy*0.8) lwidth(medium)) ///
    null(0) textsize(120) astext(75) force ///
    graphregion(fcolor(white) lcolor(white)) ///
    xlab(-0.17, -0.1, -0.05, 0) subtitle("") // 建议加上0刻度以便观察

graph export "D:\真菌aging\results_paper_final_final\sacc_age_metaupdate.pdf", as(pdf) replace

*---------------MMSE prediction compair------------------

use "D:\真菌分析\data_pre\pros_predictdata1.dta",clear  
preserve
keep SampleID MMSE_score diet_index  
save "D:\真菌分析\data_pre\2015_mmsepredict_conti_diet.dta",replace  
restore

use "D:\真菌分析\data_pre\pros_predictdata1.dta",clear  
preserve
keep SampleID MMSE_score  age sex BMI
save "D:\真菌分析\data_pre\2015_mmsepredict_agesexbmi.dta",replace  
restore

preserve
keep SampleID MMSE_score  region  age sex BMI diet_index g280	g139	g137	g186	g332	g717	g702	g564	g700 g305	g149	g347	g78
save "D:\真菌分析\data_pre\2015_mmsepredict_conti_all.dta",replace  
restore

preserve
keep SampleID MMSE_score  age sex BMI diet_index
save "D:\真菌分析\data_pre\2015_mmsepredict_agesexbmidiet.dta",replace  
restore

*-----GNHS MMSE prediction----------

use "D:\真菌aging\final_results\GNHS_F3phe.dta",clear
gen ID=substr(SampleID,3,6)
merge 1:1 ID using "D:\CNHS\data_pre\GNHS_phe_fungi.dta"
keep if _merge==3 // 1147
keep SampleID mmse g1-g204
keep SampleID mmse g41	g43	g53	g77	g80	g95	g99	g142	g189	g191	g196
drop if mmse==.
rename mmse MMSE_score
save "D:\CNHS\data_pre\GNHS_fungimmsepredictiton.dta",replace

*--------复旦MMSE and fungi dis--------

import delimited "D:\真菌aging\复旦验证\Re_ Re_ 真菌与认知功能如皋长寿与衰老队列验证结果\File3 ITS_counts.csv",  clear 
rename mmse_score mmse
gen cohort="RLAS"
keep mmse cohort
append using  "D:\真菌分析\data_pre\mmse_compare.dta"
tab cohort
save "D:\真菌分析\data_pre\mmse_compare_update.dta",replace


*-----interaction analysis----

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
drop if diet_index==. | MMSE_score==.
merge 1:1 SampleID using  "D:\真菌分析\data_pre\2015_fungi_mmse_prosid.dta"
keep if _merge==3 
xtile cat_bac = g139, nq(2)
regress MMSE_score  c.diet_index##c.g139   age i.sex BMI   

sort cat_bac
tab cat_bac
by cat_bac:regress MMSE_score  diet_index   age i.sex BMI  

by cat_bac:su MMSE_score  diet_index  

