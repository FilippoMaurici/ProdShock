clear

cd C:\Users\filip\OneDrive\Desktop\Lavori\Miei\PROD_TFP

use "C:\Users\filip\OneDrive\Desktop\Lavori\Miei\PROD_TFP\pwt1001.dta" 

gen CS=1- labsh 
* capital share

gen TFP=ctfp *rtfpna

drop rtfpna

rename TFP rtfpna

encode country, generate(country_id)
xtset country_id year

sort country_id year

drop if CS == CS[_n-1]


reg CS rtfpna , robust
est sto e1

reg CS rtfpna i.year, robust
est sto e2

*//xtreg CS rtfpna, fe
**est sto e3

**xtreg CS rtfpna i.year, fe
**est sto e4

* outreg2 [e1 e2 ] using "al_TFP.tex", replace

betareg CS rtfpna , vce(robust)
est sto betareg1
predict CS_betareg
betareg CS rtfpna i.year, vce(robust)
est sto betareg2
predict CS_betareg_Y
*betareg CS rtfpna i.country_id, vce(robust)
*est sto betareg3
betareg CS rtfpna i.year i.country_id, vce(robust)
est sto betareg4
predict CS_betareg_YC

fracreg probit CS rtfpna , vce(robust)
est sto fracreg1
predict CS_fracreg
fracreg probit CS rtfpna i.year, vce(robust)
est sto fracreg2
predict CS_fracreg_Y
*fracreg probit CS rtfpna i.country_id, vce(robust)
*est sto fracreg3
fracreg probit CS rtfpna i.year i.country_id, vce(robust)
est sto fracreg4
predict CS_fracreg_YC

outreg2 [betareg1 betareg2 betareg4 fracreg1 fracreg2 fracreg4] using "al_TFP.tex", replace

save al_TFP.dta, replace

keep year country_id CS* rtfpna

export excel using "al_TFP_pred.xlsx", sheet("Sheet1") firstrow(variables) replace


clear 
use al_TFP.dta

* Create empty variables to store coefficients (b), upper bound (u), lower bound (d), and lag number
gen b = .
gen u = .
gen d = .
gen lag = .

* Loop over the specified lags
forval h = 0/10 {

    * Generate the lagged variable for each country by year
    bysort country_id (year): gen rtfpna_lag`h' = L`h'.rtfpna

    * Run the fractional probit regression with lagged variable
    fracreg probit CS rtfpna_lag`h' i.year i.country_id, vce(robust)

    * Store the b, u, and d values for the current lag (store only for rtfpna_lag)
    replace b = _b[rtfpna_lag`h'] if _n == `h' + 1
    replace u = _b[rtfpna_lag`h'] + 2.576 * _se[rtfpna_lag`h'] if _n == `h' + 1
    replace d = _b[rtfpna_lag`h'] - 2.576 * _se[rtfpna_lag`h'] if _n == `h' + 1
    replace lag = `h' if _n == `h' + 1
}

* Remove missing values to keep only the filled rows
drop if missing(b)

* Plot the coefficients with confidence intervals
twoway (rcap u d lag, sort) ///
       (scatter b lag, sort), ///
       ytitle("Coefficient (b)") xtitle("Lag") ///
       title("Lagged Coefficients with 99% CI")

keep d u b lag

save al_TFP_lagged.dta,replace

export excel using "al_TFP_lagged.xlsx", sheet("Sheet1") firstrow(variables) replace
