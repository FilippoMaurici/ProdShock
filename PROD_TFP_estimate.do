clear

cd C:\Users\filip\OneDrive\Desktop\Lavori\Miei\PROD_TFP

use "C:\Users\filip\OneDrive\Desktop\Lavori\Miei\PROD_TFP\pwt1001.dta" 

encode country, gen(country_id)
drop country

rename country_id country

sort country year

xtset country year

drop countrycode pl_n pl_m pl_x pl_g pl_i pl_c csh_r csh_m csh_x csh_g statcap csh_c csh_i cor_exp

