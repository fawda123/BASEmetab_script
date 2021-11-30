# README

Modification of materials for running BASEmetab.  All content derived from the original repository [here](https://github.com/dgiling/BASEmetab).

Comparison of metabolism results from [Fwoxy](https://github.com/jmarriola/fwoxy), [Odum](https://github.com/fawda123/WtRegDO), BASEmetab, and [EBASE](https://github.com/fawda123/EBASE) using Fwoxy estimated DO [here](https://fawda123.github.io/BASEmetab_script/fwoxy_comp).

Comparison of BASEmetab, EBASE, and Odum metabolism results with 2012 Apalachicola generated using the scripts 1-6 [here](https://fawda123.github.io/BASEmetab_script/comp_plots).

## R scripts 

* `script1.R` Run stripped down BASEmetab code with 2012 Appalachicola data

* `script2.R` Run BASEmetab using the original package with 2012 Appalachicola data

* `script3.R` Run stripped down BASEmetab code with instantaneous K with 2012 Appalachicola data

* `script4.R` Run stripped down BASEmetab code with instantaneous K and metabolic day with 2012 Appalachicola data

* `script5.R` Run stripped down BASEmetab code with instantaneous K, metabolic day, p and theta estimated with 2012 Appalachicola data

* `script6.R` Run EBASE with with 2012 Appalachicola data

* `fwoxy_basemetab.R` Use [Fwoxy](https://github.com/jmarriola/fwoxy) estimated DO to estimate metabolism flux rates with BASEmetab, as in `script5.R`

* `fwoxy_ebase.R` Use [Fwoxy](https://github.com/jmarriola/fwoxy) estimated DO to estimate metabolism flux rates with EBASE, as in `script6.R`

* `method_comp.R` Comparisons of Odum and BASEmetab methods for the different attempts

* `dat_proc.R` Data prep for interoperability between packages

* `funcs.R` Misc. helper functions