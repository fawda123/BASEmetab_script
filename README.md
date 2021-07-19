# README

Modification of materials for running BASEmetab.  All content derived from the original repository [here](https://github.com/dgiling/BASEmetab).

Comparison of metabolism results from [Fwoxy](https://github.com/jmarriola/fwoxy), [Odum](https://github.com/fawda123/WtRegDO), and BASEmetab using Fwoxy estimated DO [here](https://fawda123.github.io/BASEmetab_script/fwoxy_comp).

Comparison of BASEmetab with Odum metabolism results with 2012 Apalachicola data for different iterations of BASEmetab generated using the scripts 1-5 [here](https://fawda123.github.io/BASEmetab_script/comp_plots).

## R scripts 

* `script1.R` Run stripped down BASEmetab code

* `script2.R` Run BASEmetab using the original package

* `script3.R` Run stripped down BASEmetab code with instantaneous K

* `script4.R` Run stripped down BASEmetab code with instantaneous K and metabolic day

* `script5.R` Run stripped down BASEmetab code with instantaneous K, metabolic day, p and theta estimated

* `fwoxy_basemetab.R` Use [Fwoxy](https://github.com/jmarriola/fwoxy) estimated DO to estimate metabolism flux rates with BASEmetab, as in `script5.R`

* `method_comp.R` Comparisons of Odum and BASEmetab methods for the different attempts

* `dat_proc.R` Data prep for interoperability between packages

* `funcs.R` Misc. helper functions