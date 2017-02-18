/*
Version 1.0.1 (Feb 2, 2017)
Authors:  Jens Hainmueller, Jonathan Mummolo, Yiqing Xu (Maintainer)
*/

program interflex, eclass
version 12

   syntax varlist(numeric)  [if] [in] [aw fw pw] [, ///
     type(string) ///
     vce(string) ///
     CLuster(varname) ///
     REPS(integer 50) ///
     FE(varlist) ///
     Nbins(integer 3) ///
     CUToffs(numlist ascending) ///
     BW(real 0) ///
     seed(integer 12345678) ///
     grid(integer 20) ///
     neval(integer 50) ///
     TItle(string) ///
     XRange(numlist ascending) ///
     YRange(numlist ascending) ///
     XDistr(string) ///
     YLABel(string) ///
     DLABel(string) ///
     XLABel(string) ///
     SAVing(string) ///
     ]

   qui: marksample touse 
   gettoken Y varlist: varlist
   gettoken D varlist: varlist
   gettoken X varlist: varlist
   local y `Y'
   local x `X'
   local d `D'
   local z `varlist'
   local var `y' `d' `z'

   /* check input and set default */
   /* type */
   local type = cond("`type'"=="","binning","`type'")
   if ("`type'"!="binning" & "`type'"!="kernel" & "`type'"!="linear") {
      dis as err "type() invalid; choose either binning or kernel"
       exit
   }
   /* vce */
   local vce = cond("`vce'"=="","homoscedastic","`vce'")
   local vce = cond("`vce'"=="cl","cluster","`vce'")
   local vce = cond("`vce'"=="r","robust","`vce'")
   local vce = cond("`vce'"=="boot","bootstrap","`vce'")
   if ("`vce'"!="bootstrap" & "`vce'"!="off" & "`cluster'"!="") {
       local vce = "cluster"
   }
   if ("`vce'"!="homoscedastic" & "`vce'"!="robust" & ///
     "`vce'"!="cluster" & "`vce'"!="bootstrap" & ///
     "`vce'"!="off") {
       dis as err "vce() invalid"
       exit
   }
   if ("`FE'" == "") {
       local vce = cond("`vce'"=="homoscedastic","ols","`vce'")
   }
   else {
       local vce = cond("`vce'"=="homoscedastic","unadjusted","`vce'")
   }
   /* nbins */
   if (`nbins'<1) {
       dis as err "nbins() must be an integer greater than 1"
       exit
   }
   /* cluster */
   if ("`fe'"!=""&"`cluster'"=="") {
       disp as res "Fixed effects included; clustered standard errors highly recommended"
   }
   if ("`vce'"=="cluser"& "`cluster'"=="") {
       disp as error "Clustering variable not specified"
       exit
   }
   if ("`vce'"=="cluster" & "`cluster'"!="") {
       local vce = "cl " + "`cluster'"
   }
   if ("`vce'"=="off") {
      if ("`FE'" == "") {
         local vce = "ols"
         }
      else {
         local vce = "unadjusted"
         }
      local CIoff = 1
      }
   else {
      local CIoff = 0
      }
   
   /* nbins */
   if (`bw'<0) {
      dis as err "bw() must be a positive number"
       exit
   }
   /* grid */
   if (`grid'<=0) {
      dis as err "grid() must be a positive integer"
       exit
   }
   /* neval */
   if (`neval'<=0) {
       dis as err "eval() must be a positive integer"
       exit
   }
   /* cutoffs */
   local cutoffs = cond("`cutoffs'"=="","","`cutoffs'")
   /* xdistr */
   local xdistr = cond("`xdistr'"=="","histogram","`xdistr'")
   local xdistr = cond("`xdistr'"=="hist","histogram","`xdistr'")
   if ("`xdistr'"!="histogram" & "`xdistr'"!="density" ) {
       dis as err "xdistr() invalid; choose either histogram or density"
       exit
   }
   /* labels */
   local xlabel = cond("`xlabel'"=="","X","`xlabel'")
   local ylabel = cond("`ylabel'"=="","Y","`ylabel'")
   local dlabel = cond("`dlabel'"=="","D","`dlabel'")
   /* saving */
   if (strpos("`saving'",".") == 0) {
       local saving = "`saving'.pdf"
   }   

   /* weight */
   local wgt
   if ("`weight'" != "") {
       local wgt [`weight'`exp']
   }
   /*check bootstrap and weighted*/
   if ("`wgt'" != "" & "`vce'"=="bootstrap"){
       dis as err "Weights not supported by bootstrap"
       exit 
   }
   if ("`wgt'"!=""&"`type'"=="kernel") {
      disp as err "Weights not supported with the kernel estimation"
      exit
      }
  /*check X more than 5 values if type = kernel*/
  qui tab `x'
  local xcount = r(r)
  if("`type'"=="kernel"&`xcount'<5){
     disp as res "Moderator has less than 5 values; consider a fully saturated model."
  }

   preserve
   
   if("`in'"!="") {
       tempvar ran
       qui: gen `ran' = _n `in'
       qui drop if `ran' ==.
       qui drop `ran'
   }

   qui drop if `touse'==0

  if("`type'"!="kernel"){ // binning or linear
       tempvar dx
       gen `dx' = `x'*`d'
       local var `y' `d' `x' `dx' `z'
       if("`vce'"!="bootstrap"){
           if("`fe'"==""){
               qui: reg `var' `wgt', vce(`vce')
           }
           else {
               qui: reghdfe `var' `wgt', absorb(`fe') vce(`vce')    
           }
       }
       else {
           disp as txt "Bootstrapping..."
           set seed `seed'
           if("`fe'"==""){
               qui: bootstrap, reps(`reps') cl(`cluster'): reg `var'
           }
           else{
               qui: bootstrap, reps(`reps') cl(`cluster'): reghdfe `var', absorb(`fe') 
           }
       }
       mat coefs = e(b)
       local coefD = coefs[1,1]
       local coefX = coefs[1,2]
       local coefDX = coefs[1,3]

       mat vcov = e(V)
       local varD = vcov[1,1]
       local varDX = vcov[3,3]
       local covDX = vcov[3,1]
       cap mat drop coefs vcov
       cap drop `dx'
   }

   qui: summarize `x'
   local xmin = `r(min)'
   local xmax = `r(max)'
   local step = (`xmax' - `xmin')/(`neval'-1)
   mat eval = J(1,`neval',.)
   
   forvalues v = 1(1)`neval'{
       mat eval[1,`v']= `xmin' + (`v'-1)*`step'
   }
   
   qui tab `d'
   
   if(`r(r)'==2) {
       mat margeff = J(`neval',7,.)
       mat colnames margeff = xlevel marg se CI_l CI_u N Ntr
       local binary = 1
   }
   else {
       mat margeff = J(`neval',6,.)
       mat colnames margeff = xlevel marg se CI_l CI_u N
       local binary = 0
   }
   /*local nc = colsof(margeff)*/
   local quant = invnormal(0.975)

   /****** cross-validation *****/
   if("`type'"=="kernel" & "`bw'"=="0") {
       local CV = 1
       mat kernel = J(`grid',1,0)
       local first = log((`xmax' - `xmin')/10)
       local last = log(`xmax' - `xmin')
       local step = (`last'-`first')/`grid'
      
       mat gs = J(`grid',1,.)
       forvalues k = 1(1)`grid'{
           mat gs[`k',1] = exp(`first' + (`k'-1)*`step')
       }

       display as txt "Cross-validating bandwidth..."
      ******************

      set seed `seed'
      if ("`fe'"==""){
         tempvar f ff fold 
         gen `f' = uniform()
         sort `f'
         gen `ff' = _n
         gen `fold' = mod(`ff',5)+1
         drop `f' `ff'
         }
      else { // with fixed effects
         tempvar f ff fold 
         gen `f' = uniform()
         qui bysort `cluster' (`f'): replace `f' = `f'[1]
         sort `f'
         egen `ff' = group(`f')
         gen `fold' = mod(`ff',5)+1
         drop `f' `ff'
         foreach v of varlist `var' {
            g _bak_`v' = `v'
            qui reghdfe `v', absorb(`fe') res(_tmp_dm_`v')
            qui replace `v' = _tmp_dm_`v'
            }
         }

       /* number of covariates (Z) */
      local con: word count `z'
      local con1 = `con'+1

      forvalues k = 1(1)`grid'{
         forvalues i = 1(1)5{
            scalar fol = `i'
            mat coef1 = J(`neval',2+`con',.)
            
            forvalues j = 1(1)`neval'{
               tempvar xx wei dx
               qui: gen `xx' = `x'- el(eval,1,`j')
               qui: gen `wei' = sqrt(normalden(`xx'/el(gs,`k',1)))
               qui: gen `dx' = `xx'*`d'
               if ("`fe'"==""){
                  cap qui reg `var' `xx' `dx' [weight= `wei'] if `fold' != `i'
                  }
               else {
                  cap qui reghdfe `var' `xx' `dx' [weight= `wei'] if `fold' != `i', absorb(`fe') 
                  }
               if(_rc==0){
                  mat coe = e(b)
                  forvalues kk = 1(1)`con1'{
                     mat coef1[`j',`kk'] = el(coe,1,`kk')    
                     }
                  if("`fe'"==""){
                     mat coef1[`j',2+`con'] = el(coe,1,colsof(coe))
                     }
                  else{
                     mat coef1[`j',2+`con'] = 0
                     } 
                  cap drop `xx' `wei' `dx'
                  cap mat drop coe
                  }
               else {
                  mat kernel[`k',1] = .  
                  }
               }
               if(el(kernel,`k',1)!=.){
                  mata:sumofse("`y' `x' `d' `z' `fold'","eval","coef1","fol")
                  mat kernel[`k',1] = el(kernel,`k',1) + r(SE)
                  cap mat drop coef1
                  }
               else {
                  mat kernel[`k',1] = .
                  }
           }
         }
      local j = 1
      local coordinate = 1
      local min = kernel[1,1]
      while `j'<`grid'{
           local ++j
           if(el(kernel,`j',1)<`min'){
               local min = kernel[`j',1]
               local coordinate = `j'
           }
         }

      local bw = exp(`first'+`step'*(`coordinate'-1))
      mat cv = gs, kernel
      mat colnames cv = bw MSPE

      /* mpse contains missing value*/
      forvalues i =1(1)`grid'{
         if(el(kernel,`i',1)==.){
            display as txt "MPSE contains missing value"
            continue, break
            }  
         }

      **mat list cv
      display as txt "The optimal bandwidth is" as res %9.4f `bw'
      cap mat drop kernel gs

      /* retreive backup data */
      if ("`fe'"!=""){
         foreach v of varlist `var' {
            qui replace `v' = _bak_`v'
            }
         }
      cap drop _tmp_dm* _bak_*
      
   }
   else {
       local CV = 0
   }

   /*******calculation *********/   
   if("`type'"=="kernel"){
   
       if ("`vce'"!="bootstrap") {
           forvalues v = 1(1)`neval'{
               tempvar xx wei dx
               qui: gen `xx' = `x'- el(eval,1,`v')
               qui: gen `dx' = `xx'*`d'
               qui: gen `wei' = normalden(`xx'/`bw')
               if("`fe'"==""){
                   cap qui: reg `var' `xx' `dx' [weight = `wei'], vce(`vce')
               }
               else {
                   cap qui: reghdfe `var' `xx' `dx' [weight = `wei'], ///
                   absorb(`fe') vce(`vce')
               }
               mat margeff[`v',1] = el(eval,1,`v')
               if(_rc==0){
                  mat coe = e(b)
                  mat vcov = e(V)
                  mat margeff[`v',2] = el(coe,1,1)
                  mat margeff[`v',3] = sqrt(el(vcov,1,1))
                  mat margeff[`v',4] = el(coe,1,1)-`quant'*sqrt(el(vcov,1,1))
                  mat margeff[`v',5] = el(coe,1,1)+`quant'*sqrt(el(vcov,1,1))
                  }
               else {
                  forvalues i =2(1)5{
                     mat margeff[`v',`i'] =.
                     }
                  }
               qui count if `x' >= margeff[`v', 1] & `x' < (margeff[`v', 1] + `step')
               mat margeff[`v', 6] = `r(N)'
               if(`binary'==1){
                   qui count if `d' == 1 & `x' >= margeff[`v', 1] & `x' < (margeff[`v', 1] + `step')
                   mat margeff[`v', 7] = `r(N)'
               }
               drop `xx' `wei' `dx'
               cap mat drop coe vcov
           }
       } // end non-bootstrap
      
       else { // bootstrap
           disp as txt "Bootstrapping..."
           set seed `seed'

           forvalues v = 1(1)`neval'{
               tempvar xx wei dx
               qui: gen `xx' = `x'- el(eval,1,`v')
               qui: gen `dx' = `xx'*`d'
               qui: gen `wei' = normalden(`xx'/`bw')

               if("`fe'"==""){
                   qui: reg `var' `xx' `dx' [weight = `wei']
               }
               else {
                   cap qui: reghdfe `var' `xx' `dx'[weight = `wei'], absorb(`fe') 
               }
               drop `xx' `wei' `dx'
               mat margeff[`v',1] = el(eval,1,`v')
               if("`fe'"==""|_rc==0){
                  mat coe = e(b)
                  mat margeff[`v',2] = el(coe,1,1)
                  mat drop coe
                  }
               else{
                  mat margeff[`v',2] = .
               }

               mat bootce = J(`reps',1,.)
               // preserve and restore
               forvalues j = 1(1)`reps' {
                   tempvar xx wei dx
                   qui: gen `xx' = `x'- el(eval,1,`v')
                   qui: gen `dx' = `xx'*`d'
                   qui: gen `wei' = normalden(`xx'/`bw')
                   cap qui myboot `var' `xx' `dx'[weight = `wei'], fe(`fe') cluster(`cluster')
                   if(_rc==0){
                      mat coe = r(coe)
                      mat bootce[`j',1] = el(coe,1,1)
                      }
                   else{
                      dis as err "bootstrap results don't converge"
                      exit 430
                      }
                   cap drop `xx' `wei' `dx'
                   cap mat drop coe

               }
               mata:bce=st_matrix("bootce")
               mata:se=sqrt(variance(bce))
               mata:st_numscalar("r(se)",se)
               local bse = r(se)
               mat margeff[`v',3] = `bse'
               mat margeff[`v',4] = el(margeff,`v',2)-`quant'*`bse'
               mat margeff[`v',5] = el(margeff,`v',2)+`quant'*`bse'
               mat drop bootce

               qui count if `x' >= margeff[`v', 1] & `x' < (margeff[`v', 1] + `step')
               mat margeff[`v', 6] = `r(N)'
               if(`binary'==1){
                   qui count if `d' == 1 & `x' >= margeff[`v', 1] & `x' < (margeff[`v', 1] + `step')
                   mat margeff[`v', 7] = `r(N)'
               }     
           } 

       } // end bootstrap 

       qui count if `x' == `xmax'
       mat margeff[`neval', 6] = margeff[`neval', 6] + `r(N)'
       if(`binary'==1){
           qui count if `d' == 1 & `x' == `xmax'
           mat margeff[`neval', 7] = margeff[`neval', 7] + `r(N)'
       } 
   } // end type == kernel

   else { //binning
       qui forvalues i = 1(1)`neval' {
           mat margeff[`i', 1] = `xmin' + (`i'-1)* `step'
           mat margeff[`i', 2] = `coefD' + `coefDX' * margeff[`i',1]
           mat margeff[`i', 3] = sqrt(`varD' + margeff[`i',1]^2 * `varDX' + `covDX' * margeff[`i',1] * 2)
           mat margeff[`i', 4] = margeff[`i', 2] - `quant' * margeff[`i', 3]
           mat margeff[`i', 5] = margeff[`i', 2] + `quant' * margeff[`i', 3]
           /* est +/- 1.96 s.e. only work for relatively large samples */
           count if `x' >= margeff[`i', 1] & `x' < (margeff[`i', 1] + `step')
           mat margeff[`i', 6] = `r(N)'
           if(`binary'==1){
               count if `d' == 1 & `x' >= margeff[`i', 1] & `x' < (margeff[`i', 1] + `step')
               mat margeff[`i', 7] = `r(N)'
           }
       }
       qui count if `x' == `xmax'
       mat margeff[`neval', 6] = margeff[`neval', 6] + `r(N)'
       if(`binary'==1){
           qui count if `d' == 1 & `x' == `xmax'
           mat margeff[`neval', 7] = margeff[`neval', 7] + `r(N)'
       } 
       if("`type'"=="binning"){
           if("`cutoffs'"==""){
               local ncuts = `nbins'-1
               mat cuts = J(`nbins'+1, 1, .)
               mat cuts[1,1] = `xmin'
               mat cuts[`nbins'+1,1] = `xmax'+ 0.01
               qui forvalues i=1/`ncuts' {
                   local pct = 100/`nbins'*`i'
                   _pctile `x', p(`pct')
                   mat cuts[(`i'+1),1] = `r(r1)'
               }
           }
           else {
               local binsize = wordcount("`cutoffs'")
               local nbins = `binsize'+1
               mat cuts = J(`nbins'+1, 1, .)
               mat cuts[1,1] = `xmin'
               local k = 1
               qui forvalues i=1/`binsize' {
                   gettoken xx cutoffs: cutoffs
                   if (`xx'>`xmin' & `xx'<`xmax') {
                       mat cuts[(`k'+1),1] = `xx'
                       local k = `k' + 1
                   }
               }
               mat cuts[(`k'+1),1] = `xmax'+ 0.01
               local j = `k'+1
               mat cuts = cuts[1..(`k'+1),1]
               local nbins = `k'
           }

           /* x level, marginal effect, SE, lower CI, upper CI, N, Ntr */
           mat MEbin = J(`nbins', 5, .) 
           mat colnames MEbin = x0 bin_marg bin_se bin_CI_l bin_CI_u
           local binvar
           qui forvalues i = 1/`nbins' {
               g _g`i' = (`x' >= cuts[`i',1] & `x' < cuts[(`i'+1),1])
               count if _g`i'==1
               if (`r(N)'>0) {
                   _pctile `x' if _g`i' ==1, p(50)
                   mat MEbin[`i', 1] = `r(r1)'
               }
               else {
                   mat MEbin[`i', 1] = (cuts[`i',1] + cuts[(`i'+1),1])/2
               }
               g _g`i'_d = _g`i' * `d'
               g _g`i'_x = _g`i' * (`x' - MEbin[`i',1])
               g _g`i'_dx = _g`i' * `d' * (`x' - MEbin[`i',1])
               local binvar `binvar' _g`i'_d _g`i' _g`i'_x _g`i'_dx
           }
         
           if("`vce'"!="bootstrap"){
               if("`fe'"==""){
                   qui: reg `y' `binvar' `z' `wgt', ///
                   noconstant vce(`vce')
               }
               else{
                   qui: reghdfe `y' `binvar' `z' `wgt', ///
                   absorb(`fe') vce(`vce')    
               }
           }    
           else {  // bootstrap
               if("`fe'"==""){
                   qui: bootstrap, reps(`reps') cl(`cluster'): reg `y' `binvar' ///
                   `z', noconstant 
               }
               else {
                   qui: bootstrap, reps(`reps') cl(`cluster'): reghdfe `y' `binvar' ///
                   `z', absorb(`fe')  
               }
           }
          /* do not forget to put in the same set of control variables */

         /* storage */
         mat coefs = e(b) 
         mat vcov = e(V)
         qui forvalues i = 1/`nbins' {
            mat MEbin[`i',2] = coefs[1,1+4*(`i'-1)]
            mat MEbin[`i',3] = sqrt(vcov[1+4*(`i'-1),1+4*(`i'-1)])
            mat MEbin[`i', 4] = MEbin[`i', 2] - `quant' * MEbin[`i', 3]
            mat MEbin[`i', 5] = MEbin[`i', 2] + `quant' * MEbin[`i', 3]
            }
         
         /* wald test */
         if (`nbins'>1) { 
            local waldvar
            qui forvalues i = 2/`nbins' {
               local waldvar `waldvar' _g`i'_d _g`i' _g`i'_x _g`i'_dx
               }
            tempvar dx
            g `dx' = `d'*`x'
            if ("`vce'"!="bootstrap"){
               if("`fe'"==""){
                  qui: reg `y' `d' `x' `dx' `waldvar' `z' `wgt', ///
                    vce(`vce')
                  }
               else{
                  qui: reghdfe `y' `d' `x' `dx' `waldvar' `z' `wgt', ///
                    absorb(`fe') vce(`vce')    
                  }
               }    
            else {  // bootstrap
               if("`fe'"==""){
                  qui: bootstrap, reps(`reps') cl(`cluster'): reg `y' `d' `x' `dx' `waldvar' ///
                    `z'
                  }
               else {
                  qui: bootstrap, reps(`reps') cl(`cluster'): reghdfe `y' `d' `x' `dx' `waldvar' ///
                    `z', absorb(`fe')  
                  }
               }
            local wcount = wordcount("`waldvar'")
            local waldtest
            qui forvalues i = 1/`wcount'{
               gettoken wv waldvar:waldvar
               local waldtest `waldtest' (`wv'=0)
               }
            qui: test `waldtest'
            if("`vce'"!="bootstrap"){
               local pwald = r(p)
               local F = r(F)
               display as text "p value of Wald test: "   as res %5.4f `pwald'
               }
            else {
               display _newline
               local pwald = r(p)
               local C = r(chi2)
               display as text "p value of Wald test: " as res %5.4f `pwald'
               }
            }
         }
      }
   


   ***************************************************************
   /* plotting */
   svmat margeff, names(col)
   if(`binary'==1){
       qui g Nco = N - Ntr
   }
   if("`type'"=="binning"){
       svmat MEbin, names(col)
   }

   /* adjustment for histogram and density plot*/
   qui sum CI_u if xlevel!=.
   local ymax = `r(max)'
   qui sum CI_l if xlevel!=.
   local ymin = `r(min)'
   if ("`type'"=="binning") {
       qui sum bin_CI_u
       local ymax = max(`ymax',`r(max)')
       qui sum bin_CI_l
       local ymin = min(`ymin',`r(min)')
   }
   if("`yrange'"!="") {
      gettoken y1 y2:yrange
      local ymax = max(`y2',`ymax')
      local ymin = min(`y1',`ymin')
      local ydiff = `ymax' - `ymin'
      local base = `ymin'
      }
   else {
      local ydiff = `ymax' - `ymin'
      local base = `ymin' - `ydiff'/5
      }
   

   if("`xdistr'"=="histogram"){
      tempvar N_adj Ntr_adj
      qui sum N if xlevel!=.
      local Nmax = `r(max)'
      qui g `N_adj' = (N/`Nmax')*`ydiff'/5 + `base'
      if(`binary'==1){
         qui g `Ntr_adj' = (Ntr/`Nmax')*`ydiff'/5 + `base'
         }
      local Xdis1 bar `N_adj' xlevel, barwidth(0.08) ///
        bfcolor(gs10) blcolor(gs8) base(`base')
      local Xdis2 bar `Ntr_adj'  xlevel, barwidth(0.08) ///
        bfcolor(red) blcolor(gs8) base(`base') 
      }
   else {  // density 
      if(`binary'==1){
         tempvar dbase Nco1 Ntr1 Nco2 Ntr2 Nco_adj Ntr_adj
         qui gen `dbase'=`base'
         qui kdensity `x' if `d'==0 & xlevel ==., kernel(gaussian) ///
           gen(`Nco1' `Nco2') nograph
         qui kdensity `x' if `d'==1 & xlevel ==., kernel(gaussian) ///
           gen(`Ntr1' `Ntr2') nograph
         qui sum `Nco2'
         local Nmax = `r(max)'
         qui g `Nco_adj' = (`Nco2'/`Nmax')*`ydiff'/5 + `base'
         qui sum `Ntr2'
         local Nmax = `r(max)'
         qui g `Ntr_adj' = (`Ntr2'/`Nmax')*`ydiff'/5 + `base'
         local Xdis1 line `Nco_adj' `dbase' xlevel, ///
           lcolor(gs8)
         local Xdis2 line `Ntr_adj' `dbase' xlevel, ///
           lcolor(red)
         }
      else {
         tempvar dbase N1 N2 N_adj
         qui: gen `dbase'=`base'
         qui kdensity `x', kernel(gaussian) gen(`N1' `N2') nograph
         qui sum `N2'
         local Nmax = `r(max)'
         qui g `N_adj' = (`N2'/`Nmax')*`ydiff'/5 + `base'
         local Xdis1 rarea `N_adj' `dbase' xlevel, ///
           fcolor(none) lcolor(gs8)
         }
      }
   
   if("`type'"=="binning"){
      if (`nbins'==3) {
         local xlow = MEbin[1,1]
         local xmid = MEbin[2,1]
         local xhigh = MEbin[3,1]
         local lpos = `ymax' - `ydiff'/10
         if(`binary'==1){
            twoway (rarea CI_u CI_l xlevel, astyle(ci)) ///
              (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              (`Xdis2') ///
              (rcap bin_CI_u bin_CI_l x0 if bin_se!=0, lcolor(red)) ///
              (scatter bin_marg x0 if bin_se!=0, mcolor(red) ms(O)) ///
              (scatter bin_marg x0 if bin_se==0, mcolor(blue) ms(Sh)) ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              text(`lpos' `xlow' "L", color(gs10) size(*1.5)) ///
              text(`lpos' `xmid' "M", color(gs10) size(*1.5)) ///
              text(`lpos' `xhigh' "H", color(gs10) size(*1.5)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid) xlabel(#5)  legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         else {
            twoway (rarea CI_u CI_l xlevel, astyle(ci)) ///
              (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              (rcap bin_CI_u bin_CI_l x0 if bin_se!=0, lcolor(red)) ///
              (scatter bin_marg x0 if bin_se!=0, mcolor(red) ms(O)) ///
              (scatter bin_marg x0 if bin_se==0, mcolor(blue) ms(Sh)) ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              text(`lpos' `xlow' "L", color(gs10) size(*1.5)) ///
              text(`lpos' `xmid' "M", color(gs10) size(*1.5)) ///
              text(`lpos' `xhigh' "H", color(gs10) size(*1.5)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid) xlabel(#5) legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         }
      else { // 2 or 4 or more bins
         if(`binary'==1){
            twoway (rarea CI_u CI_l xlevel, astyle(ci)) ///
              (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              (`Xdis2') ///
              (rcap bin_CI_u bin_CI_l x0 if bin_se!=0, lcolor(red)) ///
              (scatter bin_marg x0 if bin_se!=0, mcolor(red) ms(O)) /// 
              (scatter bin_marg x0 if bin_se==0, mcolor(blue) ms(Sh)) ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid)  xlabel(#5)  legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         else {
            twoway (rarea CI_u CI_l xlevel, astyle(ci)) ///
              (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              (rcap bin_CI_u bin_CI_l x0 if bin_se!=0, lcolor(red)) ///
              (scatter bin_marg x0 if bin_se!=0, mcolor(red) ms(O)) ///
              (scatter bin_marg x0 if bin_se==0, mcolor(blue) ms(Sh)) ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid)  xlabel(#5)  legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         }
      }
   else { // kernel
      if(`CIoff'==0){
         if(`binary'==1){
            twoway (rarea CI_u CI_l xlevel, astyle(ci)) ///
              (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              (`Xdis2') ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid)  xlabel(#5)  legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         else {
            twoway (rarea CI_u CI_l xlevel, astyle(ci)) ///
              (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid)  xlabel(#5)  legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         }
      else { // without CI
         if(`binary'==1){
            twoway (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              (`Xdis2') ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid)  xlabel(#5)  legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         else {
            twoway (line marg xlevel, lstyle(p1)) ///
              (`Xdis1') ///
              , yline(0, lcolor(gs14) lwidth(1.2)) ///
              xtitle("Moderator: `xlabel'") ///
              ytitle("Marginal Effect of `dlabel' on `ylabel'") ///
              ylabel(#6,ang(h) nogrid)  xlabel(#5) legend(off) ///
              xsc(r(`xrange')) ysc(r(`yrange')) title("`title'")
            }
         }
      }

   if("`saving'"!=""){
       cap graph export `saving',replace
   }

   if("`type'"=="kernel" & `CIoff'==1){
       mata:marge=st_matrix("margeff")
       mata:marge=marge[.,(1::3)]
       mata:st_matrix("r(m)",marge)
       mat margeff = r(m)
       mat colnames margeff =xlevel marg se
   }
   else {
       mata:marge=st_matrix("margeff")
       mata:marge=marge[.,(1::5)]
       mata:st_matrix("r(m)",marge)
       mat margeff = r(m)
       mat colnames margeff =xlevel marg se CI_l CI_u
   }

   restore
   ereturn clear
   ereturn matrix margeff = margeff
   if ("`type'"=="binning") {
      ereturn matrix estBin = MEbin
      }
   if ("`type'"=="binning" & `nbins'>1){
       ereturn scalar pwald = `pwald'   
   }
   if ("`type'"=="kernel"){
       ereturn scalar bandiwth = `bw'
       if(`CV'==1){
           ereturn matrix CVout = cv   
       }
   }
      
end

*********************************************************bootstrap
program myboot, rclass
version 12
  
  syntax varlist(numeric) [aw fw pw] [, fe(varlist) cluster(varname) constant(string)]
  local wgt
  if ("`weight'" != "") {
       local wgt [`weight'`exp']
  }
  preserve
  local constant = cond("`constant'"=="","yes","no")
  bsample ,cluster(`cluster')
  if("`fe'"==""){
      if("`constant'"=="yes"){
          qui regress `varlist' `wgt'   
      }
      else{
          qui regress `varlist' `wgt', noconstant 
      }
  }
  else{
     qui reghdfe `varlist' `wgt',absorb(`fe') 
  }
  mat coe = e(b)
  restore
  return mat coe=coe

end
*********************************************************mata
version 12
mata:
   mata set matastrict on 
   void sumofse (string scalar varname, string scalar name1, string scalar name2, string scalar fold)
   {
      f = st_numscalar(fold)
      test = st_data(.,varname)
      c = rows(test)
      p = cols(test)
      ii = 0
      for(j=1;j<=c;j++){ 
          if(test[j,p]==f){
              ii = ii + 1
          }
      }
      r = J(1,ii,.)
      jj = 0
      for(j=1;j<=c;j++){ 
          if(test[j,p]==f){
              jj = jj + 1
              r[1,jj] = j
          }
      }

      y = test[r,1]
      x = test[r,2]
      /*independent variables in test set*/ 
      dep = test[r,3::(p-1)],J(ii,1,1)
      coe = J(rows(dep),cols(dep),.)
      
      /*eveal*/
      xx = st_matrix(name1)
      /*estimated coefficents*/
      coef = st_matrix(name2)

      /*2nn*/
      mar = J(ii,2,.)
      for(i=1;i<=ii;i++){
          dif = J(cols(xx),1,.)
          for(j=1;j<=cols(xx);j++){
              dif[j,1]=abs(x[i,1]-xx[1,j])
          }
          mini = min(dif)
          for(k=1;k<=cols(xx);k++){
              if(dif[k,1]==mini){
                  mar[i,1]=k
                  dif[k,1]=.
                  break;
              }
          }
          mini2 = min(dif)
          for(k=1;k<=cols(xx);k++){
              if(dif[k,1]==mini2){
                  mar[i,2]=k
                  break;
              }
          }
      }
      
      for(i=1;i<=ii;i++){
          for(j=1;j<=cols(dep);j++){
              coe[i,j] = 0.5*(coef[mar[i,1],j]+coef[mar[i,2],j])
          }
      }

      yy = rowsum(dep :* coe)
      SE = sum((y-yy):*(y-yy))/ii
      st_numscalar("r(SE)", SE)
   }
end
