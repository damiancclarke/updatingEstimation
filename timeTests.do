/* timerTests.do                                           yyyy-mm-dd:2019-07-05
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

This file runs some tests using various file sizes with the updatedLS algorithm.
It is currently calling a not necessarily optimal mata script.  I have tried to
be careful with I/O, but not checked all steps for optimality in mata...


https://www.stata.com/stata14/huge-datasets/
 obs = (memory used)/(width+24)*(1024^3)/(1000^3)

MAX obs = (memory)/(width)*(1024^3)/(1000^3)

width is memory required to store an observation.

https://www.stata.com/manuals13/u6.pdf
https://www.stata.com/manuals/ddatatypes.pdf
- byte uses 1 byte
- int uses 2 bytes
- long uses 4 bytes
- float uses 4 bytes
- double uses 8 bytes

So, say typical observation made up of 10 float variables, then width is 40

For a computer with 16GB of RAM this is:
 (16000000000/40)*(1024^3/1000^3)
  429,500,000

On unix see memory !free | grep Mem:
 ** Automate to check if sufficient memory is available for regress

To run with no hangup on server:
nohup /usr/local/stata/stata-mp -b do timeTests2.do &

*/

*-------------------------------------------------------------------------------
*--- (0) Setup and locals
*-------------------------------------------------------------------------------    
clear all
cap log close
do updatedls.ado
do cumulativels.ado

global OUT "/root/results"

log using timeTests.txt, text replace

local 1GB    = 1000000000
local maxmem = 1*`1GB'

dis `maxmem'
local vars = 10
*local vars = 5
local maxobs = (`maxmem'/(4*(`vars'+1)))*(1024^3/1000^3)
dis `maxobs'


*-------------------------------------------------------------------------------
*--- (1) Test 1: Vary N obs, hold N var constant
*-------------------------------------------------------------------------------    
local nvars = `vars'
local XVARS
foreach k of numlist 1(1)`nvars' {
    local XVARS `XVARS' x`k'
}	

cap program drop gettime
program gettime, rclass
version 15.0
  args n k blocks 

  dis "hello"
  preserve
  clear
  set obs `n'
  gen x1 = 1
  gen y  = 5*x1
  local xvars x1
  foreach num of numlist 2(1)`k' {
      gen x`num'= runiform()
      replace y = y+2*x`num'
      local xvars `xvars' x`num'
  }	
  replace y = y + rnormal()*3
  order y 
  !free
  outsheet using "rawdata2.raw", delim(",")

  foreach num of numlist 1(1)10 {
      timer clear `num'
      timer on `num'
      cumulativels y `xvars', filename("rawdata2.raw") blocksize(`blocks')
      timer off `num'
  }
  timer list
  local times = r(t1)
  
  return local testtime = `times'
  foreach num of numlist 2(1)10 {
      qui timer list
      local times = r(t`num')
      return local testtime`num' = `times'
  }	    
  restore
  dis "good bye"
end

set obs 100
gen Nobs = .
gen time = .
foreach num of numlist 2(1)10 {
    gen time`num'=.
}
gen regs = .

local i = 1
foreach fsize of numlist 0.025(0.025)0.8 {
    dis "fsize is `fsize'"
    !date 
    cap rm "rawdata2.raw"

    !free
    local fsize = round((`fsize'*`maxobs')/10)*10
    local a = `fsize'/10
    dis `a'
    dis "here"
    **gettime `fsize' 10 `a'
    gettime `fsize' `nvars' `a'
    replace Nobs = `fsize' in `i'
    replace time = `r(testtime)' in `i'
    foreach num of numlist 2(1)10 {
    	replace time`num'= `r(testtime`num')' in `i'
    }
    replace regs = 0 in `i'
    local ++i
    
    preserve
    clear
    timer clear 1
    timer on 1
    insheet using "rawdata2.raw", names
    reg y `XVARS', nocons
    timer off 1
    timer list
    restore
    replace Nobs = `fsize' in `i'
    replace time = r(t1) in `i'
    replace regs = 1 in `i'

    foreach num of numlist 2(1)10 { 
        preserve
        clear
        timer clear `num'
        timer on `num'
        insheet using "rawdata2.raw", names
        reg y `XVARS', nocons
        timer off `num'
        timer list
        local t = r(t`num')
	restore
        replace time`num' = r(t`num') in `i'
    }	
    local ++i    

    preserve
    egen meanTime=rowmean(time time2 time3 time4 time5 time6 time7 time8 time9 time10)
    local end SE
    if c(MP)==1 local end MP
    #delimit ;
    replace Nobs = Nobs/1000;
    format %9.0gc Nobs;
    twoway connected time Nobs if regs==0, 
        || connected time Nobs if regs==1, ms(Sh) lpattern(dash)
    legend(order(1 "Cumulative LS" 2 "Stata regress")
           position(6) rows(1)) ytitle("Execution Time (in seconds)") 
    xtitle("Observations (in thousands)") scheme(plotplain);
    graph export "$OUT/NobsTime_`end'_`nvars'.eps", replace;
    
    twoway connected meanTime Nobs if regs==0, 
        || connected meanTime Nobs if regs==1, ms(Sh) lpattern(dash)
    legend(order(1 "Cumulative LS" 2 "Stata regress")
           position(6) rows(1)) ytitle("Execution Time (in seconds)") 
    xtitle("Observations (in thousands)") scheme(plotplain);
    graph export "$OUT/NobsMeanTime_`end'_`nvars'.eps", replace;
    #delimit cr
    
    keep if time!=.
    save "$OUT/times`end'_`nvars'.dta", replace
    list
    restore
}
cap rm "rawdata2.raw"




exit




























xscale(log)



**Below imports free memory to Stata
    *preserve
    *!free > mem.txt
    *insheet using mem.txt, clear
    *split v1
    *destring v14, force replace
    *sum v14 in 2
    *local realmem = r(mean)*1000
    *local vars = 11
    *local newobs  = (`realmem'/(4*`vars'))*(1024^3/1000^3) 
    *restore
    *local test = `maxobs'*`fsize'
    *dis "Required obs is: `test'. Feasible obs is: `newobs'."
    *if `test'>`newobs' continue, break










exit
cap rm "rawdata2.raw"
local i = 1
foreach size of local filesizes {
    set obs `size'
    gen x1 = 1
    gen y  = 5*x1
    local xvars x1
    foreach num of numlist 2(1)`nvars' {
        gen x`num'= runiform()
        replace y = y+2*x`num'
        local xvars `xvars' x`num'
    }
    replace y = y + rnormal()*3
    order y 
    outsheet using "rawdata2.raw", delim(",")
    clear
    
    local bsize =  `size'/10
    timer on `i'
    cumulativels y `xvars', filename("rawdata2.raw") blocksize(`bsize')
    timer off `i'
    local ++i

    timer on `i'
    insheet using "rawdata2.raw", names
    reg y `xvars', nocons
    timer off `i'
    local ++i
    clear
    
    rm "rawdata2.raw"
}
timer list

set obs 100
gen nobs = 1000 in 1/2
replace nobs = 10000 in 3/4
replace nobs = 100000 in 5/6
replace nobs = 200000 in 7/8
replace nobs = 300000 in 9/10
replace nobs = 400000 in 11/12
gen time = .
foreach num of numlist 1(1)12 { 
    replace time = r(t`num') in `num'
}
gen cumls = 1
replace cumls = 0 in 2
replace cumls = 0 in 4
replace cumls = 0 in 6
replace cumls = 0 in 8
replace cumls = 0 in 10
replace cumls = 0 in 12

twoway connected time nobs if cumls == 1 || connected time nobs if cumls==0
graph export "$OUT/tester.eps", replace

/*
THINGS TO VARY
 - N vars
 - Memory
 - partition size


*/







local filesizes 1000 10000 100000 1000000 10000000 100000000 200000000 300000000

cap rm "rawdata2.raw"
local i = 1
foreach size of local filesizes {
    clear all
    set obs `size'
    gen x1 = 1
    foreach num of numlist 2(1)6 {
        gen x`num'= runiform()
    }
    gen y = 5*x1 + 4*x2 + 3*x3 + 2*x4 + 1*x5 + 0.5*x6 + rnormal()*3
    order y x1 x2 x3 x4 x5 x6
    outsheet using "rawdata2.raw", delim(",")
    clear
    local bsize = `size'/10
    timer on `i'
    updatedls y x1 x2 x3 x4 x5 x6, filename("rawdata2.raw") blocksize(`bsize')
    timer off `i'
    local ++i

    local bsize = `size'/20
    timer on `i'
    updatedls y x1 x2 x3 x4 x5 x6, filename("rawdata2.raw") blocksize(`bsize')
    timer off `i'
    local ++i

    local bsize = 1000
    timer on `i'
    updatedls y x1 x2 x3 x4 x5 x6, filename("rawdata2.raw") blocksize(`bsize')
    timer off `i'
    local ++i

    local bsize = `size'
    timer on `i'
    updatedls y x1 x2 x3 x4 x5 x6, filename("rawdata2.raw") blocksize(`bsize')
    timer off `i'
    local ++i

    timer on `i'
    insheet using "rawdata2.raw", names
    reg y x1 x2 x3 x4 x5 x6, nocons
    timer off `i'
    local ++i

    timer list
    rm "rawdata2.raw"
}

log close
