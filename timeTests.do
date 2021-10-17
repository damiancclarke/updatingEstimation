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
global NTESTS = 5

dis `maxmem'
local vars = 5
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

  foreach num of numlist 1(1)$NTESTS {
      timer clear `num'
      timer on `num'
      cumulativels y `xvars', filename("rawdata2.raw") blocksize(`blocks')
      timer off `num'
  }

  foreach num of numlist 1(1)$NTESTS {
      qui timer list
      local times = r(t`num')
      return local testtime`num' = `times'
  }	    
  restore
  dis "good bye"
end

set obs 100
gen Nobs = .
foreach num of numlist 1(1)$NTESTS {
    gen time`num'=.
}
gen regs = .

local i = 1
foreach fsize of numlist 0.025(0.025)0.9 {
    dis "fsize is `fsize'"
    !date 
    cap rm "rawdata2.raw"

    !free
    local fsize = round((`fsize'*`maxobs')/10)*10
    local a = `fsize'/10
    dis `a'
    dis "here"
    gettime `fsize' `nvars' `a'
    replace Nobs = `fsize' in `i'
    foreach num of numlist 1(1)$NTESTS {
    	replace time`num'= `r(testtime`num')' in `i'
    }
    replace regs = 0 in `i'
    local ++i
    
    replace Nobs = `fsize' in `i'
    replace regs = 1 in `i'
    foreach num of numlist 1(1)$NTESTS { 
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
    egen meanTime=rowmean(time1-time$NTESTS)
    local end SE
    if c(MP)==1 local end MP
    #delimit ;
    replace Nobs = Nobs/1000;
    format %9.0gc Nobs;
    twoway connected time1 Nobs if regs==0, 
        || connected time1 Nobs if regs==1, ms(Sh) lpattern(dash)
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
    
    keep if time1!=.
    save "$OUT/times`end'_`nvars'.dta", replace
    list
    restore
}
cap rm "rawdata2.raw"
log close
