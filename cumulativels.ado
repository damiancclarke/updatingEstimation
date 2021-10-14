*! cumulativels: Cumulative Least Squares for large databases
*! version 0.0.0 August 1, 2019 @ 12:59:11


cap program drop cumulativels
program cumulativels, eclass
version 15.0

#delimit ;
syntax namelist(min=2)
       ,
       filename(string)
       BLOCKsize(integer)
       ;
#delimit cr
local K1 : word count `namelist'
local K = `K1'-1

mata: A=runcls("`filename'",`blocksize',`K')
mata: st_matrix("betas", A)

***USE POINTERS WHEN PASSING THINGS TO updateCLS

ereturn matrix CLSbeta = betas
end


version 15.0
mata:

struct resultsCLS {
    real scalar yy
    real colvector Xy
    real matrix XX
}

numeric matrix runcls(string scalar filename, real scalar blocksize, real scalar K)
{
    fh = fopen(filename,"r")
    names = fget(fh)
    X = J(blocksize,K,.)
    y = J(blocksize,1,.)
    i = 1
    block = 1
    struct resultsCLS scalar r1
    r1.yy    = 0
    r1.Xy    = J(K,1,0)
    r1.XX    = J(K,K,0)
    N        = 0
    
    while((line=fget(fh)) != J(0,0,"")) {
        vecline = strtoreal(tokens(subinstr(line,",", " ")))
        X[i,.] = vecline[2..K+1] 
        y[i,1] = vecline[1]
        if (i==blocksize) {
            
            Xp = &X
            yp = &y
            r1=updateCLS(Xp,yp,r1)
            i=0
            ++block
            X = J(blocksize,K,.)
            y = J(blocksize,1,.)
        }
        ++i
        ++N
    }
    if (y[1]!=.) { 
        y=y[1..i-1]
        X=X[1..i-1,1..K]
        r1=updateCLS(X,y,r1)
    };
    fclose(fh)
    beta = invsym(r1.XX)*r1.Xy
    uPu  = (r1.yy - 2*beta'*r1.Xy + beta'*r1.XX*beta)
    vcov = uPu/(N-K)*invsym(r1.XX)
    ses  = sqrt(diagonal(vcov))
    output = (beta, ses)
    output
    return(beta)
}

struct resultsCLS scalar updateCLS(pointer X, pointer y,
                             struct resultsCLS scalar r) {
    struct resultsCLS scalar a
    a.yy = r.yy+quadcross(*y,*y)
    a.Xy = r.Xy+quadcross(*X,*y)
    a.XX = r.XX+quadcross(*X,*X)
    return(a)
}
end


