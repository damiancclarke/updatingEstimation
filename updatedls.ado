*! updatedls: Updated Least Squares for large databases
*! version 0.0.0 July 1, 2019 @ 15:12:20

/*
This was useful for testing interactively -- not now.
struct results scalar myfcn2() {
    struct results scalar a
    a.beta = J(6,1,0)
    a.Sigma= J(6,6,0)
    return(a)
}
r1=myfcn2()
Omega   = I(K)-invsym(I(K)+SM1*SP1)*SM1*SP1
*/

cap program drop updatedls
program updatedls, eclass
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

mata: runuls("`filename'",`blocksize',`K')

end

version 15.0
mata:

struct results {
    real colvector beta
    real matrix Sigm1
    real scalar yy
    real colvector Xy
    real matrix XX
}

void runuls(string scalar filename, real scalar blocksize, real scalar K)
{
    fh = fopen(filename,"r")
    names = fget(fh)
    X = J(blocksize,K,.)
    y = J(blocksize,1,.)
    i = 1
    block = 1
    struct results scalar r1
    r1.beta  = J(K,1,1)
    r1.Sigm1 = J(K,K,1)
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
            r1=update(Xp,yp,block,K,r1)
            i=0
            ++block
            X = J(blocksize,K,.)
            y = J(blocksize,1,.)
        }
        ++i
        ++N
    }
    fclose(fh)
    uPu  = (r1.yy - 2*r1.beta'*r1.Xy + r1.beta'*r1.XX*r1.beta)
    vcov = uPu/(N-K)*invsym(r1.XX)
    r1.beta
    vcov
}

struct results scalar update(pointer X, pointer y, real scalar block,
                             real scalar K, struct results scalar r) {
    struct results scalar a

    if (block==1) {
        a.Sigm1 = invsym(quadcross(*X,*X))
        a.beta  = a.Sigm1*(quadcross(*X,*y))
    }
    else {
        Sigma   = quadcross(*X,*X)
        Omega   = I(K)-r.Sigm1*Sigma*luinv(I(K)+r.Sigm1*Sigma)
        a.Sigm1 = Omega*r.Sigm1
        a.beta  = Omega*(r.beta+r.Sigm1*quadcross(*X,*y))
    }
    a.yy = r.yy+quadcross(*y,*y)
    a.Xy = r.Xy+quadcross(*X,*y)
    a.XX = r.XX+quadcross(*X,*X)
    return(a)
}
end
