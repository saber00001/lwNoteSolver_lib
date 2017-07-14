module IntegrationLib
use constants
use polynomial_
implicit none

    private
    
    public:: integrate
    public:: integrateTrapezoid
    public:: integrateAdaptiveSimpson
    
!---------------------------------------------------------
    interface integrate
        procedure:: integrateAdaptiveSimpson_f1
        procedure:: integeratePolynomial
    end interface integrate

!--------------------------------------------------------
    interface integrateTrapezoid
        procedure::  integrateTrapezoid_f1
    end interface integrateTrapezoid
    
    interface integrateAdaptiveSimpson
        procedure:: integrateAdaptiveSimpson_f1
    end interface integrateAdaptiveSimpson
!---------------------------------------------------------

contains
    
    !-------------------------------------------------------------------------
    pure function integrateTrapezoid_f1(func,lo,up,tolin) result(integral)
    procedure(absf1)::                      func
    real(rp),intent(in)::                   lo,up
    real(rp),intent(in),optional::          tolin
    integer(ip)::                           n
    real(rp)::                              integral0,integral,er,tol
        
        er  = 1.d100
        n   = 101
        tol = merge(tolin,IntegralDefaultTolerance,present(tolin))
        
        integral = trapInt()
        do while(er>max(integral*tol,tol).and.n<1e9)
            integral0 = integral
            n = n * 1.618d0     !golden
            integral = trapInt()
            er = abs(integral-integral0)
        enddo
        
    contains
    
        pure function trapInt()
        real(rp)::      trapInt,x,y0,y1,dx
        integer(ip)::   i
        
            trapInt = zero
            dx  = (up-lo)/dfloat(n)
            x   = lo
            y1  = func(lo)
            
            do i=1,n
                x   = lo + i*dx
                y0  = y1
                y1  = func(x)
                trapInt = trapInt + 0.5d0*(y0+y1)*dx
            enddo
            
        end function trapInt
        
    end function integrateTrapezoid_f1
    
    
    !---------refer to wiki(adaptive simpson's method) and 
    !https: //www.zhihu.com/question/30141678
    pure function integrateAdaptiveSimpson_f1(f,lo,up,tolin) result(integral)
    procedure(absf1)::              f
    real(rp),intent(in)::           lo,up
    real(rp),intent(in),optional::  tolin
    real(rp)::                      integral,er,s,flo,fup
        er  = merge(tolin,IntegralDefaultTolerance,present(tolin))
        flo = f(lo)
        fup = f(up)
        s   = simpson_f1(f,lo,flo,up,fup)
        integral = integrateAdaptiveSimpson_f1_rcs(f,lo,flo,up,fup,er,s,Integralmaxloop)
    end function integrateAdaptiveSimpson_f1
    
    !--
    pure recursive function integrateAdaptiveSimpson_f1_rcs(f,lo,flo,up,fup,er,s,max) result(integral)
    procedure(absf1)::          f
    real(rp),intent(in)::       lo,flo,up,fup,er,s
    integer(ip),intent(in)::    max
    real(rp)::                  integral,mi,fmi,l,r
        mi = (lo+up)/2.d0
        fmi= f(mi)
        l = simpson_f1(f,lo,flo,mi,fmi)
        r = simpson_f1(f,mi,fmi,up,fup)
        if(max==0.or.abs(l+r-s)<=15.d0*er) then   !15.d0 comes from error analysis
            integral = l+r+(l+r-s)/15.d0
        else
            integral = integrateAdaptiveSimpson_f1_rcs(f,lo,flo,mi,fmi,er/2.d0,l,max-1)+ &
                        integrateAdaptiveSimpson_f1_rcs(f,mi,fmi,up,fup,er/2.d0,r,max-1)
        endif
    end function integrateAdaptiveSimpson_f1_rcs

    !--
    pure function simpson_f1(f,lo,flo,up,fup) result(r)
    procedure(absf1)::          f
    real(rp),intent(in)::       lo,up,flo,fup
    real(rp)::                  r,mi
        mi = lo + (up-lo)/2.d0
        r = ( flo + 4.d0*f(mi) + fup ) * (up-lo) / 6.d0
    end function simpson_f1
    
    !-----------------------------------------------
    pure real(rp) function integeratePolynomial(p,lo,up) result(r)
    class(polynomial),intent(in)::  p
    real(rp),intent(in)::           lo,up
        r = p%integral(lo,up)
    end function integeratePolynomial
    
    
end module IntegrationLib