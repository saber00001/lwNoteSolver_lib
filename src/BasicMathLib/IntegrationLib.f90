!numerical integeration is also denoted as cubature
!quadrature is more applied to one dimensional integral
!and cubature is referred to over more than one dimension
module IntegrationLib
use constants
use arrayOpsLib
use polynomial_
use SpecialFunctionLib
implicit none

    private
    
    public:: integrate
    
    public:: GaussLegendre
!----------------------------------------
    public:: integrateTrapezoid
    public:: integrateAdaptiveSimpson
    
    !--
    interface integrate
        procedure:: integrateAdaptiveSimpson_f1
        procedure:: integeratePolynomial
    end interface integrate
    !--

!---------------------------------------------------------    
!method for integrate one-dimensional function
    interface integrateAdaptiveSimpson
        procedure:: integrateAdaptiveSimpson_f1
    end interface integrateAdaptiveSimpson
    
    !not good
    interface integrateTrapezoid
        procedure::  integrateTrapezoid_f1
    end interface integrateTrapezoid
!---------------------------------------------------------
    
    
!---------------------------------------------------------
!Gaussian Quadrature Related
    
    
    

contains
    
    !-------------------------------------------------------------------------
    pure function integrateTrapezoid_f1(func,lo,up,epsilon) result(integral)
    procedure(absf1)::                      func
    real(rp),intent(in)::                   lo,up
    real(rp),intent(in),optional::          epsilon
    integer(ip)::                           n
    real(rp)::                              integral0,integral,epsrun,eps
        
        epsrun  = 1.d100
        n       = 101
        eps     = merge(epsilon,GlobalEps,present(epsilon))
        
        integral = trapInt()
        do while(epsrun>max(integral*eps,eps).and.n<1e9)
            integral0   = integral
            n           = n * 1.618d0     !golden
            integral    = trapInt()
            epsrun      = abs(integral-integral0)
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
    pure real(rp) function integrateAdaptiveSimpson_f1(f,lo,up,eps,maxLoop) result(integral)
    procedure(absf1)::                  f
    real(rp),intent(in)::               lo,up
    real(rp),intent(in),optional::      eps
    integer(ip),intent(in),optional::   maxLoop
        integral = integrateAdaptiveSimpson_f1_rcs(f,lo,f(lo),up,f(up), &
                                    simpson_f1(f,lo,f(lo),up,f(up)),    &
                                    merge(eps,Globaleps,present(eps)),  &
                                    merge(maxLoop,15_ip,present(maxLoop)))
    end function integrateAdaptiveSimpson_f1
    
    !--read flo and fup can reduce the computing, notice the <fmiddle>
    pure recursive real(rp) function integrateAdaptiveSimpson_f1_rcs(f,lo,flo,up,fup,s,eps,max) result(integral)
    procedure(absf1)::          f
    real(rp),intent(in)::       lo,flo,up,fup,s,eps
    integer(ip),intent(in)::    max
    real(rp)::                  mi,fmi,l,r
        mi  = (lo+up)/2.d0
        fmi = f(mi)
        l   = simpson_f1(f,lo,flo,mi,fmi)
        r   = simpson_f1(f,mi,fmi,up,fup)
        if(max==0.or.abs(l+r-s)<=15.d0*eps) then   !15.d0 comes from error analysis
            integral = l+r+(l+r-s)/15.d0
        else
            integral = integrateAdaptiveSimpson_f1_rcs(f,lo,flo,mi,fmi,l,eps/2.d0,max-1)+ &
                        integrateAdaptiveSimpson_f1_rcs(f,mi,fmi,up,fup,r,eps/2.d0,max-1)
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
    
    
    
!---------------------------------------------------
!refer to https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/legpts.m
!calculate \[ \int_{-1}^{1} f(x) dx \],(refer to wiki, \pi(x)=1 rather than \pi(x)=\frac{1}{2})
!asymptotic method unavailiable temporarily
    pure subroutine GaussLegendre(quadx,quadw)
    real(rp),dimension(:),intent(out)::         quadx,quadw
    integer(ip)::                               n
    real(rp)::                                  a,b,c,d
        !--
        n = size(quadx)
        !--
        if(n<6) then
        
            select case(n)
            case(1)
                quadx = 0._rp
                quadw = 2._rp
            case(2)
                c = 1._rp/sqrt(3._rp)
                !--
                quadx = [-c , c]
                quadw = [1._rp , 1._rp]
            case(3)
                c = sqrt(3._rp/5._rp)
                !--
                quadx = [-c , 0._rp , c]
                quadw = [5._rp/9._rp , 8._rp/9._rp , 5._rp/9._rp]
            case(4)
                a = 2._rp/7._rp * sqrt(6._rp/5._rp)
                b = sqrt(3._rp/7._rp - a)
                c = sqrt(3._rp/7._rp + a)
                d = sqrt(30._rp)/36._rp
                !--
                quadx = [-c , -b , b , c]            
                quadw = [0.5_rp-d , 0.5_rp+d , 0.5_rp+d , 0.5_rp-d]
            case(5)
                a = 2._rp * sqrt(10._rp/7._rp)
                b = 1._rp/3._rp * sqrt(5._rp - a)
                c = 1._rp/3._rp * sqrt(5._rp + a)
                d = 13._rp * sqrt(70._rp)
                !--
                quadx = [-c , -b , 0._rp , b , c]
                quadw = [ (322._rp-d)/900._rp , (322._rp+d)/900._rp , 128._rp/225._rp,  &
                            (322._rp+d)/900._rp , (322._rp-d)/900._rp ]
            end select
            
        elseif(n<80) then  !60,80,100: all ok
        
            call iterativeMethod(n)
            
        else
        
            call asymptoticMethod(n)
            
        endif
    
    contains
    
        !use newton's method to find the roots of Legendre Polynomials
        pure subroutine iterativeMethod(n)
        integer(ip),intent(in)::            n
        integer(ip)::                       k,iter
        real(rp),dimension(ishft(n+1,-1)):: x0,dx,pm1,pm2,ppm1,ppm2,p,pp

            !Asymptotic formula (Tricomi), only for positive x, as the initialization
            x0 = dfloat( [ ishft(n+1,-1) : 1 : -1 ] )
            x0 = pi * (4._rp*x0 - 1._rp) / (4._rp*n + 2._rp)
            x0 = (1._rp - (n-1._rp)/(8._rp*n**3) - &
                1._rp/(384._rp*n**4)*(39._rp-28._rp/sin(x0)**2) ) * cos(x0)
                
            !--
            !initialize
            pm1 = x0
            pm2 = 1._rp
            ppm1 = 1._rp
            ppm2 = 0._rp
            dx = 100._rp
            iter = 0
            !---
            do while(norm(dx,-1)>GlobalEps.and.iter<10)
                iter = iter + 1
                !
                do k=1,n-1
                    p = ((2._rp*k + 1._rp)*pm1*x0 - k*pm2)/(k+1._rp)
                    pm2 = pm1
                    pm1 = p
                    pp = ((2._rp*k + 1._rp)*(Pm2 + x0*PPm1) - k*PPm2 ) / (k+1._rp)
                    ppm2 = ppm1
                    ppm1 = pp
                enddo
            
                dx = - p/pp
                x0 = x0 + dx
                
                !--
                pm1 = x0
                pm2 = 1._rp
                ppm1 = 1._rp
                ppm2 = 0._rp
            enddo
            
            do k=1,n-1
                p = ((2._rp*k + 1._rp)*pm1*x0 - k*pm2)/(k+1._rp)
                pm2 = pm1
                pm1 = p
                pp = ((2._rp*k + 1._rp)*(Pm2 + x0*PPm1) - k*PPm2 ) / (k+1._rp)
                ppm2 = ppm1
                ppm1 = pp
            enddo
            
            quadx = [-x0(ishft(n+1,-1):1+ibits(n,0,1):-1),x0]
            !deriavtives
            quadw = [ pp(ishft(n+1,-1):1+ibits(n,0,1):-1),pp]
            !--refer to wiki(GausssQuadrature)
            quadw = 2._rp /( (1._rp-quadx**2) * quadw**2 )
            
        end subroutine iterativeMethod
        
        !--
        pure subroutine asymptoticMethod(n)
        integer(ip),intent(in)::        n
        integer(ip)::                   m
        real(rp),dimension(ishft(n+1,-1)):: jk
            
            !give up, too many, waiting for some day
            call disableprogram
            
            m = ishft(n+1,-1)
            jk = besseljn_roots(0,m)
            
        
        end subroutine asymptoticMethod
        
    end subroutine GaussLegendre
    
end module IntegrationLib