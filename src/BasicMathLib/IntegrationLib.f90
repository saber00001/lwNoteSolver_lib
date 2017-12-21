!numerical integeration is also denoted as cubature
!quadrature is more applied to one dimensional integral
!and cubature is referred to over more than one dimension
module IntegrationLib
use constants
use arrayOpsLib
use SpecialFunctionLib
use fftWrapperLib
use LAwrapperLib
use odelib
implicit none

    private
    
    public:: integrate
    !--
    public:: GaussQuadrature
    public:: GaussLegendre
    public:: GaussHermite
    public:: ClenshawCurtis
    public:: ClenshawCurtisNative
    !--
    public:: integrateTrapezoid
    public:: integrateAdaptiveSimpson
    
    
    
!------------------------------------------------
    interface integrate
        procedure:: integrateAdaptiveSimpson_f1
        procedure:: integrateGaussRule
    end interface integrate
    !--

!--------------------------------------------------
    interface GaussHermite
        procedure:: GaussHermitePhys
    end interface

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
    
    
    
    !-------------------------------------------------
    pure real(rp) function integrateGaussRule(f,ruleType,Order,normalCoef) result(r)
    procedure(absf1)::              f
    character(*),intent(in)::       ruleType
    integer(ip),intent(in)::        Order
    real(rp),optional,intent(in)::  normalCoef
    real(rp),dimension(ishft(Order,-1)+1):: x,w
    
        call GaussQuadrature(ruleType,x,w)
        r = merge(normalCoef,1._rp,present(normalCoef)) * sum(f(x)*w)
        
    end function integrateGaussRule
    
    

!---------------------------------------------------
    pure subroutine GaussQuadrature(ruletype,quadx,quadw,quadv)
    character(*),intent(in)::                       ruleType
    real(rp),dimension(:),intent(out)::             quadx,quadw
    real(rp),dimension(:),intent(out),optional::    quadv
    
        select case(adjustl(ruleType))
        case('Legendre','legendre')
            call GaussLegendre(quadx,quadw)
        case('Hermite','hermite')
            if(.not.present(quadv)) then
                call GaussHermite(quadx,quadw)
            else
                call GaussHermite(quadx,quadw,quadv)
            end if
        case default
            quadx = nanrp; quadw = nanrp
        end select
        
    end subroutine GaussQuadrature
    
    
!refer to https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/legpts.m
!maybe a better choice https://github.com/Pazus/Legendre-Gauss-Quadrature/blob/master/legzo.m
!calculate \[ \int_{-1}^{1} f(x) dx \],(refer to wiki, \pi(x)=1 rather than \pi(x)=\frac{1}{2})
!asymptotic method unavailiable temporarily
!n points reach 2n-1 order accuracy | N+1 points reach 2N+1 order accuracy
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
            
        elseif(n<1e6) then  !60,80,100: all ok. temperarialy relax the restriction
        
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
            call disableProgram
            
            m = ishft(n+1,-1)
            jk = besseljn_roots(0,m)
            
        
        end subroutine asymptoticMethod
        
    end subroutine GaussLegendre
!------------------------------------------------------------------------------------
!refer to https://github.com/chebfun/chebfun/blob/development/hermpts.m
!calculate \[ \int_{-inf}^{inf}exp(-x^2)f(x) dx \],(refer to wiki)
!hermite asymptotic method unavailiable temporarily
!GaussHermite returns hermite points and gauss-hermite quadrature weights.
!by default these are roots of the 'physicist'-type hermite polynomials, which are orthogonal with respect to the weight exp(-x.^2).
!------------------------------------------------------------------------------------
    pure subroutine GaussHermitePhys(quadx,quadw,quadv)
    real(rp),dimension(:),intent(out)::                   quadx,quadw
    real(rp),dimension(:),intent(out),optional::          quadv  
    integer(ip)::                                           n
    logical(lp)::                                           ve
        
        ve = present(quadv)
        n=size(quadx)
        if(n==0) then
            quadx = 0._rp
            quadw = 0._rp
            if(ve) quadv = 0._rp
        elseif(n==1) then
            quadx = 0._rp 
            quadw = spi 
            if(ve) quadv = 1._rp
        else if(n<21) then
            call GolubWelsch_eigenvalueMethod(n)
        else if(n<200) then
            call hermpts_rec(n)
        else
            !normally,we won't use this subroutine which is not accomplished
            call hermpts_asy(n)
        end if
        
        !Normalise so that sum(w) = spi
        quadw = (spi/sum(quadw))*quadw     
        
        if(ve) call CaculateQuadv(n)
            
    contains
    
        pure subroutine GolubWelsch_eigenvalueMethod(n)
        integer(ip),intent(in)::                        n
        integer(ip)::                                   i,info
        integer(ip),dimension(n)::                      index_
        integer(ip),dimension(ishft(n,-1))::            ii
        real(rp)::                                      vmid
        real(rp),dimension(n-1)::                       beta,D
        real(rp),dimension(n)::                         Value1,D_G
        real(rp),dimension(n,n)::                       T,Z_real
        real(rp),dimension(ishft(n,-1))::               x_half,w_half,v_half
        complex(rp),dimension(n,n)::                    Z
        
            call GolubWelschParameter(n,D_G,Z_real,index_)
            quadx = D_G
            quadw = spi*Z_real(1,index_)**2
            !Enforce symmetry
            ii = [(i,i=1,ishft(n,-1))]                                      
            x_half = quadx(ii)
            w_half = quadw(ii)       
            if (ibits(n,0,1)) then                                          
                quadx = [x_half(:),0._rp,-x_half(ishft(n,-1):1:-1)]
                quadw = [w_half(:),spi-sum(2._rp*w_half(:)),w_half(ishft(n,-1):1:-1)]
            else
                quadx = [x_half(:),-x_half(ishft(n,-1):1:-1)]
                quadw = [w_half(:),w_half(ishft(n,-1):1:-1)]
            end if    
            
        end subroutine GolubWelsch_eigenvalueMethod
        !--------------------------------------------------------------------------
        pure subroutine GolubWelschParameter(n,D_G,Z_real,index_)
        integer(ip),intent(in)::                            n
        real(rp),dimension(:),intent(out)::                 D_G
        real(rp),dimension(:,:),intent(out)::               Z_real
        integer(ip),dimension(:),intent(out)::              index_
        integer(ip)::                                       i
        real(rp),dimension(n-1)::                           beta,D
        real(rp),dimension(n)::                             Value1
        real(rp),dimension(n,n)::                           T
        complex(rp),dimension(n,n)::                        Z
        
            beta = [(sqrt(0.5_rp*i),i=1,n-1)]
            T = diag(beta,1)+diag(beta,-1)                  
            Value1 = [(T(i,i),i=1,n)]
            D = [(T(i,i+1),i=1,n-1)]
            call eigenTriDiagonal(Value1,D,Z)
            Z_real = real(Z,kind=rp)
            D_G = Value1
            call sort(D_G,index_)
            
        end subroutine GolubWelschParameter

        !--------------------------------------------------------------------------
        pure subroutine hermpts_rec(n)
        integer(ip),intent(in)::                                n
        real(rp),dimension(:),allocatable::                     x1,w1,v1    
        integer(ip)::                                           m,kk,i,j  
        
            if(ibits(n,0,1)==1) then
                m=ishft(n,-1)
                allocate(x1(m+1),w1(m+1),v1(m+1))
            else
                m=ishft(n,-1)
                allocate(x1(m),w1(m),v1(m))
            end if  
            
            call Hermpts_recParameter(n,m,x1,w1,v1)
            
            if (ibits(n,0,1)==1) then                                   !fold out
                quadx = [-x1(m+1:1:-1),x1(2:m+1)]
                quadw = [w1(m+1:1:-1),w1(2:m+1)]
                quadw = (spi/sum(quadw))*quadw
            else
                quadx = [-x1(m:1:-1),x1(:)]
                quadw = [w1(m:1:-1),w1(:)]
                quadw = (spi/sum(quadw))*quadw
            end if
        
        end subroutine hermpts_rec
        !-------------------------------------------------------
        pure subroutine Hermpts_recParameter(n,m,x1,w1,v1)
        integer(ip),intent(in)::                                n
        integer(ip),intent(inout)::                             m
        real(rp),dimension(:),intent(out)::                     x1,w1,v1    
        real(rp),dimension(:),allocatable::                     x0,val,dval,dx
        integer(ip)::                                           kk,i,j,size_para  
        
            size_para=size(x1)
            allocate(x0(size_para),val(size_para),dval(size_para),dx(size_para))
            
            call HermiteInitialGuesses(m,n,x0)
            x0 = x0*sqrt(2._rp)
            do kk = 1,10
                call hermpoly_rec(n,x0,val,dval)  
                dx = val/dval
                do i = 1,size(dx)
                    if(isnan(dx(i))) dx(i)=0._rp 
                end do
                x0 = x0 - dx
                if(norm(dx,-1)<sqrt(GlobalEps)) exit
            end do
        
            x1 = x0/sqrt(2._rp)
            w1 = (exp(-x1**2)/dval**2)                                  !quadrature weights
            v1 = exp(-x1**2/2._rp)/dval                                 !Barycentric weights
        end subroutine Hermpts_recParameter
        !-------------------------------------------------------
        pure subroutine HermiteInitialGuesses(m,n,x)
        integer(ip),intent(in)::                            m,n
        real(rp),dimension(:),intent(out)::                 x
        integer::                                           i,k
        real(rp)::                                          a,nu,p,tin
        real(rp),dimension(m)::                             airyrts,x_init,x_init_airy,Tnk0,rhs,val,dval,dTnk0,tnk,x_init_sin 
        real(rp),dimension(10)::                            airyrts_exact
            ! hermiteintitialguesses(n), Initial guesses for Hermite zeros.
            if(ibits(n,0,1)==1) then
                a =  0.5_rp
            else
                a = -0.5_rp
            end if  
            nu = 4*m + 2*a + 2   
            do i =1,m
                airyrts(i) = -T(3._rp/8._rp*pi*(4._rp*i-1._rp))
            end do
        
            airyrts_exact = [-2.338107410459762_rp,-4.087949444130970_rp,-5.520559828095555_rp,&               ! Exact Airy roots.
                             -6.786708090071765_rp,-7.944133587120863_rp,-9.022650853340979_rp,&
                            -10.040174341558084_rp,-11.008524303733260_rp,-11.936015563236262_rp,&
                            -12.828776752865757_rp]
            
            airyrts(1:10) = airyrts_exact(1:10)                                         ! correct first 10.
        
            x_init = real(sqrt(cmplx(nu + 2**(2._rp/3._rp)*airyrts*nu**(1._rp/3._rp) + &
            1._rp/5._rp*2**(4._rp/3._rp)*airyrts**2*nu**(-1._rp/3._rp) + &
            (11._rp/35._rp-a**2-12._rp/175._rp*airyrts**3)/nu + &
            (16._rp/1575._rp*airyrts+92._rp/7875._rp*airyrts**4)*2**(2._rp/3._rp)*nu**(-5._rp/3._rp) - &
            (15152._rp/3031875._rp*airyrts**5+1088._rp/121275._rp*airyrts**2)*2**(1._rp/3._rp)*nu**(-7._rp/3._rp),kind=rp)),kind=rp)  
            x_init_airy(1:m) = x_init(m:1:-1)
     
            ! Tricomi initial guesses. Equation (2.1) in [1]. Originally in [2].
            ! These initial guesses are good near x = 0 . Note: zeros of besselj(+/-.5,x)
            ! are integer and half-integer multiples of pi.
            ! x_init_bess =  bess/sqrt(nu).*sqrt((1+ (bess.^2+2*(a^2-1))/3/nu^2) );
    
            Tnk0 = [(pi/2,i=1,m)]
            nu = 4*m+2*a+2
            rhs = [((4*m-4*i+3)/nu*pi,i=1,m)]
            do k = 1,7
                val = Tnk0 - sin(Tnk0) - rhs 
                dval = 1 - cos(Tnk0)
                dTnk0 = val/dval 
                Tnk0 = Tnk0 - dTnk0 
            end do
      
            tnk = cos(Tnk0/2._rp)**2
            x_init_sin = sqrt(nu*tnk - (5._rp/(4._rp*(1._rp-tnk)**2) - 1._rp/(1._rp-tnk)-1._rp+3._rp*a**2)/3._rp/nu)
            !Patch together
            p = 0.4985_rp+GlobalEps
            x_init = [x_init_sin(1:floor(p*n)),x_init_airy(ceiling(p*n):m)]                   
            if (ibits(n,0,1)==1) then
                x = [0._rp,x_init]                                                             
                x = x(1:m+1)
            else
                x = x_init(1:m)
            end if
            
        end subroutine HermiteInitialGuesses   
        !-------------------------------------------------------
        pure real(rp) function T(x) result(y)
        real(rp),intent(in)::               x
        
            y = x**(2._rp/3._rp)*(1._rp+5._rp/48._rp*x**(-2)-5._rp/36._rp*x**(-4)+(77125._rp/82944._rp)*x**(-6)-108056875._rp/6967296._rp*x**(-8)+162375596875._rp/334430208._rp*x**(-10))
        
        end function T
        !--------------------------------------------------------------------------
        pure subroutine hermpoly_rec(n,x0,val,dval)
        integer(ip),intent(in)::                        n
        real(rp),dimension(:),intent(in)::              x0
        real(rp),dimension(:),intent(out)::             val,dval
        integer(ip)::                                   k
        real(rp),dimension(size(x0))::                  Hold,H,Hnew      
            !HERMPOLY_rec evaluation of scaled Hermite poly using recurrence
            !evaluate:
            Hold = exp(-x0**2/4._rp) 
            H = x0*exp(-x0**2/4._rp)
            do k = 1,n-1
                Hnew = (x0*H/sqrt((k+1)*1._rp) - Hold/sqrt((1+1._rp/k)*1._rp))
                Hold = H 
                H = Hnew
            end do
            !evaluate derivative:
            val = Hnew
            dval = (-x0*Hnew + n**(1._rp/2._rp)*Hold)
        end subroutine hermpoly_rec
        !--------------------------------------------------------------------------
        pure subroutine hermpts_asy(n)
        integer(ip),intent(in)::                        n
        
        !HERMPTS_ASY, fast algorithm for computing Gauss-Hermite nodes and weights 
        !using Newton's method with polynomial evaluation via asymptotic expansions.   
        !x = Gauss-Hermite nodes, w = quad weights, v = bary weights.
        call disableProgram    
        end subroutine hermpts_asy
        !--------------------------------------------------------------------------
        pure subroutine CaculateQuadv(n)
        integer(ip),intent(in)::                        n
        real(rp),dimension(n)::                         D_G
        real(rp),dimension(n,n)::                       Z_real
        integer(ip),dimension(n)::                      index_
        integer(ip)::                                   i
        real(rp)::                                      vmid
        real(rp),dimension(ishft(n,-1))::               v_half
        integer(ip),dimension(ishft(n,-1))::            ii
        real(rp),dimension(n)::                         ders
        real(rp),dimension(:),allocatable::             x1,w1,v1      
        integer(ip)::                                   m
      
            if(n==0) then
                quadv = 0._rp
                
            elseif(n==1) then 
                quadv = 1._rp
            
            else if(n<21) then
                call GolubWelschParameter(n,D_G,Z_real,index_)
                quadv = abs(Z_real(1,index_))
                quadv = quadv/maxval(quadv)
                quadv(2:n:2) = -quadv(2:n:2)
                !Enforce symmetry
                ii = [(i,i=1,ishft(n,-1))]                              
                vmid   = quadv(ishft(n,-1)+1)
                v_half = quadv(ii)
                if (ibits(n,0,1)) then                                          
                    quadv = [v_half(:),vmid,v_half(ishft(n,-1):1:-1)]
                else
                    quadv = [v_half(:),-v_half(ishft(n,-1):1:-1)]
                end if     
            
            else if(n<200) then
                m=ishft(n,-1)
                if(ibits(n,0,1)==1) then
                    allocate(x1(m+1),w1(m+1),v1(m+1))
                else
                    m=ishft(n,-1)
                    allocate(x1(m),w1(m),v1(m))
                end if  
                call Hermpts_recParameter(n,m,x1,w1,v1)
                if (ibits(n,0,1)==1) then                                       !fold out
                    quadv = [v1(m+1:1:-1),v1(2:m+1)]
                    quadv = quadv/maxval(abs(quadv))
                else
                    quadv = [v1(m:1:-1),-v1(:)]
                    quadv = quadv/maxval(abs(quadv))
                end if   
            else
                !normally,we won't use this subroutine which is not accomplished
                call hermpts_asy(n)
            end if
        end subroutine  CaculateQuadv
        
    end subroutine GaussHermitePhys
    !----------------------------------------------------------------------------------------------------------------------------------
    pure subroutine GaussHermiteProb(quadx,quadw,quadv)
    real(rp),dimension(:),intent(out)::                     quadx,quadw
    real(rp),dimension(:),intent(out),optional::            quadv
   
        if(present(quadv)) then
            call GaussHermitePhys(quadx,quadw,quadv)
            quadx = quadx*sqrt(2._rp)
            quadw = quadw*sqrt(2._rp)
        else
            call GaussHermitePhys(quadx,quadw)
            quadx = quadx*sqrt(2._rp)
            quadw = quadw*sqrt(2._rp)
        end if
        
    end subroutine GaussHermiteProb
    
!----------------------------------------------------------------------------------------------------------------------------------
    pure subroutine GaussHermitePhys_GLR(quadx,quadw,quadv)
    real(rp),dimension(:),intent(out)::                     quadx,quadw
    real(rp),dimension(:),intent(out),optional::            quadv  
    integer(ip)::                                           n
    logical(lp)::                                           vexist
        vexist=present(quadv)
        n=size(quadx)
        call GlaserLiuRokhlinAlgorithm(n)
        if(vexist) call GaussHermitePhys_GLR_Caculate_Quadv(n)
        
    contains
        
        pure subroutine GlaserLiuRokhlinAlgorithm(n)
        integer(ip),intent(in)::                                n
        real(rp),dimension(n)::                                 ders
        
            !'GLR' uses Glaser-Liu-Rokhlin fast algorithm which is much faster for large N
            call alg0_Herm(n,quadx,ders)
            quadw = (2._rp*exp(-quadx**2)/ders**2)                  ! Quadrature weights
            
        end subroutine GlaserLiuRokhlinAlgorithm
        
        !-------------------------------------------------------
        pure subroutine alg0_Herm(n,x,ders)
        integer(ip),intent(in)::                                n
        real(rp),dimension(:),intent(out)::                     x,ders    
            ! Compute coefficients of H_m(0), H_m'(0), m = 0,..,N
            call first_root(n,x,ders)                           !find first root
            call alg1_Herm(n,x,ders)                            !compute roots and derivatives
        
        end subroutine alg0_Herm
        !-------------------------------------------------------
        pure subroutine first_root(n,x,ders)
        integer(ip),intent(in)::                                n
        real(rp),dimension(:),intent(out)::                     x,ders
        real(rp)::                                              Hm2, Hm1, Hpm2, Hpm1,H,Hp
        integer(ip)::                                           k                    
        
            Hm2 = 0._rp
            Hm1 = pi**(-1._rp/4._rp)
            Hpm2 = 0._rp
            Hpm1 = 0._rp       
            do k = 0,n-1
                H = -sqrt(1._rp*k/(k+1))*Hm2
                Hp = sqrt(2._rp/(k+1))*Hm1-sqrt(1._rp*k/(k+1))*Hpm2
                Hm2 = Hm1
                Hm1 = H
                Hpm2 = Hpm1
                Hpm1 = Hp
            end do
            x=0._rp
            ders=0._rp                      
            if (ibits(n,0,1)) then
                !zero is a root:
                x(ishft((n-1),-1)) = 0._rp                          !x((n-1)/2) = 0 
                ders(ishft((n+1),-1)) = Hp                          !ders((n+1)/2) = Hp        
            else
                !find first root
                call alg2_Herm(H,n,x(ishft(n,-1)+1),ders(ishft(n,-1)+1))
            end if    
        
        end subroutine first_root
        !-------------------------------------------------------
        pure subroutine alg2_Herm(H,n,x1,d1) 
        real(rp),intent(in)::                               H
        integer(ip),intent(in)::                            n
        real(rp),intent(out)::                              x1,d1
        integer(ip)::                                       k,i,l
        integer(ip),parameter::                             m = 30      !number of terms in Taylor expansion
        real(rp)::                                          m1,c,step
        real(rp),dimension(m)::                             z,z1,z2
        real(rp),dimension(m+1)::                           u,up,x1k
            ! find the first root (note H_n'(0) = 0)
            ! advance ODE via Runge-Kutta for initial approx
            x1 = odeRK2(0._rp,-pi/2,0._rp,n)  
            !scaling
            m1 = 1/x1
            
            !initialise
            u = 0._rp
            up =0._rp
            
            !recurrence relation for Legendre polynomials
            u(1) = H
            u(3) = -0.5_rp*(2._rp*n+1)*u(1)/m1**2                       !attention
            up(1) = 0._rp 
            up(2) = 2._rp*u(3)*m1
            do k = 2,m-2,2
                u(k+3) = (-(2._rp*n+1)*u(k+1)/m1**2 + u(k-1)/m1**4)/((k+1)*(k+2))
                up(k+2) = (k+2)*u(k+3)*m1
            end do
            
            !flip for more accuracy in inner product calculation
            u = [(u(i) ,i = m+1,1,-1)]
            up =[(up(i),i = m+1,1,-1)]
            z = 0._rp                                                  
            z1 = m1*x1+z
            z2 = cumprod(z1)                                          
            x1k = [m1,(z2(i),i=1,m)]
            step = maxrp                                                
            l = 0
            
            !Newton iteration
            do while ((abs(step)>GlobalEps).and.(l<10))
                l = l + 1
                step = (u.ip.x1k)/(up.ip.x1k)
                x1 = x1 - step
                !powers of h (This is the fastest way!)
                z1 = m1*x1+z
                z2 = cumprod(z1)
                x1k = [1._rp,(z2(i),i=1,m)]
                x1k = [(x1k(i),i=m+1,1,-1)]
            end do  
            
            !Update derivative
            d1 = (up.ip.x1k)
        end subroutine alg2_Herm
        !-------------------------------------------------------
        pure subroutine alg1_Herm(n,roots,ders)
        integer(ip),intent(in)::                            n
        real(rp),dimension(:),intent(out)::                 roots,ders
        integer(ip)::                                       s,N1,j,l,i,k
        integer(ip),parameter::                             m = 30  !number of terms in Taylor expansion 
        real(rp)::                                          x,h,m1,c1,c2,c3,step 
        real(rp),dimension(m)::                             z,z1,z2
        real(rp),dimension(m+1)::                           hh1,u,up,hh
        
            s=ibits(n,0,1)
            N1=ishft(n,-1)
            !initialise
            hh1= 1._rp
            u  = 0._rp
            up = 0._rp
        
            do j =(N1+1),(n-1)
                !previous root
                x = roots(j) 
                !initial approx
                h = odeRK2(pi/2._rp,-pi/2._rp,x,n)-x
                !scaling
                m1 = 1._rp/h            
                !recurrence relation for Hermite polynomials
                c1 = -(2._rp*n+1-x**2)/(m1**2) 
                c2 = 2._rp*x/(m1**3)
                c3 = 1._rp/(m1**4)
                u(1) = 0._rp
                u(2) = ders(j)/m1 
                u(3) = 0.5_rp*c1*u(1)
                u(4) = (c1*u(2)+c2*u(1))/6._rp
                up(1) = u(2)
                up(2) = 2._rp*u(3)*m1
                up(3) = 3._rp*u(4)*m1
                up(m+1) = 0._rp
                do k = 2,m-2
                    u(k+3) = (c1*u(k+1) + c2*u(k) + c3*u(k-1))/((k+1)*(k+2))
                    up(k+2) = (k+2)*u(k+3)*m1
                end do
                ! flip for more accuracy in inner product calculation
                u = [(u(i) ,i = m+1,1,-1)]
                up =[(up(i),i = m+1,1,-1)]
                !Newton iteration
                hh = hh1
                hh(m+1) = m1
                step = maxrp   
                l = 0
                z = 0._rp
                do while ((abs(step)>GlobalEps).and.(l<10))
                    l = l + 1              
                    step = (u.ip.hh)/(up.ip.hh)
                    h = h - step
                    !powers of h (This is the fastest way!)
                    z1 = m1*h+z
                    z2 = cumprod(z1)
                    hh = [m1,(z2(i),i=1,m)]           
                    !flip for more accuracy in inner product calculation
                    hh(1:m+1) = hh(m+1:1:-1)                
                end do
                roots(j+1) = x + h
                !ders(j+1) =  dot_product(up,hh)
                ders(j+1) = (up.ip.hh)
            end do   
            !nodes are symmetric
            roots(1:N1+s) = -roots(n:N1+1:-1)
            ders(1:N1+s) = ders(n:N1+1:-1)        
        end subroutine alg1_Herm
        
        pure subroutine GaussHermitePhys_GLR_Caculate_Quadv(n)
        integer(ip),intent(in)::                                n
        integer(ip),dimension(ishft(n,-1))::                    ii          
        integer(ip)::                                           i
        real(rp),dimension(n)::                                 ders
             
                call alg0_Herm(n,quadx,ders)
                quadv = exp(-quadx**2/2._rp)/ders                       ! Barycentric weights
                quadv = quadv/maxval(abs(quadv))                        ! Normalize
                if (.not.ibits(n,0,1)) then
                   ii = [(i,i=ishft(n,-1)+1,n)]
                   quadv(ii) = -quadv(ii) 
                end if   
        end subroutine GaussHermitePhys_GLR_Caculate_Quadv
        
    end subroutine GaussHermitePhys_GLR
    
!--------------------------------------------------------------------------------------------------
    !https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/%40chebtech2/chebpts.m
    !https://github.com/chebfun/chebfun/blob/34f92d12ca51003f5c0033bbeb4ff57ac9c84e78/%40chebtech2/quadwts.m
    !also refer to <lwReference>/<Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules>
    !calculate $\int_{-1}^{1} f(x) dx$
    !n points reach n-1 order accuracy | N+1 points reach N orderaccuracy
    subroutine ClenshawCurtis(quadx,quadw)
    real(rp),dimension(:),intent(out)::         quadx,quadw
    integer(ip)::                               i,n,nm1,nhc
    complex(rp),dimension(size(quadx)-1)::      c
    
        n = size(quadx)
        nm1 = n - 1
        nhc = ishft(size(quadx)+1,-1)   !n_half_ceil
        
        !--
        if(n==1) then
            quadx(1) = 0._rp
            quadw(1) = 2._rp
            return
        endif
        
        !--
        forall(i=0:nm1) quadx(i+1) = - cos( i * pi / real(nm1,kind=rp) )
        
        !--Exact integrals of $\int_{-1}^{1} T_k (x) dx (k \in even number of [0,n-1])$
        !Warning...
        !this array constructor leads to stack overflow, if necessary, modify it.
        c(1:nhc) = cmplx( 2._rp/[1._rp , 1._rp - real([2:n-1:2]**2,kind=rp)] , kind=rp )
        c(1:nm1) = [c(1:nhc),c(ishft(n,-1):2:-1)]
        
        call ifft(c)
        
        quadw(1:nm1) = real(c,kind=rp)
        quadw(1) = quadw(1)/2._rp
        quadw(n) = quadw(1)
    
    end subroutine ClenshawCurtis
    
    !a quite slower method, as a reference
    !refer to http://people.sc.fsu.edu/~jburkardt/f_src/f_src.html %sparse_grid_cc
    pure subroutine ClenshawCurtisNative(quadx,quadw)
    real(rp),dimension(:),intent(out)::         quadx,quadw
    integer(ip)::                               n,nm1,i,j
    real(rp)::                                  b,theta(size(quadx))
          
        n = size(quadx)
        nm1 = n-1
        forall(i=1:n) theta(i) = (i-1) * pi / nm1
        
        !--
        if(n==1) then
            quadx(1) = 0._rp
            quadw(1) = 2._rp
            return
        endif

        !--
        forall(i=1:n) quadx(i) = - cos(theta(i))
        
        !--
        do i = 1 , n
            quadw(i) = 1._rp
            do j = 1 , nm1 / 2
                b = merge(1._rp , 2._rp , 2*j==nm1)
                quadw(i) = quadw(i) - b * cos(2._rp * j * theta(i)) / (4._rp * j**2 - 1._rp)
            end do
        end do

        quadw = quadw / nm1
        quadw(2:nm1) = 2._rp * quadw(2:nm1)
    
    end subroutine ClenshawCurtisNative
    
end module IntegrationLib