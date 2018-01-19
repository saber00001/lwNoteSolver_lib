module polynomial_
use constants
use arrayOpsLib
use SpecialFunctionLib
implicit none
    
    private
    public:: polynomial
    
    !--add special integrate method for polynomial
    public:: integrate
    
    !some related polynomials and function
    !--zero
    public:: zeroPolynomial
    !--Legendre, use recursive method rather than explict expression to avoid failure of large n for binomialCoef/factorial method
    public:: LegendrePolynomial
    public:: LegendrePolynomialSet
    public:: normalLegendrePolynomial
    public:: normalLegendrePolynomialSet
    !--Hermite, use recursive method rather than explict expression
    public:: HermitePolynomial
    public:: HermitePolynomialSet
    public:: normalHermitePolynomial
    public:: normalHermitePolynomialSet
    !--chebyshev
    public:: ChebyshevPolynomialT
    public:: ChebyshevPolynomialTset
    public:: ChebyshevPolynomialT_Clenshaw
    !--
    type:: polynomial
    
        private
        
        real(rp),allocatable,dimension(:):: coef_
        
    contains
    
        !--
        generic::               init    => init_degree,init_ar,init_ply
        procedure,private::     init_degree
        procedure,private::     init_ar
        procedure,private::     init_ply
        !--
        procedure::             degree
        procedure::             scoef
        procedure::             coefadd
        procedure::             contract
        procedure::             funcval
        procedure::             integral
        procedure::             differential
        !procedure::             realRoot
        !--
        generic::               coef => coef_i,coef_ptr
        procedure,private::     coef_i
        procedure,private::     coef_ptr
        !--
        generic::           assignment(=)   => paEq
        generic::           operator(+)     => psPlus,spPlus,ppPlus
        generic::           operator(-)     => ppMinus,negativePoly
        generic::           operator(*)     => ppMultiply,spmultiply,psMultiply
        generic::           operator(/)     => psDivide
        generic::           operator(==)    => ppjdEq
        !a strange override
        !https: //software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/713637
        procedure,pass(lhs),private::   paEq
        procedure,pass(lhs),private::   psPlus
        procedure,pass(rhs),private::   spPlus
        procedure,pass(lhs),private::   ppPlus
        procedure,pass(rhs),private::   negativePoly
        procedure,pass(lhs),private::   ppMinus
        procedure,pass(lhs),private::   ppMultiply
        procedure,pass(lhs),private::   psMultiply
        procedure,pass(rhs),private::   spMultiply
        procedure,pass(lhs),private::   psDivide
        procedure,pass(lhs),private::   ppjdEq

    end type polynomial

    
    
!---------------------------------------------------------------------
    interface integrate
        procedure:: integratePolynomial
    end interface integrate   
    !---------------------------
    interface HermitePolynomial
        procedure:: HermitePhysicistsPolynomial
    end interface HermitePolynomial
    !---------------------------
    interface HermitePolynomialSet
        procedure:: HermitePhysicistsPolynomialSet
    end interface HermitePolynomialSet
    !---------------------------
    interface normalHermitePolynomial
        procedure:: normalHermitePhysicistsPolynomial
    end interface normalHermitePolynomial
    !---------------------------
    interface normalHermitePolynomialSet
        procedure:: normalHermitePhysicistsPolynomialSet
    end interface normalHermitePolynomialSet
    !---------------------------
contains


    !--
    elemental subroutine init_degree(this,n)
    class(polynomial),intent(out)::         this
    integer(ip),intent(in)::                n
        allocate(this%coef_(0:n))
        this%coef_  = zero
    end subroutine init_degree
    
    !--
    pure subroutine init_ar(this,ar)
    class(polynomial),intent(out)::         this
    real(rp),dimension(0:),intent(in)::     ar
        allocate(this%coef_,source = ar)
    end subroutine init_ar
    
    !--
    elemental subroutine init_ply(this,that)
    class(polynomial),intent(out)::     this
    class(polynomial),intent(in)::      that
        allocate(this%coef_,source = that%coef_)
    end subroutine init_ply

    
    !---
    function coef_ptr(this)
    class(polynomial),target,intent(in)::this
    real(rp),dimension(:),pointer::     coef_ptr
        coef_ptr => this%coef_
    end function coef_ptr
    
    !---
    elemental real(rp) function coef_i(this,i)
    class(polynomial),intent(in)::  this
    integer(ip),intent(in)::        i
        coef_i = this%coef_(i)
    end function coef_i
    
    !---
    elemental subroutine scoef(this,i,coef)
    class(polynomial),intent(inout)::  this
    integer(ip),intent(in)::           i
    real(rp),intent(in)::              coef
        this%coef_(i) = coef
    end subroutine scoef
    
    !--
    elemental subroutine coefadd(this,i,v)
    class(polynomial),intent(inout)::   this
    integer(ip),intent(in)::            i
    real(rp),intent(in)::               v
        this%coef_(i) = this%coef_(i) + v
    end subroutine coefadd
    
    !---
    elemental integer(ip) function degree(this)
    class(polynomial),intent(in)::      this
    type(polynomial)::                  cthis
        cthis = this%contract()
        degree = ubound(cthis%coef_,dim=1)
    end function degree
    
    !--
    elemental type(polynomial) function contract(this) result(cp)
    class(polynomial),intent(in)::      this
    integer(ip)::                       i,n
        n = ubound(this%coef_,dim=1)
        do i=n,1,-1
            if(this%coef_(i)==zero) then
                n = n - 1
            else
                exit
            endif
        enddo
        allocate(cp%coef_(0:n),source=this%coef_(0:n))
    end function contract
    
    !--
    elemental real(rp) function funcval(this,x) result(y)
    class(polynomial),intent(in)::  this
    real(rp),intent(in)::           x
        y = polyval(this%coef_,x)
    end function funcval
    
    !--
    elemental real(rp) function integral(this,lo,up)
    class(polynomial),intent(in)::  this
    real(rp),intent(in)::           lo,up
    real(rp)::                      ui,li,coef
    integer(ip)::                   i
        ui = zero; li = zero
        do i=0,this%degree()
            coef = this%coef(i) / dfloat(i+1)
            ui = ui + up**(i+1) * coef
            li = li + lo**(i+1) * coef
        enddo
        integral = ui - li  !less minus better
    end function integral    
    
    !--
    elemental type(polynomial) function differential(this) result(d)
    class(polynomial),intent(in)::  this
    integer(ip)::                   i,n
        n = this%degree()
        if(n==0) then
            d = [0._rp]
        else
            call d%init( n - 1 )
            do i=0,n-1
                d%coef_(i) = this%coef_(i+1) * real(i+1,kind=rp)
            enddo
        endif
    end function differential
    
    !--Constructing
    pure subroutine realRoot(this,root)
    class(polynomial),intent(in)::                  this
    real(rp),dimension(:),allocatable,intent(out):: root
    real(rp),dimension(:),allocatable::             extremum
    integer(ip)::                                   n
    type(polynomial)::                              w

        root = 0._rp
        call disableprogram    
        !--
        w = this%contract()
        n = ubound(w%coef_,dim=1)
        select case(n)
        case(0)
            return
        case(1)
            call p1zero(w%coef_,root)
        case(2)
            call p2zero(w%coef_,root)
        case default
            call extremumSearch(w,root)
        endselect
        
    contains
    
        pure subroutine p1zero(a,root)
        real(rp),dimension(0:),intent(in)::             a
        real(rp),dimension(:),allocatable,intent(out):: root
            allocate(root,source=[-a(0)/a(1)])
        end subroutine p1zero
        !--
        pure subroutine p2zero(a,root)
        real(rp),dimension(0:),intent(in)::             a
        real(rp),dimension(:),allocatable,intent(out):: root
        real(rp)::                                      d
            d = a(1)**2 - 4._rp *a(0)*a(2)
            if(d<0._rp) then
                return
            elseif(d==0._rp) then
                allocate(root,source=[-a(1)/a(2)/2._Rp])
            else
                allocate(root,source=[(-a(1)-sqrt(d))/a(2)/2._rp,(-a(1)+sqrt(d))/a(2)/2._rp])
            endif
        end subroutine p2zero
        !--hard to compute
        pure recursive subroutine extremumSearch(w,root)
        class(polynomial),intent(in)::                  w
        real(rp),dimension(:),allocatable,intent(out):: root
        type(polynomial)::                              dif
        real(rp),dimension(:),allocatable::             dz
        integer(ip)::                                   n
            if(ubound(w%coef_,dim=1)==2) then
                call p2zero(w%coef_,root)
            else
                dif = w%differential()
                call extremumSearch(dif,dz)
                if(allocated(dz)) then
                    n = size(dz)
                else
                    return
                endif
            endif
        end subroutine extremumSearch
        
    end subroutine realRoot
    
!--------operator
    pure subroutine paEq(lhs,rhs)
    class(polynomial),intent(out)::     lhs
    real(rp),dimension(0:),intent(in):: rhs
        allocate(lhs%coef_,source=rhs)
    end subroutine paEq
    !--
    elemental type(polynomial) function psPlus(lhs,rhs) result(p)
    class(polynomial),intent(in)::      lhs
    real(rp),intent(in)::               rhs
        allocate(p%coef_ , source = lhs%coef_)
        p%coef_(0) = p%coef_(0) + rhs
    end function psPlus
    !--
    elemental type(polynomial) function spPlus(lhs,rhs) result(p)
    real(rp),intent(in)::               lhs
    class(polynomial),intent(in)::      rhs
        p = rhs + lhs
    end function spPlus
    !--
    elemental type(polynomial) function ppPlus(lhs,rhs) result(p)
    class(polynomial),intent(in)::      lhs,rhs
    integer(ip)::                       lu,ru
        lu = ubound(lhs%coef_,dim=1)
        ru = ubound(rhs%coef_,dim=1)
        if(lu > ru) then
            p = lhs
            p%coef_(0:ru) = p%coef_(0:ru) + rhs%coef_(0:ru)
        else
            p = rhs
            p%coef_(0:lu) = p%coef_(0:lu) + lhs%coef_(0:lu)
            if(lu==ru) p = p%contract()
        endif
    end function ppPlus
    !--
    elemental type(polynomial) function negativePoly(rhs) result(p)
    class(polynomial),intent(in)::      rhs
        allocate(p%coef_,source=rhs%coef_)
        p%coef_(:) = - p%coef_(:)
    end function negativePoly
    !--
    elemental type(polynomial) function ppMinus(lhs,rhs) result(p)
    class(polynomial),intent(in)::      lhs,rhs
        p = lhs + (-rhs)
    end function ppMinus
    !--
    elemental type(polynomial) function ppMultiply(lhs,rhs) result(p)
    class(polynomial),intent(in)::      lhs
    type(polynomial),intent(in)::       rhs
    integer(ip)::                       i,j,ld,rd
        ld = lhs%degree();  rd =rhs%degree()
        call p%init(ld+rd)
        do j=0,ld
            do i=0,rd
                p%coef_(i+j) = p%coef_(i+j) + lhs%coef_(j)*rhs%coef_(i)
            enddo
        enddo
    end function ppMultiply
    !--
    elemental type(polynomial) function psMultiply(lhs,rhs) result(p)
    class(polynomial),intent(in)::      lhs
    real(rp),intent(in)::               rhs
        if(rhs/=zero) then
            allocate(p%coef_,source=lhs%coef_)
            p%coef_ = rhs * p%coef_
        else
            p = zeroPolynomial()
        endif
    end function psMultiply
    !--
    elemental type(polynomial) function spMultiply(lhs,rhs) result(p)
    real(rp),intent(in)::               lhs
    class(polynomial),intent(in)::      rhs
        p = rhs * lhs
    end function spMultiply
    !--
    elemental type(polynomial) function psDivide(lhs,rhs) result(p)
    class(polynomial),intent(in)::      lhs
    real(rp),intent(in)::               rhs
        allocate(p%coef_,source=lhs%coef_)
        p%coef_ = p%coef_ / rhs
    end function psDivide
    !--
    elemental logical(lp) function ppjdeq(lhs,rhs) result(l)
    class(polynomial),intent(in)::  lhs,rhs
    integer(ip)::                   n,i
        l = .true.
        n = lhs%degree()
        if(n==rhs%degree()) then
            do i=0,n
                if(lhs%coef(i)==rhs%coef(i)) cycle
                l = .false.; exit
            enddo
        else
            l = .false.
        endif
    end function ppjdeq
    
    
    !-----------------------------------------------
    pure real(rp) function integratePolynomial(p,lo,up) result(r)
    class(polynomial),intent(in)::  p
    real(rp),intent(in)::           lo,up
        r = p%integral(lo,up)
    end function integratePolynomial
    
    
    
!------------------------------------------------------------------
!specified polynomials
    
    !---------------------
    elemental type(polynomial) function zeroPolynomial() result(z)
        allocate(z%coef_(0:0)); z%coef_ = zero
    end function zeroPolynomial
    
    
    
    !---------------------
    !n P_n = (2n-1) x P_{n-1} - (n-1) P_{n-2}
    elemental type(polynomial) function LegendrePolynomial(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp,1._rp]
        else
            x   = [0._rp,1._rp]
            tm2 = [1._rp]
            tm1 = x
            do i = 2 , n
                poly = (2._rp*i - 1._rp) * x * tm1 - (i - 1._rp) * tm2
                poly = poly / real(i,kind=rp)
                tm2 = tm1
                tm1 = poly
            enddo
        endif
    end function LegendrePolynomial
    
    !--
    pure function LegendrePolynomialSet(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp,1._rp]
        else
            x   = [0._rp,1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp,1._rp]
            do i = 2 , n
                poly(i) = (2._rp*i - 1._rp) * x * poly(i-1) - (i - 1._rp) * poly(i-2)
                poly(i) = poly(i) / real(i,kind=rp)
            enddo
        endif
    end function LegendrePolynomialSet
    
    !--
    elemental type(polynomial) function normalLegendrePolynomial(n) result(poly)
    integer(ip),intent(in)::                n
        poly = sqrt(dfloat(2*n+1)/2.d0) * LegendrePolynomial(n)
    end function normalLegendrePolynomial
    
    !--
    pure function normalLegendrePolynomialSet(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial),dimension(0:n)::       poly
    integer(ip)::                           i
        poly = LegendrePolynomialSet(n)
        do i=0,n
            poly(i) = sqrt(dfloat(2*i+1)/2.d0) * poly(i)
        enddo
    end function normalLegendrePolynomialSet
    

!wiki(chebyshev polynomial) T(n) for first kind and U(n) for second kind
    elemental type(polynomial) function ChebyshevPolynomialT(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp,1._rp]
        else
            x   = [0._rp,1._rp]
            tm2 = [1._rp]
            tm1 = x
            do i = 2 , n
                poly = 2._rp*x*tm1 - tm2
                tm2 = tm1
                tm1 = poly
            enddo
        endif
    end function ChebyshevPolynomialT
    
    !--
    pure function ChebyshevPolynomialTset(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp,1._rp]
        else
            x       = [0._rp,1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp,1._rp]
            do i = 2 , n
                poly(i) = 2._rp*x*poly(i-1) - poly(i-2)
            enddo
        endif
    end function ChebyshevPolynomialTset
    
    
!--compute sum_0^n(c*ChebPoly), sum = c(0)*T(0)%funcval(x)+......+c(n)*T(n)%funcval(x)
!use Clenshaw algorithm, see wiki
!https://github.com/chebfun/chebfun/blob/development/%40chebtech/clenshaw.m
    pure real(rp) function ChebyshevPolynomialT_Clenshaw(x,c) result(s)
    real(rp),intent(in)::               x
    real(rp),dimension(0:),intent(in):: c
    real(rp)::                          bk1,bk2,b,x2
    integer(ip)::                       n,k
        bk1 = 0._rp
        bk2 = bk1
        x2 = 2._rp * x   !double
        n = ubound(c,dim=1)
        do k=n,2,-2
            bk2 = c(k) + x2*bk1 - bk2
            bk1 = c(k-1) + x2*bk2 - bk1
        enddo
        if(ibits(n,0,1)==1) then
            b = bk1
            bk1 = c(1) + x2*bk1 - bk2
            bk2 = b
        endif
        s = c(0) + x * bk1 - bk2
    end function ChebyshevPolynomialT_Clenshaw
    
    
    !---------------------
    !H_n = 2x H_{n-1} - 2(n-1) H_{n-2}
    elemental type(polynomial) function HermitePhysicistsPolynomial(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp,2._rp]
        else
            x   = [0._rp,1._rp]
            tm2 = [1._rp]
            tm1 = 2._rp*x
            do i = 2 , n
                poly = 2._rp* x * tm1 - 2._rp*(i - 1._rp) * tm2
                tm2 = tm1
                tm1 = poly
            enddo
        endif
    end function HermitePhysicistsPolynomial
    
    !--
    pure function HermitePhysicistsPolynomialSet(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp,2._rp]
        else
            x   = [0._rp,1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp,2._rp]
            do i = 2 , n
                poly(i) = 2._rp* x * poly(i-1) - 2._rp*(i - 1._rp) * poly(i-2)
            enddo
        endif
    end function HermitePhysicistsPolynomialSet
    
    !--
    elemental type(polynomial) function normalHermitePhysicistsPolynomial(n) result(poly)
    integer(ip),intent(in)::                n
        poly = ( 1._rp / sqrt( spi * 2**n * factorial(n)) ) * HermitePhysicistsPolynomial(n)
    end function normalHermitePhysicistsPolynomial
    
    !--
    pure function normalHermitePhysicistsPolynomialSet(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial),dimension(0:n)::       poly
    integer(ip)::                           i
        poly = HermitePhysicistsPolynomialSet(n)
        do i=0,n
            poly(i) = ( 1._rp / sqrt( spi * 2**n * factorial(n)) ) * poly(i)
        enddo
    end function normalHermitePhysicistsPolynomialSet
    
     !---------------------
    !H_n = x H_{n-1} - (n-1) H_{n-2}
    elemental type(polynomial) function HermiteProbabilistsPolynomial(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial)::                  tm2,tm1,x
    integer(ip)::                       i
        if(n<=0) then
            poly = [1._rp]
        elseif(n==1) then
            poly = [0._rp,1._rp]
        else
            x   = [0._rp,1._rp]
            tm2 = [1._rp]
            tm1 = x
            do i = 2 , n
                poly = x * tm1 - (i - 1._rp) * tm2
                tm2 = tm1
                tm1 = poly
            enddo
        endif
    end function HermiteProbabilistsPolynomial
    
    !--
    pure function HermiteProbabilistsPolynomialSet(n) result(poly)
    integer(ip),intent(in)::            n
    type(polynomial),dimension(0:n)::   poly
    type(polynomial)::                  x
    integer(ip)::                       i
        if(n<=0) then
            poly(0) = [1._rp]
        elseif(n==1) then
            poly(0) = [1._rp]
            poly(1) = [0._rp,1._rp]
        else
            x   = [0._rp,1._rp]
            poly(0) = [1._rp]
            poly(1) = [0._rp,1._rp]
            do i = 2 , n
                poly(i) = x * poly(i-1) - (i - 1._rp) * poly(i-2)
            enddo
        endif
    end function HermiteProbabilistsPolynomialSet
    
    !--
    elemental type(polynomial) function normalHermiteProbabilistsPolynomial(n) result(poly)
    integer(ip),intent(in)::                n
        poly = ( 1._rp / sqrt( sqrt(2._rp) * spi * factorial(n)) ) * HermiteProbabilistsPolynomial(n)
    end function normalHermiteProbabilistsPolynomial
    
    !--
    pure function normalHermiteProbabilistsPolynomialSet(n) result(poly)
    integer(ip),intent(in)::                n
    type(polynomial),dimension(0:n)::       poly
    integer(ip)::                           i
        poly = HermiteProbabilistsPolynomialSet(n)
        do i=0,n
            poly(i) = ( 1._rp / sqrt( sqrt(2._rp) * spi * factorial(n)) ) * poly(i)
        enddo
    end function normalHermiteProbabilistsPolynomialSet
    
    
end module polynomial_