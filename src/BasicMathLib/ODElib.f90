!i want to write some basic operations for ordinary differetial equations
!as the basic lib for Type<timedriver> or others
    

!nothing temporarily
module ODElib
use constants
implicit none

    
    private
    public::    ODEeuler
    public::    absdydx
    
    
    
!-------------------------------------------------------------------    
    interface ODEeuler
        procedure:: ODEeuler_1step
        procedure:: ODEeuler_nstep
    end interface
    
!-------------------------------------------------------------------
    abstract interface
        pure real(rp) function absdydx(x,y) result(dydx)
        import:: rp
        real(rp),intent(in)::   x,y 
        end function absdydx
    end interface

    
contains

    
    pure real(rp) function ODEeuler_1step(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
        y = y0 + dydx(x0,y0)*dx
    end function ODEeuler_1step
    
    pure function ODEeuler_nstep(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1) = y0
        do i=2,n
            y(i) = ODEeuler(dydx,dx,x0+(i-1)*dx,y(i-1))
        enddo
    end function ODEeuler_nstep
    
end module ODElib