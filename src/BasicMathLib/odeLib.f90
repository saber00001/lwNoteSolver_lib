!i want to write some basic operations for ordinary differetial equations
!as the basic lib for Type<timedriver> or others
    

!Here the odelib contains Euler solver and Runge-Kutta solver that uses the fourth order TVD format.
module odelib
use constants
implicit none

    
    private
    public::    odeEuler
    public::    odeRK4,odeSystemRK4
    public::    odeRK2
    
    
!-------------------------------------------------------------------    
    interface odeEuler
        procedure:: odeEuler_1step
        procedure:: odeEuler_nstep
    end interface odeEuler
    
!-------------------------------------------------------------------
    interface odeRK4
        procedure:: odeRK4_TVD_1step
        procedure:: odeRK4_TVD_nstep
    end interface odeRK4

!-------------------------------------------------------------------
    interface odeSystemRK4
        procedure:: odeSystemRK4_TVD_1step
    end interface odeSystemRK4
    
    
    
!-------------------------------------------------------------------
    abstract interface
        pure real(rp) function absdydx(x,y) result(dydx)
        import:: rp
        real(rp),intent(in)::   x,y 
        end function absdydx
        
        pure real(rp) function absdydxSystem(i,x,y) result(dydx)
        import:: rp,ip
        integer(ip),intent(in)::    i
        real(rp),intent(in)::       x,y
        end function absdydxSystem
    end interface

    
contains

    
    pure real(rp) function odeEuler_1step(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
        y = y0 + dydx(x0,y0)*dx
    end function odeEuler_1step
    
    pure function odeEuler_nstep(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeEuler(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeEuler(dydx , dx , x0+dfloat(i)*dx , y(i-1))
        end do
    end function odeEuler_nstep
    
!------------------------------------------------------------------- 
    pure real(rp) function odeRK4_TVD_1step(dydx,dx,x0,y0) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    real(rp)::                  k1,k2,k3,k4
        k1 = dydx(x0,y0)
        k2 = dydx(x0 + 0.5d0*dx,y0 + 0.5d0*dx*k1)
        k3 = dydx(x0 + 0.5d0*dx,y0 + 0.5d0*dx*k2)
        k4 = dydx(x0 + dx,y0 + dx*k3)
        y = y0 + (1.d0/6.d0)*dx*(k1 + 2.d0*k2 + 2.d0*k3 + k4)
    end function odeRK4_TVD_1step

    pure function odeRK4_TVD_nstep(dydx,dx,x0,y0,n) result(y)
    procedure(absdydx)::        dydx
    real(rp),intent(in)::       dx,x0,y0
    integer(ip),intent(in)::    n
    real(rp),dimension(n)::     y
    integer(ip)::               i
        y(1)=odeRK4(dydx,dx,x0,y0)
        do i=2,n
            y(i) = odeRK4(dydx , dx , x0+dfloat(i)*dx , y(i-1))
        end do
    end function odeRK4_TVD_nstep
    
    !------------ 
    pure subroutine odeSystemRK4_TVD_1step(dydx,dx,x0,y0,y)
    procedure(absdydxSystem)::          dydx
    real(rp),intent(in)::               dx
    real(rp),dimension(:),intent(in)::  x0,y0
    real(rp),dimension(:),intent(out):: y
    real(rp)::                          k1,k2,k3,k4
    integer(ip)::                       i
        do i=1,size(x0)
            y(i) = odeRK4(d,dx,x0(i),y0(i))
        enddo
    contains
        pure real(rp) function d(x,y)
        real(rp),intent(in)::   x,y
            d = dydx(i,x,y)
        end function d
    end subroutine odeSystemRK4_TVD_1step
    
    
!------------------------------------------------------------------- 
    pure real(rp) function odeRK2(t, tn, x, n) result(y)
    real(rp),intent(in)::           t,tn,x   
    integer(ip),intent(in)::        n   
    integer(ip)::                   m,j
    real(rp)::                      h,k1,k2,t1,x1
        t1 = t
        x1 = x
        m = 10
        h = (tn-t)/m
        do j = 1,m
            k1 = -h/(sqrt(2.*n+1-x1**2) - .5*x1*sin(2.*t1)/(2.*n+1-x1**2))
            t1 = t1 + h
            k2 = -h/(sqrt(2.*n+1-(x1+k1)**2) - .5*x1*sin(2.*t1)/(2.*n+1-(x1+k1)**2))
            x1 = x1 + .5*(k1 + k2)
        end do
        y = x1  
    end function odeRK2
    
end module odelib
