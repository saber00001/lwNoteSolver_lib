module arrayOpsLib
use constants, only:ip,rp,zero
implicit none

    private
    public:: operator(.ip.),operator(.op.),operator(.cpv.),operator(.cps.)
    public:: operator(-),operator(+),operator(*)
    public:: magSqr,mag
    public:: diag,trace
    
!---------------------------------------------
    interface operator(*)
        procedure:: mvinnerproduct
        procedure:: vminnerproduct
    end interface 
    
    interface operator(.ip.)
        procedure:: vvinnerproduct
        procedure:: mminnerproduct
        procedure:: mvinnerproduct
        procedure:: vminnerproduct
    end interface
    
    interface operator(.op.)
        procedure:: vvouterproduct
    end interface
    
    interface operator(.cpv.)
        procedure:: vvcrossproduct3
    end interface
    
    interface operator(.cps.)
        procedure:: vvcrossproduct2
    end interface
    
    interface operator(+)
        procedure:: mvplus
        procedure:: vmplus
    end interface

    interface operator(-)
        procedure:: mvminus
        procedure:: vmminus
    end interface
    
    
    
    
contains


    !---inner product
    !--
    pure real(rp) function vvinnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:),intent(in)::  lhs,rhs
        p = dot_product(lhs,rhs)
    end function vvinnerproduct
    
    !--
    pure function mminnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:,:),intent(in)::lhs,rhs
    real(rp),dimension(size(lhs,dim=1),size(rhs,dim=2))::p
        p = matmul(lhs,rhs)
    end function mminnerproduct
    
    !--
    pure function mvinnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:,:),intent(in)::lhs
    real(rp),dimension(:),intent(in)::  rhs
    real(rp),dimension(size(lhs,dim=1))::p
        p = matmul(lhs,rhs)
    end function mvinnerproduct
    
    !--
    pure function vminnerproduct(lhs,rhs) result(p)
    real(rp),dimension(:),intent(in)::  lhs
    real(rp),dimension(:,:),intent(in)::rhs
    real(rp),dimension(size(rhs,dim=2))::p
        p = matmul(lhs,rhs)
    end function vminnerproduct
    
    
    
    !---outer product
    pure function vvouterproduct(lhs,rhs) result(p)
    real(rp),dimension(:),intent(in)::  lhs,rhs
    real(rp),dimension(size(lhs),size(rhs)):: p
    integer(ip)::                       i,j
        do j=1,size(rhs)
            do i=1,size(lhs)
                p(i,j) = lhs(i) * rhs(j)
            enddo
        enddo
    end function vvouterproduct
    
    
    
    !---cross product
    pure function vvcrossproduct3(lhs,rhs) result(p)
    real(rp),dimension(3),intent(in)::  lhs,rhs
    real(rp),dimension(3)::             p
        p(1)  =   lhs(2) * rhs(3)   -   lhs(3) * rhs(2)
        p(2)  =   lhs(3) * rhs(1)   -   lhs(1) * rhs(3)
        p(3)  =   lhs(1) * rhs(2)   -   lhs(2) * rhs(1)
    end function vvcrossproduct3
    
    pure real(rp) function vvcrossproduct2(lhs,rhs) result(p)
    real(rp),dimension(2),intent(in)::  lhs,rhs
        p = lhs(1) * rhs(2) - lhs(2) * rhs(1)
    end function vvcrossproduct2
    
    
    !some special plus
    !--
    pure function mvplus(lhs,rhs) result(m)
    real(rp),dimension(:,:),intent(in)::lhs
    real(rp),dimension(:),intent(in)::  rhs
    real(rp),dimension(size(lhs,dim=1),size(lhs,dim=2))::m
    integer(ip)::                       i,j,d
        m = lhs
        d = size(rhs)
        forall(i=1:d,j=1:d,i==j) m(i,j) = m(i,j) + rhs(i)
    end function mvplus
    
    !--
    pure function vmplus(lhs,rhs) result(m)
    real(rp),dimension(:),intent(in)::  lhs
    real(rp),dimension(:,:),intent(in)::rhs
    real(rp),dimension(size(rhs,dim=1),size(rhs,dim=2))::m
    integer(ip)::                       i,j,d
        m = rhs + lhs
    end function vmplus
    
    !--
    pure function mvminus(lhs,rhs) result(m)
    real(rp),dimension(:,:),intent(in)::lhs
    real(rp),dimension(:),intent(in)::  rhs
    real(rp),dimension(size(lhs,dim=1),size(lhs,dim=2))::m
        m = lhs + ( - rhs )
    end function mvminus
    
    !--
    pure function vmminus(lhs,rhs) result(m)
    real(rp),dimension(:),intent(in)::  lhs
    real(rp),dimension(:,:),intent(in)::rhs
    real(rp),dimension(size(rhs,dim=1),size(rhs,dim=2))::m
        m = ( - rhs ) + lhs
    end function vmminus
    
    
    !--
    pure real(rp) function magSqr(v) result(ms)
    real(rp),dimension(:),intent(in)::  v
        ms = sum( v**2 )
    end function magSqr
    
    pure real(rp) function mag(v)
    real(rp),dimension(:),intent(in)::  v
        mag = sqrt( magSqr(v) )
    end function mag
    
    pure function diag(m)
    real(rp),dimension(:,:),intent(in)::m
    real(rp),dimension(min(size(m,dim=1),size(m,dim=2))):: diag
    integer(ip)::                       i
        do i = 1,size(diag)
            diag(i) = m(i,i)
        enddo
    end function diag
    
    pure real(rp) function trace(m)
    real(rp),dimension(:,:),intent(in)::m
    integer(ip)::                       i
        trace = zero
        do i=1,min(size(m,dim=1),size(m,dim=2))
            trace = trace + m(i,i)
        enddo
    end function trace
    
end module arrayOpsLib