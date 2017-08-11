module arrayOpsLib
use constants, only:ip,rp,lp,zero
implicit none

    private
    public:: operator(.ip.),operator(.op.),operator(.cpv.),operator(.cps.)
    public:: operator(-),operator(+),operator(*),operator(.eqvl.)
    public:: magSqr,mag,angle,unit,para,orth,norm
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
    
    interface operator(.eqvl.)
        procedure:: realiseq
        procedure:: integeriseq
    end interface
    
contains

!-------------------------------------------------------
    pure real(rp) function mag(v)
    real(rp),dimension(:),intent(in)::  v
        mag = norm2(v)
    end function mag
    
    pure real(rp) function magSqr(v) result(ms)
    real(rp),dimension(:),intent(in)::  v
        ms = sum(v**2)
    end function magSqr
    
    pure real(rp) function angle(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
        angle = acos((v1.ip.v2)/mag(v1)/mag(v2))
    end function angle
    
    pure function unit(v)
    real(rp),dimension(:),intent(in)::  v
    real(rp),dimension(size(v))::       unit
        unit = v / mag(v)
    end function unit
    
    pure function para(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
    real(rp),dimension(size(v1))::      para
        para = (v1.ip.v2)*unit(v2)
    end function para
    
    pure function orth(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
    real(rp),dimension(size(v1))::      orth
        orth = v1 - para(v1,v2)
    end function orth
    
    !Lp norm
    pure real(rp) function norm(v,p)
    real(rp),dimension(:),intent(in)::  v
    integer(ip),intent(in)::            p
    integer(ip)::                       i
        if(p>0) then
            norm = sum(abs(v)**p)**(1._rp/p)
        elseif(p==0) then
            norm = count( v /= [(0._rp,i=1,size(v))] )
        else!let minus p be the infinity norm
            norm = maxval(abs(v)) 
        endif
    end function norm
    
!---------------------------------------------    
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







!---------------------inner product--------------------------------
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
    pure logical(lp) function realiseq(lhs,rhs) result(r)
    real(rp),dimension(:),intent(in)::  lhs,rhs
        r = .false.
        if(size(lhs)==size(rhs)) r = all(lhs==rhs)
    end function realiseq
    
    !--
    pure logical(lp) function integeriseq(lhs,rhs) result(r)
    integer(ip),dimension(:),intent(in)::lhs,rhs
        r = .false.
        if(size(lhs)==size(rhs)) r = all(lhs==rhs)
    end function integeriseq
    
end module arrayOpsLib