module arrayOpsLib
use constants, only:ip,rp,lp,zero,swap
implicit none

    private
    public:: operator(.ip.),operator(.op.),operator(.cpv.),operator(.cps.)
    public:: operator(-),operator(+),operator(*),operator(.eqvl.)
    public:: magSqr,mag,angle,normal,para,orth,norm,polyval
    public:: trace,diag,sort,cumprod
    
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
    
    interface diag
        procedure:: diagCreateMatrix
        procedure:: diagExtractElement
    end interface
    
    interface sort
        procedure:: sortOneDimension
        procedure:: sortTwoDimension
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
    
    pure function normal(v)
    real(rp),dimension(:),intent(in)::  v
    real(rp),dimension(size(v))::       normal
        normal = v / mag(v)
    end function normal
    
    pure function para(v1,v2)
    real(rp),dimension(:),intent(in)::  v1,v2
    real(rp),dimension(size(v1))::      para
        para = (v1.ip.v2)*normal(v2)
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
        if(p==2) then
            norm = norm2(v)
        elseif(p==1) then
            norm = sum(abs(v))
        elseif(p==0) then
            norm = count(v /= [(0._rp,i=1,size(v))])
        elseif(p>2) then
            norm = sum(abs(v)**p)**(1._rp/p)
        else
            norm = maxval(abs(v))
        endif
    end function norm
    
    !---------
    pure real(rp) function polyval(v,x)
    real(rp),dimension(0:),intent(in):: v
    real(rp),intent(in)::               x
    integer(ip)::                       i
        polyval = 0._rp
        do i=0,ubound(v,dim=1)
            polyval = polyval + v(i) * x**i 
        enddo
    end function polyval
    
    !---------
    pure function diagCreateMatrix(m,k)
    real(rp),dimension(:),intent(in)::                    m
    integer(ip),intent(in)::                              k
    real(rp),dimension(:,:),allocatable::                 diagCreateMatrix
    integer(ip)::                                         i,j,n,p
        p = abs(k); n = size(m)
        allocate(diagCreateMatrix(p+n,p+n))
        diagCreateMatrix = 0._rp
        if(k>=0) then
            do i=1,n
                diagCreateMatrix(i,i+k) = m(i)
            enddo
        elseif(k<0) then
            do j=1,n
                diagCreateMatrix(j+p,j) = m(j)
            enddo
        endif 
    end function diagCreateMatrix
    !--
    pure function diagExtractElement(m)
    real(rp),dimension(:,:),intent(in)::    m
    real(rp),dimension(min(size(m,dim=1),size(m,dim=2))):: diagExtractElement
    integer(ip)::                           i
        forall(i=1:size(diagExtractElement)) diagExtractElement(i)=m(i,i)
    end function diagExtractElement
    
    !----
    pure real(rp) function trace(m)
    real(rp),dimension(:,:),intent(in)::    m
    integer(ip)::                           i
        trace = zero
        do i=1,min(size(m,dim=1),size(m,dim=2))
            trace = trace + m(i,i)
        enddo
    end function trace
    
    !---
    pure subroutine sortOneDimension(array,location)
    real(rp),dimension(:),intent(inout)::       array
    integer(ip),dimension(:),intent(out)::      location
    integer(ip)::                               i,k,n
        n = size(array) 
        do i=1,n
            location(i) = i
        end do
        do i=n-1,1,-1
            do k=1,i
                if(array(k)>array(k+1)) then
                    call swap(array(k),array(k+1))
                    call swap(location(k),location(k+1))
                end if
            end do
        end do
    end subroutine sortOneDimension
    !--
    pure subroutine sortTwoDimension(array,location)
    real(rp),dimension(:,:),intent(inout)::     array
    integer(ip),dimension(:,:),intent(out)::    location
    integer(ip)::                               j
        do j=1,size(array,dim=2)
            call sort(array(:,j),location(:,j))
        end do
    end subroutine sortTwoDimension
    
    !-----
    pure function cumprod(x) 
    real(rp),dimension(:),intent(in)::  x
    real(rp),dimension(size(x))::       cumprod
    integer(ip)::                       i
        cumprod(1) = x(1)
        do i=2,size(x)
            cumprod(i) = cumprod(i-1)*x(i)  
        end do   
    end function cumprod
    
    
    
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