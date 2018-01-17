!reconstruction method is mainly derived from conservation equation
!but we still can consider reconstruction as a special interpolation: 
!<see the Shu's manualscripts about weno reconstruction and weno interpolation>
!here we put some kernel of reconstruction method for general use
module interpolationLib
use constants
use arrayOpsLib
use vector_
implicit none

    private
    !lgrInt = Lagrange
    public:: lgrInt,centerlgrInt,extralgrInt
    public:: muscl2c,muscl2c_SpOrthoNormalCoefsArray
    public:: weno5z
    
    
    !---------------------------------------------------
    interface lgrInt
        procedure:: lgrInt1d
    end interface lgrInt
    
    interface extralgrInt
        !(x--1--s1--1--s2)
        procedure:: extralgrInt2_scalar
        procedure:: extralgrInt2_vector
        procedure:: extralgrInt2_array
        !(x--1--s1--1--s2--1--s3)
        procedure:: extralgrInt3_scalar
        procedure:: extralgrInt3_vector
        procedure:: extralgrInt3_array
    end interface extralgrInt
    
    interface centerlgrInt
        !(l2--1--l1--1--x--1--r1--1--r2)
        procedure:: centerlgrInt2_scalar
        procedure:: centerlgrInt2_vector
        procedure:: centerlgrInt2_array
    end interface centerlgrInt
    
    interface muscl2c
        !(p(1)--1--p(2)--1/2--x--1/2--p(3))
        procedure muscl2c_scalar
        procedure muscl2c_array
        procedure muscl2c_specArrayDiscret
    end interface muscl2c
    
    interface weno5z
        !(p(1)--1--p(2)--1--p(3)--1/2--x--1/2--p(4)--1--p(5))
        procedure:: weno5z_scalar
        procedure:: weno5z_array
    end interface weno5z
    

contains
    
    pure real(rp) function lgrInt1d(x,vals,valsx) result(v)
    real(rp),intent(in)::               x
    real(rp),dimension(:),intent(in)::  vals,valsx
    integer(ip)::                       i,j,n
    real(rp)::                          cu,cl
        n = min(size(vals),size(valsx))
        v = zero
        do j=1,n
            cu = 1._rp
            cl = 1._rp
            do i=1,j-1
                cu = cu * (x-valsx(i))
                cl = cl * (valsx(j) - valsx(i))
            enddo
            do i=j+1,n
                cu = cu * (x-valsx(i))
                cl = cl * (valsx(j) - valsx(i))
            enddo
            v = v + vals(j) * cu / cl
        enddo
    end function lgrInt1d

    
    !-----------------------------------------------------------
    pure real(rp) function extralgrInt2_scalar(s1,s2) result(x)
    real(rp),intent(in):: s1,s2
        x = 2._rp * s1 - s2
    end function extralgrInt2_scalar
    
    pure type(vector) function extralgrInt2_vector(s1,s2) result(x)
    type(vector),intent(in):: s1,s2
        x = 2._rp * s1 - s2
    end function extralgrInt2_vector
    
    pure function extralgrInt2_array(s1,s2) result(x)
    real(rp),dimension(:),intent(in)::  s1,s2
    real(rp),dimension(size(s1))::      x
        x = 2._rp * s1 - s2
    end function extralgrInt2_array
    
    pure real(rp) function extralgrInt3_scalar(s1,s2,s3) result(x)
    real(rp),intent(in):: s1,s2,s3
        x = (11._rp * s1 - 7._rp * s2 + 2._rp * s3) / 6._rp
    end function extralgrInt3_scalar
    
    pure type(vector) function extralgrInt3_vector(s1,s2,s3) result(x)
    type(vector),intent(in):: s1,s2,s3
        x = (11._rp * s1 - 7._rp * s2 + 2._rp * s3) / 6._rp
    end function extralgrInt3_vector
    
    pure function extralgrInt3_array(s1,s2,s3) result(x)
    real(rp),dimension(:),intent(in)::  s1,s2,s3
    real(rp),dimension(size(s1))::      x
        x = (11._rp * s1 - 7._rp * s2 + 2._rp * s3) / 6._rp
    end function extralgrInt3_array
    
    
    !-------
    pure real(rp) function centerlgrInt2_scalar(l2,l1,r1,r2) result(x)
    real(rp),intent(in):: l2,l1,r1,r2
        x = ( 4._rp * ( l1 + r1) - l2 - r2 ) / 6._rp
    end function centerlgrInt2_scalar
    
    pure type(vector) function centerlgrInt2_vector(l2,l1,r1,r2) result(x)
    type(vector),intent(in):: l2,l1,r1,r2
        x = ( 4._rp * ( l1 + r1) - l2 - r2 ) / 6._rp
    end function centerlgrInt2_vector
    
    pure function centerlgrInt2_array(l2,l1,r1,r2) result(x)
    real(rp),dimension(:),intent(in)::  l2,l1,r1,r2
    real(rp),dimension(size(l2))::      x
        x = ( 4._rp * ( l1 + r1) - l2 - r2 ) / 6._rp
    end function centerlgrInt2_array
    
     
    !---------------------------------------------------------------
    !refer to <Semi-discrete central-upwind scheme for hyperbolic conservation laws and Hamilton-Jacobi equations>
    pure real(rp) function muscl2c_scalar(f,theta) result(x)
    real(rp),dimension(:),intent(in)::  f
    real(rp),intent(in)::               theta
    real(rp),dimension(3)::             delta
    real(rp)::                          r,beta
        delta(1) = f(2) - f(1)
        delta(2) = f(3) - f(2)
        delta(3) = ( f(3) - f(1) ) / 2._rp
        r = delta(2) / (delta(1) + GlobalEps)
        beta = delta(3) / (delta(1) + GlobalEps)
        x = f(2) + 0.5_rp * max(zero,min(theta,theta*r,beta)) * delta(1)
    end function  muscl2c_scalar
    !--
    pure function muscl2c_array(f,theta) result(x)
    real(rp),dimension(:,:),intent(in)::    f
    real(rp),intent(in)::                   theta
    real(rp),dimension(size(f,dim=1))::     x,r,beta
    real(rp),dimension(size(f,dim=1),3)::   delta
        delta(:,1) = f(:,2) - f(:,1)
        delta(:,2) = f(:,3) - f(:,2)
        delta(:,3) = ( f(:,3) - f(:,1) ) / 2._rp
        r = delta(:,2) / (delta(:,1) + GlobalEps)
        beta = delta(:,3) / (delta(:,1) + GlobalEps)
        x = f(:,2) + 0.5_rp * max(zero,min(theta,theta*r,beta)) * delta(:,1)
    end function  muscl2c_array
    !--
    pure function muscl2c_specArrayDiscret(f,theta) result(x)
    real(rp),dimension(:,:,:),intent(in)::              f
    real(rp),intent(in)::                               theta
    real(rp),dimension(size(f,dim=1),size(f,dim=2))::   x,r,beta
    real(rp),dimension(size(f,dim=1),size(f,dim=2),3):: delta
        delta(:,:,1) = f(:,:,2) - f(:,:,1)
        delta(:,:,2) = f(:,:,3) - f(:,:,2)
        delta(:,:,3) = ( f(:,:,3) - f(:,:,1) ) / 2._rp
        !discreted approach
        r = delta(:,:,2) / (delta(:,:,1) + GlobalEps)
        beta = delta(:,:,3) / (delta(:,:,1) + GlobalEps)
        x = f(:,:,2) + 0.5_rp * max(zero,min(theta,theta*r,beta)) * delta(:,:,1)
    end function  muscl2c_specArrayDiscret
    !--
    pure function muscl2c_SpOrthoNormalCoefsArray(f,theta) result(x)
    real(rp),dimension(:,:,:),intent(in)::              f
    real(rp),intent(in)::                               theta
    real(rp),dimension(size(f,dim=1),size(f,dim=2))::   x,r,beta
    real(rp),dimension(size(f,dim=1),size(f,dim=2),3):: delta
    integer(ip)::                                       i,loc
        delta(:,:,1) = f(:,:,2) - f(:,:,1)
        delta(:,:,2) = f(:,:,3) - f(:,:,2)
        delta(:,:,3) = ( f(:,:,3) - f(:,:,1) ) / 2._rp
        !functional approah 1. check direction; 2. check variation
        do i=1,size(f,dim=2)
            if( ( delta(:,i,1) .ip. delta(:,i,2) ) < 0._rp ) then
                x(:,i) = f(:,i,2)
            else
                loc = minloc([theta*norm2(delta(:,i,1)),theta*norm2(delta(:,i,2)),norm2(delta(:,i,3))],dim=1)
                if(loc==3) then
                    x(:,i) = f(:,i,2) + 0.5_rp * delta(:,i,loc)
                else
                    x(:,i) = f(:,i,2) + 0.5_rp * delta(:,i,loc) * theta
                endif
            endif
        enddo
    end function  muscl2c_SpOrthoNormalCoefsArray

    
    !---------
    !refer to <jcp - An improved WENO-Z scheme>
    !refer to <High OrderWeighted Essentially Nonoscillatory Schemes for Convection Dominated Problems>
    !this is reconstruction rather than interpolation which accounts for the spatial average in conservative PDF
    !reconstruction is a special interpolation based on averaged information rather than point information
    pure function weno5z_scalar(f) result(x)
    real(rp),dimension(:),intent(in)::f
    real(rp)::      x,s0,s1,s2,w0,w1,w2,w
    real(rp),parameter::    eps = GlobalEps * 10._rp
    integer(ip),parameter:: p   = 2
    real(rp),parameter::    alpha1 = 0.25_rp, alpha2 = 13._rp/12._rp
    real(rp),parameter::    a0 = 1._rp/3._rp, b0 = -7._rp/6._rp, c0 = 11._rp/6._rp
    real(rp),parameter::    a1 =-1._rp/6._rp, b1 =  5._rp/6._rp, c1 =  1._rp/3._rp
    real(rp),parameter::    a2 = 1._rp/3._rp, b2 =  5._rp/6._rp, c2 = -1._rp/6._rp
    
        s0 = alpha1 * (f(1) - 4._rp * f(2) + 3._rp * f(3)) ** 2     &
           + alpha2 * (f(1) - 2._rp * f(2) + f(3)  ) ** 2
                                
        s1 = alpha1 * (f(2) - f(4)) ** 2 &
           + alpha2 * (f(2) - 2._rp * f(3) + f(4)) ** 2     

        s2 = alpha1 * ( 3._rp * f(3) - 4._rp * f(4) + f(5)) ** 2    &
           + alpha2 * (f(3) - 2._rp * f(4) + f(5)) ** 2

        w0 = 0.1_rp * (1._rp + ((s0-s2)/(s0+eps))**p)
        w1 = 0.6_rp * (1._rp + ((s0-s2)/(s1+eps))**p)
        w2 = 0.3_rp * (1._rp + ((s0-s2)/(s2+eps))**p)    
        w = w0 + w1 + w2                    
        w0 = w0/w
        w1 = w1/w
        w2 = w2/w

        x  = w0 * (a0*f(1) + b0*f(2) + c0*f(3))    &
           + w1 * (a1*f(2) + b1*f(3) + c1*f(4))    &
           + w2 * (a2*f(3) + b2*f(4) + c2*f(5))
           
    end function weno5z_scalar
    
    !--dim(f,1) is field dimension and dim(f,2) is the node index
    pure function weno5z_array(f) result(x)
    real(rp),dimension(:,:),intent(in)::f
    real(rp),dimension(size(f,dim=1)):: x
    integer(ip)::                       i
        forall(i=1:size(f,dim=1)) x(i) = weno5z(f(i,:))
    end function weno5z_array
    
end module interpolationLib