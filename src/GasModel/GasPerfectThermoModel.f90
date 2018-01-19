!physical model branches into an extreme scale with various corrections.
!create an nested object to organize these models is a smarter way to refactor
!use the instatiated object always
module GasPerfectThermoModel_
use constants
implicit none

    private
    !instantiation
    public:: GasPerfectThermoModel
    !-----------------------------
    !assembly of all of models in gas
    !-----------------------------
    type GasPerfectThermoModel
        
        !private
        real(rp),allocatable,dimension(:)::     Coef_High_T_Range
        real(rp),allocatable,dimension(:)::     Coef_Low_T_Range
        
    contains
        generic::                   init    =>  init_n,     &
                                                init_n_val, &
                                                init_vals    
        procedure::                 init_n
        procedure::                 init_n_val
        procedure::                 init_vals
        !Sutherland Equation to calculate viscocity
        generic::                   Sutherland => Sutherland_custom,Sutherland_air
        procedure,nopass::          Sutherland_custom
        procedure,nopass::          Sutherland_air
        !-------------
        !idealGasModel
        !-------------
        !ideal Gas state Equation T(p,rho,R)
        procedure,nopass::          idealGasT

        !ideal gas isoEntropic flow
        procedure,nopass::          isoEntropicPressureRatio
        procedure,nopass::          isoEntropicTemperatureRatio
        
        procedure,nopass::          Cp
        procedure,nopass::          H
        procedure,nopass::          S
        
    end type GasPerfectThermoModel
        
contains
    
    !-------------
    pure subroutine init_n(this,ncoef)
    class(GasPerfectThermoModel),intent(out)::      this
    integer(ip),intent(in)::                        ncoef
    
        allocate(this%Coef_High_T_Range(ncoef))
        allocate(this%Coef_Low_T_Range(ncoef))
        this%Coef_High_T_Range          = 0._rp
        this%Coef_Low_T_Range           = 0._rp
    
    end subroutine init_n
    !--
    pure subroutine init_n_val(this,ncoef,Coef_H,Coef_L)
    class(GasPerfectThermoModel),intent(out)::     this
    integer(ip),intent(in)::                ncoef
    real(rp),intent(in)::                   Coef_H,Coef_L
    
        allocate(this%Coef_High_T_Range(ncoef))
        allocate(this%Coef_Low_T_Range(ncoef))
        this%Coef_High_T_Range          = Coef_H
        this%Coef_Low_T_Range           = Coef_L
    
    end subroutine init_n_val
    !--
    pure subroutine init_vals(this,Coef_H,Coef_L)
    class(GasPerfectThermoModel),intent(out):: this
    real(rp),dimension(:),intent(in)::  Coef_H,Coef_L
        
        allocate(this%Coef_High_T_Range,source=Coef_H)
        allocate(this%Coef_Low_T_Range,source =Coef_L)
        
    end subroutine init_vals
    
    !---------------
    !--
    elemental real(rp) function Sutherland_custom(T,T0,T1,Miu0) result(miu)
    real(rp),intent(in)::   T,T0,T1,Miu0
        miu = (T0 + T1) / (T + T1)
        miu = miu * Miu0 * sqrt( (T / T0)**3 )
    end function Sutherland_custom
    
    !--
    elemental real(rp) function Sutherland_air(T) result(miu)
    real(rp),intent(in)::   T
        miu = Sutherland_custom(T,273.15_rp,110.4_rp,1.7161_rp)
    end function Sutherland_air
    
!------------------------------------------------------------------------
    elemental real(rp) function idealGasT(rho,p,R) result(T)
    real(rp),intent(in)::   rho,p,R
        T = p / rho / R
    end function idealGasT
    
    !p0/p for ideal isoentropic flow
    elemental real(rp) function isoEntropicPressureRatio(ma,gm) result(pr)
    real(rp),intent(in)::   ma,gm
        pr = (1._rp + (gm-1._rp)/2._rp * ma**2 ) ** (gm/(gm-1._rp))
    end function isoEntropicPressureRatio
    
    !T0/T for ideal isoentropic flow
    elemental real(rp) function isoEntropicTemperatureRatio(ma,gm) result(TR)
    real(rp),intent(in)::   ma,gm
        TR = 1._rp + (gm-1._rp)/2._rp * ma**2
    end function isoEntropicTemperatureRatio    
    
!------------------------------------------------------------------------
    !specific heat capacity polynomial fitting
    pure real(rp) function Cp(this,t,r) result(shcp)
    class(GasPerfectThermoModel),intent(in)::  this
    real(rp),intent(in)::               t,r
    real(rp),parameter::                temperature_boundary=1000._rp
    
        if(t > temperature_boundary) then
            shcp = r*(this%Coef_High_T_range(1)*t**(-2) + this%Coef_High_T_range(2)*t**(-1) + this%Coef_High_T_range(3)&
                  + this%Coef_High_T_range(4)*t + this%Coef_High_T_range(5)*t**2&
                  + this%Coef_High_T_range(6)*t**3 + this%Coef_High_T_range(7)*t**4)
        else
            shcp = r*(this%Coef_Low_T_range(1)*t**(-2) + this%Coef_Low_T_range(2)*t**(-1)&
                  + this%Coef_Low_T_range(3) + this%Coef_Low_T_range(4)*t + this%Coef_Low_T_range(5)*t**2&
                  + this%Coef_Low_T_range(6)*t**3 + this%Coef_Low_T_range(7)*t**4)
        end if
    
    end function Cp
!------------------------------------------------------------------------
    !Specific heat enthalpy polynomial fitting
    pure real(rp) function H(this,t,r) result(enthalpy)
    class(GasPerfectThermoModel),intent(in)::  this
    real(rp),intent(in)::               t,r
    real(rp),parameter::                temperature_boundary=1000._rp
        
        if(t > temperature_boundary) then
            enthalpy= r*(-this%Coef_High_T_Range(1)*t**(-1) + this%Coef_High_T_Range(2)*dlog(t)&
                    + this%Coef_High_T_Range(3)*t + this%Coef_High_T_Range(4)*t**2/2._rp + this%Coef_High_T_Range(5)*t**3/3._rp&
                    + this%Coef_High_T_Range(6)*t**4/4._rp + this%Coef_High_T_Range(7)*t**5/5._rp + this%Coef_High_T_Range(8))
        else
            enthalpy= r*(-this%Coef_Low_T_Range(1)*t**(-1) + this%Coef_Low_T_Range(2)*dlog(t) + this%Coef_Low_T_Range(3)*t&
                    + this%Coef_Low_T_Range(4)*t**2/2._rp + this%Coef_Low_T_Range(5)*t**3/3._rp&
                    + this%Coef_Low_T_Range(6)*t**4/4._rp + this%Coef_Low_T_Range(7)*t**5/5._rp + this%Coef_Low_T_Range(8))
        end if

    end function H
!------------------------------------------------------------------------
    !specific entropy polynomial fitting
    pure real(rp) function S(this,t,r) result(entropy)
    class(GasPerfectThermoModel),intent(in)::  this
    real(rp),intent(in)::               t,r
    real(rp),parameter::                temperature_boundary=1000._rp
        
        if(t > temperature_boundary) then
            entropy = r*(-this%Coef_High_T_Range(1)*t**(-2)/2._rp - this%Coef_High_T_Range(2)*t**(-1)&
                    + this%Coef_High_T_Range(3)*dlog(t) + this%Coef_High_T_Range(4)*t + this%Coef_High_T_Range(5)*t**2/2._rp&
                    + this%Coef_High_T_Range(6)*t**3/3._rp + this%Coef_High_T_Range(7)*t**4/4._rp + this%Coef_High_T_Range(9))
        else
            entropy = r*(-this%Coef_Low_T_Range(1)*t**(-2)/2._rp - this%Coef_Low_T_Range(2)*t**(-1)&
                    + this%Coef_Low_T_Range(3)*dlog(t) + this%Coef_Low_T_Range(4)*t + this%Coef_Low_T_Range(5)*t**2/2._rp&
                    + this%Coef_Low_T_Range(6)*t**3/3._rp + this%Coef_Low_T_Range(7)*t**4/4._rp + this%Coef_Low_T_Range(9))
        end if
        
    end function S
!------------------------------------------------------------------------
    
end module GasPerfectThermoModel_