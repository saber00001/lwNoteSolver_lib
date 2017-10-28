!physical model branches into a extreme scale with various correlations. 
!create an nested object to organize these models is a smarter way to refactor
module thermoLib
use constants
implicit none

    private
    public:: Sutherlandmiu
    public:: ArrheniusK
    public:: isoEntropicPressureRatio
    public:: isoEntropicTemperatureRatio
    !--
    public:: stateEquation
    
    
!--------------------------------------------
    type stateEquation
    contains
        procedure,nopass::     idealGasT
    end type stateEquation
    
    
!-------------------------------------------------------------------
    interface SutherlandMiu
        procedure:: SutherlandMiu_custom
        procedure:: SutherlandMiu_air
    end interface SutherlandMiu
    

contains


    !--
    elemental real(rp) function SutherlandMiu_custom(T,T0,T1,Miu0) result(miu)
    real(rp),intent(in)::   T,T0,T1,Miu0
        miu = (T0 + T1) / (T + T1)
        miu = miu * Miu0 * sqrt( (T / T0)**3 )
    end function SutherlandMiu_custom
    
    !--
    elemental real(rp) function SutherlandMiu_air(T) result(miu)
    real(rp),intent(in)::   T
        miu = SutherlandMiu(T,273.15_rp,110.4_rp,1.7161_rp)
    end function SutherlandMiu_air
    
!-------------------------------------------------------------------------
    elemental real(rp) function ArrheniusK(T,Ea,a,b,R) result(k)
    real(rp),intent(in)::   T,Ea,a,b,R
        k = a * T**b * exp(-Ea/R/T)
    end function ArrheniusK
    
!------------------------------------------------------------------------
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

!-------------------------------------------------------------------------
    elemental real(rp) function idealGasT(rho,p,R) result(T)
    real(rp),intent(in)::   rho,p,R
        T = p / rho / R
    end function idealGasT
    
    
end module thermoLib