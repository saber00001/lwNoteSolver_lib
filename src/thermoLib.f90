module thermoLib
use constants
implicit none

    private
    public:: Sutherlandmiu
        

    
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
    
    
    !--
    elemental real(rp) function idealEquationOfStateT(rho,p,R) result(T)
    real(rp),intent(in)::   rho,p,R
        T = p / rho / R
    end function idealEquationOfStateT
    
    
end module thermoLib