!physical model branches into an extreme scale with various corrections.
!create an nested object to organize these models is a smarter way to refactor
!use the instatiated object always
module GasChemReactionModel_
use constants
use GasPerfectThermoModel_
implicit none

    !private::IdealGasComponentThermo
    !instantiation
    !public:: IdealGasComponentThermo
    private          
    public:: GasChemReactionModel
    !public:: GasComponentModel
    !--
    type GasChemReactionModel
    
    contains
        
        !chemistry production of basic model, see chemKin of Ansys
        generic::                   ProductionRate => ProductionRate_Basic
        procedure,nopass::          ProductionRate_Basic
        
        !
        procedure,nopass::          ProgressRate
        
        
        !calculate rate constant => k = a T^b \exp{- \frac{Ea}{R T}}
        generic::                   RateConstant => Arrhenius
        procedure,nopass::          Arrhenius
        
        !calculate reverse Rate Constant based on Equilibrium assumption
        generic::                   reverseRateConstant => reverseRateConstant_equilibrium!,
        procedure,nopass::          reverseRateConstant_equilibrium
        !procedure,nopass::          reverseRateConstant_Saha
        
        !calculate Equilibrium Constant based on experience formula
        procedure,nopass::          EquilibriumConstant
        
        !procedure,nopass::          Temperature
        procedure,nopass::          ThreeBodyConcentration
        
        procedure,nopass::          Temperature
        
    end type GasChemReactionModel
    !----------------------------------------------------------------
    type GasComponentModel
        
        private
            real(rp),allocatable,dimension(:)::     Coef_High_T_Range
            real(rp),allocatable,dimension(:)::     Coef_Low_T_Range
    contains  
          
        generic::                   init    =>      init_n,     &
                                                    init_n_val, &
                                                    init_vals    
        
        procedure::                                 init_n
        procedure::                                 init_n_val
        procedure::                                 init_vals
    !
    end type GasComponentModel
    
    !type(GasPerfectThermoModel)::IdealGasComponentThermo
        
contains
!-------------------------------------------------------------------------
    !-------------
    pure subroutine init_n(this,ncoef)
    class(GasComponentModel),intent(out)::      this
    integer(ip),intent(in)::                    ncoef
    
        allocate(this%Coef_High_T_Range(ncoef))
        allocate(this%Coef_Low_T_Range(ncoef))
        this%Coef_High_T_Range          = 0._rp
        this%Coef_Low_T_Range           = 0._rp
    
    end subroutine init_n
    !--
    pure subroutine init_n_val(this,ncoef,Coef_H,Coef_L)
    class(GasComponentModel),intent(out)::      this
    integer(ip),intent(in)::                    ncoef
    real(rp),intent(in)::                       Coef_H,Coef_L
    
        allocate(this%Coef_High_T_Range(ncoef))
        allocate(this%Coef_Low_T_Range(ncoef))
        this%Coef_High_T_Range          = Coef_H
        this%Coef_Low_T_Range           = Coef_L
    
    end subroutine init_n_val
    !--
    pure subroutine init_vals(this,Coef_H,Coef_L)
    class(GasComponentModel),intent(out):: this
    real(rp),dimension(:),intent(in)::  Coef_H,Coef_L
        
        allocate(this%Coef_High_T_Range,source=Coef_H)
        allocate(this%Coef_Low_T_Range,source =Coef_L)
        
    end subroutine init_vals
    
    !---------------
    
!-------------------------------------------------------------------------
    !this is a corrected version, details refer to wiki Arrhenius method
    elemental real(rp) function Arrhenius(a,b,Ea,R,T) result(k)
    real(rp),intent(in)::                                   a,b,Ea,R,T
        k = a * T**b * exp(-Ea/R/T)
    end function Arrhenius
    
    !production rate of concentraction[mol], Species concentraction[mol] = rho_i / MolecularWeight
    pure function ProductionRate_Basic(Species,stoichiometricF,stoichiometricB,&
                                                            RateConstantF,RateConstantR,&
                                                            ThreeBodyConcentration) result(pr)
    real(rp),dimension(:),intent(in)::                      Species,RateConstantF,RateConstantR,ThreeBodyConcentration
    real(rp),dimension(:,:),intent(in)::                    stoichiometricF,stoichiometricB
    real(rp),dimension(size(Species))::                     pr
    integer(ip)::                                           si
        do si=1,size(pr)
            pr(si) = sum((stoichiometricB(:,si) - stoichiometricF(:,si)) * &
                            ProgressRate( Species,stoichiometricF,stoichiometricB,&
                            RateConstantF,RateConstantR &
                            ) * ThreeBodyConcentration(:) &    
                        )
        enddo
    end function ProductionRate_Basic
    
    !rate of progress variables
    pure function ProgressRate(Species,stoichiometricF,stoichiometricB,&
                                                            RateConstantF,RateConstantR) result(pr)
    real(rp),dimension(:),intent(in)::                      Species,RateConstantF,RateConstantR
    real(rp),dimension(:,:),intent(in)::                    stoichiometricF,stoichiometricB
    real(rp),dimension(size(RateConstantF))::               pr
    integer(ip)::                                           ri        
        do ri=1,size(pr)
            pr(ri) = RateConstantF(ri) * product(Species**stoichiometricF(ri,:))
            pr(ri) = pr(ri) - RateConstantR(ri) * product(Species**stoichiometricB(ri,:))
        enddo
    end function ProgressRate

    !reverseRateConstant
    pure function reverseRateConstant_equilibrium(rateConstant,EquilibriumConstant) result(rrc)
    real(rp),dimension(:),intent(in)::                      rateConstant,EquilibriumConstant
    real(rp),dimension(size(rateConstant))::                rrc
        rrc = rateConstant / EquilibriumConstant
    end function reverseRateConstant_equilibrium
    
    !stochiometric = stochiometricb - stochiometricf | (nr,ns)
    pure function EquilibriumConstant(stoichiometric,entropy,enthalpy,R,T) result(ec)
    real(rp),dimension(:),intent(in)::                      entropy,enthalpy
    real(rp),dimension(:,:),intent(in)::                    stoichiometric
    real(rp),intent(in)::                                   R,T
    real(rp),dimension(size(stoichiometric,dim=1))::        ec
    integer(ip)::                                           ri
    real(rp),parameter::                                    Patm=P_atm * 1.e-6_rp !pa -> 1 dynes/cm^2 = 10-6 bar
        do ri = 1,size(ec)
            !equilibrium constant constant
            ec(ri) = exp(sum(stoichiometric(ri,:)*entropy/R) - sum(stoichiometric(ri,:)*enthalpy/R/T))
            !equilibrium constant
            ec(ri) = ec(ri) * (Patm / R/T)**(sum(stoichiometric(ri,:)))
        enddo
    end function EquilibriumConstant
    
    !Three body concentration = sum(molality * Three body coefficient)
    pure function ThreeBodyConcentration(MC,ThreeBodyCoefficient) result(tbc)
    real(rp),dimension(:),intent(in)::                      MC
    real(rp),dimension(:,:),intent(in)::                    ThreeBodyCoefficient
    real(rp),dimension(size(ThreeBodyCoefficient,dim=1))::  tbc
    integer(ip)::                                           re,nr
        
        nr = size(ThreeBodyCoefficient,dim=1)
        do re =1,nr
            tbc(re) = sum(ThreeBodyCoefficient(re,:) * MC(:))
            !not three body reaction
            if(sum(ThreeBodyCoefficient(re,:))==0._rp) tbc(re)=1._rp
        end do
          
    end function ThreeBodyConcentration
    
    !caculate gas temperature
    pure subroutine Temperature(t0,ee,ri,yi,Total_Coef_High_T_range,Total_Coef_Low_T_range,temp)
    real(rp),intent(in)::                                           t0,ee
    real(rp),dimension(:),intent(in)::                              ri,yi
    real(rp),dimension(:,:),intent(in)::                            Total_Coef_High_T_range,Total_Coef_Low_T_range
    real(rp),intent(out)::                                          temp
    integer(ip)::                                                   count,sp,i,j,ns
    real(rp)::                                                      t1,t2,fff,ff
    real(rp),dimension(:),allocatable::                             cpi,hi
    type(GasComponentModel),dimension(:),allocatable::              GasComponent
    type(GasPerfectThermoModel)::                                   GasPerfectThermo

        ns =size(ri)
        allocate(cpi(ns),hi(ns))
        allocate(GasComponent(ns)) 
        do i=1,ns
            call GasComponent(i)%init(Total_Coef_High_T_range(i,:),Total_Coef_Low_T_range(i,:))
        end do
        count=0
        do count=0,10
            t1=t0
            t2=t1
            do while((t2==t1))
                do sp=1,ns
                    call GasPerfectThermo%init(GasComponent(sp)%Coef_High_T_range,GasComponent(sp)%Coef_Low_T_range)
                    cpi(sp) =   GasPerfectThermo%Cp(GasPerfectThermo,t1,ri(sp))
                    hi(sp)  =   GasPerfectThermo%H(GasPerfectThermo,t1,ri(sp))  
                end do
                fff=(sum(yi(:)*hi(:))-ee*1.e-3_rp)-sum(yi(:)*ri(:))*t1
                ff=(sum(yi(:)*cpi(:)))-sum((yi(:))*ri(:))
                t2=t1-fff/ff
                if(dabs(t2-t1).gt.1.e-6_rp) then
                    t1=t2
                else
                    exit
                end if
            end do
        end do
        temp=t2
    end subroutine Temperature
    
    
end module GasChemReactionModel_