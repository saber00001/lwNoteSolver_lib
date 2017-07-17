!refer to h ttps://gcc.gnu.org/onlinedocs/gfortran/STRUCTURE-and-RECORD.html
!period(.) is an old component access, and we use percent(%)
module constants
use,intrinsic:: ieee_arithmetic
implicit none
    
    integer,parameter ::            isp = selected_int_kind(9)
    integer,parameter ::            idp = selected_int_kind(13)
    integer,parameter ::            rsp = selected_real_kind(p=6,r=37)
    integer,parameter ::            rdp = selected_real_kind(p=15,r=307)
    
    !this is standard(like *) output of different types, as a reference here
    character(10),  parameter::     fmt_rdp  = '(E23.15E3)'
    character(9),   parameter::     fmt_rsp  = '(E13.6E2)'
    character(5),   parameter::     fmt_idp  = '(I20)'
    character(5),   parameter::     fmt_isp  = '(I11)'
    
    integer,parameter::             ip  = isp
    integer,parameter::             rp  = rdp
    integer,parameter::             lp  = isp
    integer,parameter::             cl  = 32
    real(rp),parameter::            zero = 0.d0
    
    integer,parameter::             lowercase_a = ichar('a')
    integer,parameter::             lowercase_z = ichar('z')
    integer,parameter::             uppercase_a = ichar('A')
    integer,parameter::             uppercase_z = ichar('Z')


    
    !--------------------------------physical and mathmatic constant    
    real(rdp),parameter::           pi      = 4.d0*atan(1.d0)
    real(rdp),parameter::           spi     = sqrt(pi)
    real(rdp),parameter::           pi2     = pi**2
    real(rdp),parameter::           e       = dexp(1.d0)
    real(rdp),parameter::           k_b     = 1.38067852d-23    !boltzmann constant
    real(rdp),parameter::           R_c     = 8.3144598d0       !gas constant   ( J * [K^-1] * [mol^-1] )
    real(rdp),parameter::           R_air   = 287.058d0         !specific gas constant for ideal gas ( J * [kg^-1] * [mol^-1] )
    real(rdp),parameter::           gm_diatomic = 1.4d0         !specific gas
    
    !--------------------------------numeric constant
    integer(isp),parameter::        Minisp  = -2147483648_isp
    integer(idp),parameter::        Minidp  = -9223372036854775808_idp
    real(rdp),parameter::           minrdp  = transfer(-1_rdp,0._rdp)
    real(rsp),parameter::           minrsp  = transfer(-1_rsp,0._rsp)
    
    
    !--------------------------------computing control parameters
    real(rdp),parameter::           gradLimitingk           = 1.2d0
    real(rdp),parameter::           GlobalEps               = 1.d-12
    
    
    !---here we offer two methods for disabling program [disablesub]&[disablenumber]
!----------------------------------
    !a special ==, because of its unlimited polymorphism, don't override the intrinsic ==
    interface operator(.lweq.)
        procedure:: anyiseq
    end interface
    
    !don't use module procedure refer to
    !h ttps://software.intel.com/en-us/forums/
    !intel-visual-fortran-compiler-for-windows/topic/721674
    interface disableprogram
        procedure:: disableprogram_
    end interface disableprogram

    !--
    interface disablenumber
        procedure::  rsp_nan
        procedure::  rdp_nan
        procedure::  isp_nan
        procedure::  idp_nan
    end interface disablenumber
    
    !--
    interface disabled_stat
        procedure::  disabled_stat_rsp
        procedure::  disabled_stat_rdp
        procedure::  disabled_stat_isp
        procedure::  disabled_stat_idp
    end interface disabled_stat
    

    !--this is a basic procedure type, i put it here for the common use
    abstract interface
        pure real(rp) function absf1(x) result(y)
        import:: rp
        real(rp),intent(in)::   x
        end function absf1
    end interface
    
contains

    !--compare two scalar byte by byte
    elemental logical(lp) function anyiseq(lhs,rhs) result(r)
    class(*),intent(in)::               lhs,rhs
    integer(1),dimension(sizeof(lhs)):: lb
    integer(1),dimension(sizeof(rhs)):: rb
    integer(ip)::                       i
        r   = .false.
        if(.not.same_type_as(lhs,rhs))  return
        if(.not.size(lb)==size(rb))     return
        lb  = transfer(lhs,mold=lb)
        rb  = transfer(rhs,mold=rb)
        r   = all(lb==rb)
    end function anyiseq
    
    !-------------------------------
    pure subroutine disableprogram_
    integer(ip),dimension(:),allocatable:: n
        n(1) = 0
    end subroutine disableprogram_

    elemental subroutine rsp_nan(r)
    real(rsp),intent(out)::     r
        r = minrsp
    end subroutine rsp_nan

    elemental subroutine rdp_nan(r)
    real(rdp),intent(out)::     r
        r = minrdp
    end subroutine rdp_nan
    
    elemental subroutine isp_nan(i)
    integer(isp),intent(out)::  i
        i = minisp
    end subroutine isp_nan
    
    elemental subroutine idp_nan(i)
    integer(idp),intent(out)::  i
        i = minidp
    end subroutine idp_nan
    
    
    !--
    elemental function disabled_stat_rsp(r) result(l)
    real(rsp),intent(in)::      r
    logical(lp)::               l
        l = isnan(minrsp)
    end function disabled_stat_rsp
    
    elemental function disabled_stat_rdp(r) result(l)
    real(rdp),intent(in)::      r
    logical(lp)::               l
        l = isnan(minrdp)
    end function disabled_stat_rdp
    
    elemental function disabled_stat_isp(i) result(l)
    integer(isp),intent(in)::   i
    logical(lp)::               l
        l = i == minisp
    end function disabled_stat_isp
    
    elemental function disabled_stat_idp(i) result(l)
    integer(idp),intent(in)::   i
    logical(lp)::               l
        l = i == minidp
    end function disabled_stat_idp
    
end module constants