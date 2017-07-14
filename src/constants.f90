!refer to h ttps://gcc.gnu.org/onlinedocs/gfortran/STRUCTURE-and-RECORD.html
!period(.) is an old component access, and we use percent(%)
module constants
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
    real(rdp),parameter::           GM_diatomic = 1.4d0         !specific gas
    
    !--------------------------------numeric constant
    real(rdp),parameter::           minrdp  = -1.d50000         !-> -infinity
    real(rsp),parameter::           minrsp  = minrdp
    integer(isp),parameter::        Minisp  = -2147483648
    integer(idp),parameter::        Minidp  = -9223372036854775808
    
    
    !--------------------------------computing control parameters
    real(rdp),parameter::           gradLimitingk           = 1.2d0
    real(rdp),parameter::           IntegralDefaultTolerance= 1.d-9
    integer(ip),parameter::         Integralmaxloop         = 15
    
    
    !---here we offer two methods for disabling program [disablesub]&[disablenumber]
!----------------------------------
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
    interface inquirenumberstat
        procedure::  inquirenumberstat_rsp
        procedure::  inquirenumberstat_rdp
        procedure::  inquirenumberstat_isp
        procedure::  inquirenumberstat_idp
    end interface inquirenumberstat
    
    !--
    interface readkeyval
        procedure::  readkeyival
        procedure::  readkeyiar1val
        procedure::  readkeyiar2val
        procedure::  readkeycval
        procedure::  readkeyrval
    end interface readkeyval
    
    
    !--this is a basic procedure type, i put it here for the common use
    abstract interface
        pure function absf1(x) result(y)
        import:: rp
        real(rp),intent(in)::   x
        real(rp)::              y
        end function absf1
    end interface
    
contains

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
    elemental function inquirenumberstat_rsp(r) result(l)
    real(rsp),intent(in)::      r
    logical(lp)::               l
        l = merge(.true. , .false. , r == minrsp)
    end function inquirenumberstat_rsp
    
    elemental function inquirenumberstat_rdp(r) result(l)
    real(rdp),intent(in)::      r
    logical(lp)::               l
        l = merge(.true. , .false. , r == minrdp)
    end function inquirenumberstat_rdp
    
    elemental function inquirenumberstat_isp(i) result(l)
    integer(isp),intent(in)::   i
    logical(lp)::               l
        l = merge(.true. , .false. , i == minisp)
    end function inquirenumberstat_isp
    
    elemental function inquirenumberstat_idp(i) result(l)
    integer(idp),intent(in)::   i
    logical(lp)::               l
        l = merge(.true. , .false. , i == minidp)
    end function inquirenumberstat_idp
    
    
!-------------------------------------------------
    !--
    elemental subroutine lowerString(str)
    character(*),intent(inout)::    str
    integer::                       i,n
        do i=1,len(str)
            n = ichar(str(i:i))
            if(n>=uppercase_a.and.n<=uppercase_z) then
                str(i:i) = char ( ichar(str(i:i)) - uppercase_a + lowercase_a )
            endif
        enddo
    end subroutine lowerString
    
    !--
    elemental subroutine upperString(str)
    character(*),intent(inout)::    str
    integer::                       i,n
        do i=1,len(str)
            n = ichar(str(i:i))
            if(n>=lowercase_a.and.n<=lowercase_z) then
                str(i:i) = char ( ichar(str(i:i)) - lowercase_a + uppercase_a )
            endif
        enddo
    end subroutine upperString
    
    !--
    elemental function countsubstring(string,substring) result(n)
    character(*),intent(in)::   string,substring
    integer(ip)::               n,loc,ls,i
        n   = 0
        ls  = len(substring)
        if( ls <= 0 .or. ls > len(string)) return

        loc = 1
        do while(.true.)
            i = index(string(loc:),substring)
            if(i > 0) then
                n   = n + 1
                loc = loc + i - 1 + ls
            else
                exit
            endif
        enddo
    end function countsubstring
    
    
    
!-----------------------------------------------------------
    !--
    pure subroutine readkeycval(string,key,val)
    character(*),intent(in)::   string,key
    character(*),intent(out)::  val
    integer(ip)::               st
    
        st = index(string,key)
        
        if(st == 0 ) then 
            val = 'null'
            return
        endif
        
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            val = 'null'
            return
        endif
        st = st + 1
        
        !st = index(string(st:),'=') + st
        read(string(st:),*) val
        
    end subroutine readkeycval

    !--    
    pure subroutine readkeyrval(string,key,val)
    character(*),intent(in)::   string,key
    real(rp),intent(out)::      val
    integer(ip)::               st
        st = index(string,key)
        if(st == 0 ) then 
            call disablenumber(val)
            return
        endif
        
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disablenumber(val)
            return
        endif
        st = st + 1
        !st = index(string(st:),'=') + st
        
        read(string(st:),*) val
        
    end subroutine readkeyrval
    
    !--
    pure subroutine readkeyival(string,key,val)
    character(*),intent(in)::   string,key
    integer(ip),intent(out)::   val
    integer(ip)::               st
    
        st = index(string,key)
        if(st == 0 ) then 
            call disablenumber(val)
            return
        endif
        
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disablenumber(val)
            return
        endif
        st = st + 1
        
        !st = index(string(st:),'=') + st
        
        read(string(st:),*) val
        
    end subroutine readkeyival
    
    !--
    pure subroutine readkeyiar1val(string,key,val)
    character(*),intent(in)::   string,key
    integer(ip),dimension(:),intent(out)::   val
    integer(ip)::               st,ed
    
        st = index(string,key)
        if(st == 0 ) then 
            call disablenumber(val)
            return
        endif
        
        st = st + len(key)
        
        !check '='
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disablenumber(val)
            return
        endif
        st = st + 1

        !check '('
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='(') then
            call disablenumber(val)
            return
        end if
        st = st + 1
        
        read(string(st:),*) val
        
    end subroutine readkeyiar1val
    
    !--
    pure subroutine readkeyiar2val(string,key,val)
    character(*),intent(in)::   string,key
    integer(ip),dimension(:,:),intent(out)::   val
    integer(ip)::               st,n,i
    
        st = index(string,key)
        if(st == 0 ) then 
            call disablenumber(val)
            return
        endif
        
        st = st + len(key)
        st = verify(string(st:),' ') + st - 1
        if(string(st:st)/='=') then
            call disablenumber(val)
            return
        endif
        st = st + 1
        
        n  = size(val,dim=2)
        do i =1,n
            st = verify(string(st:),' ') + st - 1
            if(string(st:st)/='(') then
                call disablenumber(val)
                return
            end if
            st = st + 1
            read(string(st:),*) val(:,i)
            st = index(string(st:),')') + st
        enddo
        
    end subroutine readkeyiar2val
    
end module constants