module stringOpsLib
use constants
implicit none

    private
    !--
    public:: readKeyVal
    public:: countsubstring,lowerString,upperString
    
    
    
    
    !--------------------------------------------------
    interface readKeyVal
        procedure::  readkeyival
        procedure::  readkeyiar1val
        procedure::  readkeyiar2val
        procedure::  readkeycval
        procedure::  readkeyrval
    end interface readKeyVal
    
contains

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


end module stringOpsLib