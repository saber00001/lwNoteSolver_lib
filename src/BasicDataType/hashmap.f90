!for now i just need a integer hashmap
!i will refer to:
!https://en.wikipedia.org/wiki/List_of_hash_functions
!src\OpenFOAM\primitives\hashes\Hasher // 
!https://github.com/vkandy/jenkins-hash-java/blob/master/src/JenkinsHash.java
!https://www.cnblogs.com/napoleon_liu/archive/2010/12/29/1920839.html
!https://github.com/nncarlson/petaca/issues/2
!https://www.cnblogs.com/napoleon_liu/archive/2010/12/29/1920839.html
!
!some tisps from fortran on bit: no unsigned type in fortran
!b'1' -> binary constants
!o'7' -> octal constants
!z'f' -> Hexadecimal Constants
!n#number -> different interpreted number
!=> 16#number == z'number'
!--------------------------------------------------------------------
!some bit operations
!iand   ==  &
!ior    ==  |
!ieor   ==  ^
!ishft  ==  <<
!--
!func<same_type_as>
module hashmap_
use constants
use arrayOpsLib
implicit none

    private
    public:: hashmap
    public:: hashfunc
    
    !----------------------------------------------------------------------------------------------------
    integer(isp),parameter::        MapSize(28) = (/                    &
            17,         37,         79,         163,        331,        &
            673,        1361,       2729,       5471,       10949,      &
            21911,      43853,      87719,      175447,     350899,     &
            701819,     1403641,    2807303,    5614657,    11229331,   &
            22458671,   44917381,   89834777,   179669557,  359339171,  &
            718678369,  1437356741, 2147483647/)
            
    integer(isp),parameter::        dByte = sizeof(0_ip)
    
!--------------------------------------------------------------------
    type:: hashEntry
        private
        !any data can be transfered to the integer by byte
        !only character and integer is considerred, integer supports the complex structure
        !and useful data is alway a big real array
        integer(ip),dimension(:),allocatable::  key_
        integer(ip),dimension(:),allocatable::  val_
        class(hashEntry),pointer::  next_   => null()
    contains
        !--
        procedure::                 depth   => depth_Entry
        !--
        procedure::                 delet   => delet_Entry
        !--
        procedure::                 clear   => clear_Entry
        !--
        procedure::                 set     => set_Entry
        !--
        procedure::                 get     => get_Entry
        !--
        procedure::                 iseq    => iseq_Entry
    end type hashEntry
    
    
!--------------------------------------------------------------------
    type:: hashmap
        private
        class(hashEntry),dimension(:),allocatable::  buckets_
    contains
        generic::                   init    => init_n
        procedure,private::         init_n
        !-----------
        procedure::                 selectMapSize
        procedure::                 mapsize => mapsize_f
        procedure::                 bid
        procedure::                 maxcol
        procedure::                 ncol
        !procedure::                set
    end type hashmap

    
!--------------------------------------------------------------------
    interface hashfunc
        procedure:: javahash
        procedure:: bkdrhash
        procedure:: b3hs_hash_key_jenkins
    end interface hashfunc
!--------------------------------------------------------------------
    
    
    
contains
    
!-------------------------------------------------------------------------
!space for hash entry
!-------------------------------------------------------------------------
    integer(ip) function depth_Entry(this) result(depth)
    class(hashEntry),target,intent(in)::    this
    class(hashEntry),pointer::              ptr
        depth = 1;  ptr => this
        do while(associated(ptr%next_))
            depth = depth + 1
            ptr => ptr%next_
        enddo
    end function depth_Entry
    
    !--
    pure recursive subroutine set_Entry(this,key,val)
    class(hashEntry),intent(inout)::        this
    integer(ip),dimension(:),intent(in)::   key,val
        if(.not.allocated(this%key_)) then
            allocate(this%key_,source=key)
            allocate(this%val_,source=val)
        else
            if(all(this%key_==key)) then
                deallocate(this%val_)
                allocate(this%val_,source=val)
            else
                if(.not.associated(this%next_)) allocate(this%next_)
                call this%next_%set(key,val)
            endif
        endif
    end subroutine set_Entry
    
    !--
    recursive function get_Entry(this,key) result(val)
    class(hashEntry),target,intent(in)::    this
    integer(ip),dimension(:),intent(in)::   key
    integer(ip),dimension(:),pointer::      val
        if(.not.allocated(this%key_)) then
            val => null()
        elseif(all(this%key_==key)) then
            val => this%val_
        elseif(associated(this%next_)) then
            val => this%next_%get(key)
        endif
    end function get_Entry
    
    !--not allocated key, no next
    pure recursive subroutine delet_Entry(this,key,previous)
    class(hashEntry),intent(inout)::            this
    integer(ip),dimension(:),intent(in)::       key
    class(hashEntry),optional,intent(inout)::   previous
    class(hashEntry),pointer::                  tp
        !means the origrin Entry in Buckets
        if(.not.allocated(this%key_)) return
        if(all(this%key_ == key)) then  !need to delet
            if(present(previous)) then  !means it's the son list
                tp => previous%next_
                previous%next_ => this%next_
                deallocate(tp)
            else! means origin Entry in Buckets
                if(associated(this%next_)) then
                    tp => this%next_
                    this%key_   = tp%key_
                    this%val_   = tp%val_
                    this%next_  => tp%next_
                    deallocate(tp)
                else
                    deallocate(this%key_,this%val_)
                endif
            endif
        else
            if(associated(this%next_)) call this%next_%delet(key,this)
        endif
    end subroutine delet_Entry
    
    !--
    pure recursive subroutine clear_Entry(this)
    class(hashEntry),intent(inout)::    this
        if(associated(this%next_)) then
            call this%next_%clear
            deallocate(this%next_)
            this%next_ => null()
        endif
        deallocate(this%key_,this%val_)
    end subroutine clear_Entry
    
    !--
    pure logical(lp) function iseq_Entry(this,key,val) result(is)
    class(hashEntry),intent(in)::           this
    integer(ip),dimension(:),intent(in)::   key,val
        is = .false.
        if(.not.allocated(this%key_))   return
        if(size(key)/=size(this%key_))  return
        if(size(val)/=size(this%val_))  return
        if(.not.all(key==this%key_))    return
        if(.not.all(val==this%val_))    return
        is = .true.
    end function iseq_Entry

    
    
!------------------------------------------------------------------------
!space for class hashmap
!--------------------------------------------------------------------------
    pure subroutine init_n(this,n)
    class(hashmap),intent(out)::    this
    integer(ip),intent(in)::        n
        allocate(this%buckets_( this%selectMapSize(n)))
    end subroutine init_n

    !---------------------------------
    pure integer(ip) function mapsize_f(this)
    class(hashmap),intent(in)::     this
        mapsize_f = size(this%buckets_)    
    end function mapsize_f
    
    pure integer(ip) function selectMapSize(this,n)
    class(hashmap),intent(in)::     this
    integer(ip),intent(in)::        n
    integer(ip)::                   i
        do i=1,28
            if(mapsize(i)>n) exit
        enddo
        selectMapSize = mapsize(i)
    end function selectMapSize
    
    !-----------------------------------
    pure integer(ip) function bid(this,key)
    class(hashmap),intent(in)::     this
    integer(isp),intent(in)::       key
        !no matter about the sign of key, see <modulo>
        bid = modulo(hashfunc(key),int(this%mapsize(),kind=isp)) + 1_ip
    end function bid
    
    pure integer(ip) function maxcol(this) result(mc)
    class(hashmap),intent(in)::     this
    integer(ip)::                   i,ic
        mc = 0_ip
        do i=1_ip,this%mapsize()
            ic = this%buckets_(i)%depth() - 1_ip
            if(ic>mc) mc = ic
        enddo
    end function maxcol
    
    pure integer(ip) function ncol(this) result(nc)
    class(hashmap),intent(in)::     this
    integer(ip)::                   i,ic
        nc = 0_ip
        do i=1_ip,this%mapsize()
            ic = this%buckets_(i)%depth() - 1_ip
            nc = nc + ic
        enddo
    end function ncol
    !pure subroutine set(this,key,val)
    !
    !
    !end subroutine set
    
    
    
    
    
    
    
    
    
    
    
    
    

!--------------------------------------------------------------
!space for different hash function
!---------------------------------------------------------------
    !https://www.zhihu.com/question/51784530
    pure integer(isp) function javahash(key)
    integer(ip),intent(in)::            key
    character(len=sizeof(key),kind=1):: ckey
        ckey = transfer(key,mold=ckey)
        javahash = hashfunc(adjustl(ckey))
        javahash = ieor(key,ieor(ishft(key,-20),ishft(key,-12)))
        javahash = ieor(ieor(javahash,ishft(javahash,-7)),ishft(javahash,-4))
    end function javahash
    
    !https://www.zhihu.com/question/20507188
    pure integer(isp) function bkdrhash(key)
    character(*),intent(in)::   key
    integer(isp),parameter::    s = 131
    integer(isp)::              n,bt,i
        n = len_trim(key)
        do i=1,n
            bt = ichar(key(i:i))
            bkdrhash = bkdrhash*s + bt
        enddo
        bkdrhash = iand(bkdrhash,z'7fffffff')
    end function bkdrhash
    
!---------------------------------------------------------------
    !http://computer-programming-forum.com/49-fortran/0596e59d0fa2e5e4.htm
    pure function b3hs_hash_key_jenkins (key, range) result (code) 
    character(*), intent(in)    :: key 
    integer(isp), intent(in)    :: range 
    integer(isp)                :: code 
    integer(isp)                :: len_key 
    integer(isp)                :: a 
    integer(isp)                :: b 
    integer(isp)                :: c 
    integer(isp)                :: k 
! hash the key into a code, using the algorithm 
! described by bob jenkins at: 
! http://burtleburtle.net/bob/hash/doobs.html 
! 
! note that range should be a power of 2, and 
! that the 32-bit algorithm is used
        len_key = len_trim(key)
        a = z'9e3779b9'
        b = a 
        c = z'12345678'
        k = 1 
        char_loop : do while(len_key>=12) 
! pack the key into 32 bits 
            a = a + ichar(key(k+0:k+0))  + ishft(ichar(key(k+1:k+1)), 8) + & 
            &       ishft(ichar(key(k+2:k+2)), 16) + ishft(ichar(key(k+3:k+3)), 24) 
            b = b + ichar(key(k+4:k+4))  + ishft(ichar(key(k+5:k+5)), 8) + & 
            &       ishft(ichar(key(k+6:k+6)), 16) + ishft(ichar(key(k+7:k+7)), 24) 
            c = c + ichar(key(k+8:k+8))  + ishft(ichar(key(k+9:k+9)), 8) + & 
            &       ishft(ichar(key(k+10:k+10)), 16) + ishft(ichar(key(k+11:k+11)), 24) 
! mix it up 
            call b3hs_hash_key_jenkins_mix_(a,b,c)
            k = k + 12 
            len_key = len_key - 12 
        end do char_loop 
        c = c + len_key
        
! process remaining bits
        select case(len_key)
        case(11) 
            c = c + ishft(ichar(key(k+10:k+10)), 24) + ishft(ichar(key(k+9:k+9)), 16) + & 
            &       ishft(ichar(key(k+8:k+8)), 8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(10) 
            c = c + ishft(ichar(key(k+9:k+9)), 16) + ishft(ichar(key(k+8:k+8)), 8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(9) 
            c = c + ishft(ichar(key(k+8:k+8)), 8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(8) 
            b = b + ishft(ichar(key(k+7:k+7)), 24) + ishft(ichar(key(k+6:k+6)), 16) + & 
            &       ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(7) 
            b = b + ishft(ichar(key(k+6:k+6)), 16) + ishft(ichar(key(k+5:k+5)), 8) + & 
            &       ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(6) 
            b = b + ishft(ichar(key(k+5:k+5)), 8) + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(5) 
            b = b + ichar(key(k+4:k+4)) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(4) 
            a = a + ishft(ichar(key(k+3:k+3)), 24) + ishft(ichar(key(k+2:k+2)), 16) + & 
            &       ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(3) 
            a = a + ishft(ichar(key(k+2:k+2)), 16) + ishft(ichar(key(k+1:k+1)), 8) + & 
            &       ichar(key(k:k)) 
        case(2) 
            a = a + ishft(ichar(key(k+1:k+1)), 8) + ichar(key(k:k)) 
        case(1) 
            a = a + ichar(key(k:k))
        end select 
        call b3hs_hash_key_jenkins_mix_(a,b,c) 
        code = iand(c, range - 1) + 1

    contains
        pure subroutine b3hs_hash_key_jenkins_mix_(a,b,c)
        integer(isp),intent(inout)::    a,b,c 
            a = ieor(a - b - c, ishft(c, -13))
            b = ieor(b - c - a, ishft(a, 8)) 
            c = ieor(c - a - b, ishft(b, -13)) 
            a = ieor(a - b - c, ishft(c, -12)) 
            b = ieor(b - c - a, ishft(a, 16)) 
            c = ieor(c - a - b, ishft(b, -5)) 
            a = ieor(a - b - c, ishft(c, -3)) 
            b = ieor(b - c - a, ishft(a, 10)) 
            c = ieor(c - a - b, ishft(b, -15)) 
        end subroutine b3hs_hash_key_jenkins_mix_
    end function b3hs_hash_key_jenkins 

end module hashmap_