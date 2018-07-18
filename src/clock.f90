module clock_
use constants
!dir$ if defined(_openmp)
use omp_lib,only: omp_get_wtime
!dir$ end if
implicit none
    
    private
    public:: wclock,fclock
    
    !--wall clock
    type:: wclock
    
        private
        real(rp),dimension(2)::         wtime_ = 0.d0
        
    contains
    
        procedure::                     start
        procedure::                     pauseReco
        procedure::                     wtimeCost
        procedure::                     clear
        procedure::                     printWallClock
        
    end type wclock
    
    !--field clock
    type:: fclock
        private
        real(rp)::                      dt
        real(rp)::                      cfl
    end type fclock
    
contains

    !--
    subroutine start(this)
    class(wclock),intent(inout)::       this
        !dir$ if defined(_openmp)
        this%wtime_(1) = omp_get_wtime()
        !dir$ else
        call cpu_time(this%wtime_(1))
        !dir$ end if
    end subroutine start
    
    
    !--
    subroutine pauseReco(this)
    class(wclock),intent(inout)::       this
    real(rp)::                          t
    
        !dir$ if defined(_openmp)
        t = omp_get_wtime()
        !dir$ else
        call cpu_time(t)
        !dir$ end if
        this%wtime_(2) = this%wtime_(2) + t - this%wtime_(1)
        this%wtime_(1) = 0._rp
        
    end subroutine pauseReco
    
    !--
    pure real(rp) function wtimeCost(this)
    class(wclock),intent(in)::          this
        wtimeCost = this%wtime_(2)
    end function wtimeCost
    
    !--
    pure subroutine clear(this)
    class(wclock),intent(out)::         this
    end subroutine clear
    
    !--
    subroutine printWallClock(this)
    class(wclock),intent(in)::          this
    character(9),parameter,dimension(12)::  month = [   &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' ]
    character(8)::                      ampm
    integer(ip),dimension(8)::          dati
    integer(ip)::                       y,m,d,h,n,s,mm
    
        call date_and_time(values=dati)
        y = dati(1);    m = dati(2);    d = dati(3)
        h = dati(5);    n = dati(6);    s = dati(7)
        mm = dati(8)

        if ( h < 12 ) then
            ampm = 'AM'
        else if ( h == 12 ) then
            ampm = merge('Noon','PM  ',n==0.and.s==0)
        else
            h = h - 12
            ampm = merge('PM      ',merge('Midnight','AM      ',n==0.and.s==0),h<12)
        end if

        print'(i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)',    &
            d, trim(month(m)), y, h, ':', n, ':', s, '.', mm, trim(ampm)

    end subroutine printWallClock
    
end module clock_