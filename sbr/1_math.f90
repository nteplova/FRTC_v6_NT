module math_module
    use kind_module
    implicit none
    
    abstract interface
          function I_func(ptet,pa) result(f)
            !! user function f(theta, pa)
            import :: wp
            implicit none
            real(wp), intent(in)    :: ptet
            real(wp), intent(in)    :: pa
            real(wp)                :: f
        end function 
    end  interface

contains
    subroutine diff(x,y,n,dy)
        implicit real*8 (a-h,o-z)
        dimension y(*),x(*),dy(*)
        integer :: k, n
        dy(1)=(y(2)-y(1))/(x(2)-x(1))
        do k=2,n-1
            dy(k)=(y(k+1)-y(k-1))/(x(k+1)-x(k-1))
        end do
            dy(n)=(y(n)-y(n-1))/(x(n)-x(n-1))
        return
    end

    subroutine integral(ibeg,iend,x,y,fout)
        implicit real*8 (a-h,o-z)
        integer :: ibeg, iend
        dimension x(*),y(*)
        integer :: i, n1, n2, ie
        fout=0.d0
        if(ibeg.eq.iend) return
        znak=1.d0
        n1=ibeg
        n2=iend
        if(n2.lt.n1) then
            znak=-1.d0
            ie=n1
            n1=n2
            n2=ie
        end if
        sum=0.d0
        do i=n1+1,n2
            dx=x(i)-x(i-1)
            dsum=y(i)+y(i-1)
            sum=sum+.5d0*dsum*dx
        end do
        fout=znak*sum
    end    

    subroutine fsmoth4(x,y,n,ys)
        use constants, only : zero
        use approximation
        implicit real*8 (a-h,o-z)
        integer, intent(in) :: n
        !external polin2
        integer, parameter :: np=10, imax=601
        integer, parameter :: m0=1, ndp=1
        ! m0,ndp - parameters of smoothing procedure
        dimension y(n),x(n),ys(n)
        dimension yy(imax),xx(imax)
        dimension coeffs(np),cffs(np)
        dimension dys(imax)
        integer :: i, j, k, id
        integer :: m, m2, nmax, jlast
        if(n.gt.imax) stop 'small imax in subroutine fsmoth4()'
        call diff(x,y,n,dys)
        do k=1,n
            ys(k)=y(k)
        end do
        m=m0
        m2=m+2
        id=m+ndp
        nmax=n-id
        xs=x(1)
        do j=1,nmax
            do i=1,id
                xx(i)=x(j+i)-xs
                yy(i)=y(j+i)-ys(j)-dys(j)*xx(i)
            end do
            call approx(xx,yy,id,polin2,m,coeffs)
            cffs(1)=ys(j)
            cffs(2)=dys(j)
            do k=1,m
                cffs(k+2)=coeffs(k)
            end do
            xs=x(j+1)
            ys(j+1)=fdf(xx(1),cffs,m2,dys(j+1))
        end do
  
        j=nmax+1
  1     continue
        jlast=j
        id=n-jlast
        m=id-1
        if(m.eq.0) then
            j=j+1
            xs=x(j)
            ys(j)=fdf(xx(2),cffs,m2,dys(j))
            return
        end if
        m2=m+2
        do i=1,id
            xx(i)=x(jlast+i)-xs
            yy(i)=y(jlast+i)-ys(jlast)-dys(jlast)*xx(i)
        end do
        call approx(xx,yy,id,polin2,m,coeffs)
        cffs(1)=ys(jlast)
        cffs(2)=dys(jlast)
        do k=1,m
            cffs(k+2)=coeffs(k)
        end do
        j=j+1
        xs=x(j)
        ys(j)=fdf(xx(1),cffs,m2,dys(j))
        go to 1
    end    



    real(wp) function gaussint(f,a,b,r,eps)
    !! Compute integral of f using 12-point Gauss quadrature
        implicit none
        !real(wp), intent(in) :: f
        procedure(I_func) :: f
        !! f(ptet,pa)
        real(wp), intent(in) :: a, b
        !! limits of integration
        real(wp), intent(in) :: r, eps
        !!
        real(wp) w(12), x(12)
        real(wp) aa, bb, c1, c2, s8, s16, u, y, delta
        integer i
        !!      save w,x,const !sav#
        real(wp), parameter :: const = 1.0d-12
        data w &
        /0.101228536290376, 0.222381034453374, 0.313706645877887, &
         0.362683783378362, 0.027152459411754, 0.062253523938648, &
         0.095158511682493, 0.124628971255534, 0.149595988816577, &
         0.169156519395003, 0.182603415044924, 0.189450610455069/
        data x &
        /0.960289856497536, 0.796666477413627, 0.525532409916329, &
         0.183434642495650, 0.989400934991650, 0.944575023073233, &
         0.865631202387832, 0.755404408355003, 0.617876244402644, &
         0.458016777657227, 0.281603550779259, 0.095012509837637/
        delta=const*dabs(a-b)
        gaussint=0d0
        aa=a
  5     y=b-aa
        if (dabs(y).le.delta) return
  2     bb=aa+y
        c1=0.5d0*(aa+bb)
        c2=c1-aa
        s8=0d0
        s16=0d0
        do i = 1,4
            u=x(i)*c2
            s8=s8+w(i)*(f(c1+u,r)+f(c1-u,r))
        end do 
        do i = 5,12
            u=x(i)*c2
            s16=s16+w(i)*(f(c1+u,r)+f(c1-u,r))
        end do
        s8=s8*c2
        s16=s16*c2
        if(dabs(s16-s8) .gt. eps*(1d0+dabs(s16))) go to 4
        gaussint=gaussint+s16
        aa=bb
        go to 5
  4     y=0.5d0*y
        if(dabs(y) .gt. delta) go to 2
        write(*,7)
        gaussint=0d0
        return
  7     format(1x,'gaussint ... too high accuracy required')
    end

end module math_module