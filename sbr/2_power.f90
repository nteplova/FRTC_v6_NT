module power
    use kind_module
    implicit none

    real(wp) :: psum4
    !!common /vvv2/ psum4
    real(wp) :: plost, pnab
    !!common /a0a4/ plost,pnab
    
contains

subroutine refresh_vzmax_vzmin(v, i)
    !! refresh vzmax and vzmin
    !! This needs for protect grid collapse
    use constants, only: clt
    use plasma, only: cltn
    use nr_grid, only: vzmin, vzmax
    implicit none
    real(wp), intent(in) :: v
    integer,  intent(in) :: i
    if(v.lt.vzmin(i)) vzmin(i)=v
    if(v.gt.vzmax(i)) vzmax(i)=v
end subroutine

subroutine dfind(j, i, v, powpr, pil,pic,pia,df,decv,refr,vlf,vrt,ifast)
    use constants, only: clt, zero
    use plasma, only: cltn, zza, vk, valfa, vperp ! def vperp(50,100) ????
    use rt_parameters, only: pchm, itend0, kv
    use nr_grid, only: ppv1, ppv2
    use nr_grid, only: dncount, dq1, dq2, dql
    use nr_grid, only: pdl, pdc, pda, pdfast, eta, dens, dqi0, fcoll
    implicit none
    integer,  intent(in) :: i, j, ifast
    real(wp), intent(in) :: v, powpr, pil, pic, pia, df, decv, refr, vlf, vrt
    integer k
    real(wp) :: pchgl, pchgc, pchga, denom, powlandau, powdamped
    real(wp) :: fff, dd, domin, parn, dvz, dnpar, weight, addd
    real(wp) :: arg, hevis, adda
    real(wp),parameter :: absorption_tiny = 1.d-20

    if(v.gt.cltn) return
    if(pil.gt.zero) then
        call refresh_vzmax_vzmin(v, j)
    end if
    pchgl = zero
    pchgc = zero
    pchga = zero
    denom = pil + pic + pia
    powlandau = 1.d0 - exp(-2.d0*pil)
    powdamped = 1.d0 - exp(-2.d0*denom)
    domin = powpr * powdamped
    if(denom.ne.zero) then
        !!       pchgl=powpr*(1.d0-dexp(-2d0*pil))
        !!       pchgc=powpr*dexp(-2d0*pil)*dabs(-2d0*pic)
        !!       pchga=powpr*dexp(-2d0*pil)*dabs(-2d0*pia)
        fff = domin/denom
        pchgl = abs(pil*fff)
        pchgc = abs(pic*fff)
        pchga = abs(pia*fff)
    end if
    dd=zero
    if(pil.gt.zero) then  !Landau absorption
        if(powlandau.gt.pchm) then !strong absorption
            ppv1=ppv1+pchgl
            if(dabs(df).gt.absorption_tiny) then
                dd = abs(-pchgl/df/vk(j)) * 1.d-10
                dncount(i,j)=dncount(i,j)+1.d0
            else
                dd=zero
            end if
            dq1(i,j)=dq1(i,j)+dd
        else  ! weak absorption
            ppv2=ppv2+pchgl
            dd = abs(2.d0*decv*powpr/vk(j)) * 1.d-10
            dncount(i,j)=dncount(i,j)+1.d0
            dq2(i,j)=dq2(i,j)+dd
        end if
    endif
    
    dql(i,j)=dql(i,j)+dd
    pdl(j)=pdl(j)+pchgl
    pdc(j)=pdc(j)+pchgc
    pda(j)=pda(j)+pchga

    if(ifast.eq.-1) pdfast(j)=pdfast(j)+pchgl+pchgc+pchga
    if(itend0.gt.0) then
        parn=cltn/v
        dvz=vrt-vlf
        dnpar=cltn*dvz/v**2
        weight=(refr**2-eta(j))**2/(refr**2*parn**3)
        !!!        adde=zze*(dd/dens(j))*weight
        !!!        e2perp(i,j)=e2perp(i,j)+adde
        addd=zza*(dd/dens(j))*weight/fcoll(j)/refr**3
        arg=clt/(refr*valfa)
        do k=1,kv
            if(vperp(k,j).gt.arg) then
                hevis=dsqrt((vperp(k,j)-arg)*(vperp(k,j)+arg))
                adda=addd*hevis
                dqi0(k,j)=dqi0(k,j)+adda*dnpar
            end if
        end do
     end if
     return
     end    
end module power

