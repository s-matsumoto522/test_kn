module p_error
    implicit none
    double precision, parameter :: pi = acos(-1.0d0)
    integer, parameter :: NXmin = 1, NXmax = 160     !x方向の計算領域の形状
    integer, parameter :: NYmin = 1, NYmax = 160     !y方向の計算領域の形状
    integer, parameter :: NZmin = 1, NZmax = 160     !z方向の計算領域の形状
    double precision, parameter :: Xmax = 2.0d0*pi, Ymax = 2.0d0*pi, Zmax = 2.0d0*pi  !各方向の計算領域の最大値
    double precision, parameter :: dt = 1.0d-2      !時間刻み幅
    double precision, parameter :: Re = 10d0        !レイノルズ数
    integer, save :: Ng                             !格子点数
    double precision, save :: dX, dY, dZ            !各方向の刻み幅
    double precision, save :: ddX, ddY, ddZ    !各方向の刻み幅の逆数
    double precision, save :: ddX2, ddY2, ddZ2 !各方向の刻み幅の逆数の二乗
    double precision, save :: Xmin, Ymin, Zmin      !各方向の計算領域の最小値
contains
!********************************
!   格子を設定するサブルーチン  *
!********************************
subroutine set_grid(X, Y, Z)
    double precision, intent(out) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    integer iX, iY, iZ
    !---各方向の刻み幅を計算---
    dX = Xmax / dble(NXmax)
    dY = Ymax / dble(NYmax)
    dZ = Zmax / dble(NZmax)
    !---計算用数値の設定---
    ddX = 1.0d0 / dX
    ddY = 1.0d0 / dY
    ddZ = 1.0d0 / dZ
    ddX2 = 1.0d0 / (dX**2)
    ddY2 = 1.0d0 / (dY**2)
    ddZ2 = 1.0d0 / (dZ**2)
    Ng = (NXmax-NXmin+1)*(NYmax-NYmin+1)*(NZmax-NZmin+1)
    !---各方向の計算領域の最小値を計算---
    Xmin = dX*dble(NXmin)
    Ymin = dY*dble(NYmin)
    Zmin = dZ*dble(NZmin)
    !---格子点の設定---
    do iZ = NZmin-1, NZmax+1
        do iY = NYmin-1, NYmax+1
            do iX = NXmin-1, NXmax+1
                X(iX, iY, iZ) = dX*dble(iX)
                Y(iX, iY, iZ) = dY*dble(iY)
                Z(iX, iY, iZ) = dZ*dble(iZ)
            enddo
        enddo
    enddo
end subroutine set_grid
!********************************
!   各格子点における圧力の計算  *
!********************************
subroutine set_pressure(p, X, Y, Z)
    double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    integer iX, iY, iZ
    do iZ = NZmin-1, NZmax+1
        do iY = NYmin-1, NZmax+1
            do iX = NXmin-1, NZmax+1
                p(iX,iY,iZ) = cos(2.0d0*(X(iX,iY,iZ) + Y(iX,iY,iZ) + Z(iX,iY,iZ)))
            enddo
        enddo
    enddo
end subroutine set_pressure
!****************************************
!   圧力項をx方向に微分するサブルーチン *
!****************************************
subroutine dif_pressure_x(p, dpx)
    double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: dpx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    integer iX, iY, iZ
    do iZ = NZmin-1, NZmax
        do iY = NYmin-1, NYmax
            do iX = NXmin-1, NXmax
                dpx(iX, iY, iZ) = ddX*(-p(iX, iY, iZ) + p(iX+1, iY, iZ))
            enddo
        enddo
    enddo
end subroutine dif_pressure_x
!****************************************
!   圧力項をy方向に微分するサブルーチン *
!****************************************
subroutine dif_pressure_y(p, dpy)
    double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: dpy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    integer iX, iY, iZ
    do iZ = NZmin-1, NZmax
        do iY = NYmin-1, NYmax
            do iX = NXmin-1, NXmax
                dpy(iX, iY, iZ) = ddY*(-p(iX, iY, iZ) + p(iX, iY+1, iZ))
            enddo
        enddo
    enddo
end subroutine dif_pressure_y
!****************************************
!   圧力項をz方向に微分するサブルーチン *
!****************************************
subroutine dif_pressure_z(p, dpz)
    double precision, intent(in) :: p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision, intent(out) :: dpz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    integer iX, iY, iZ
    do iZ = NZmin-1, NZmax
        do iY = NYmin-1, NYmax
            do iX = NXmin-1, NXmax
                dpz(iX, iY, iZ) = ddZ*(-p(iX, iY, iZ) + p(iX, iY, iZ+1))
            enddo
        enddo
    enddo
end subroutine dif_pressure_z
!********************************************
!   圧力項の解析解を計算するサブルーチン    *
!********************************************
    subroutine th_pressure(X, Y, Z, dpx_th, dpy_th, dpz_th)
        double precision, intent(in) :: X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(in) :: Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
        double precision, intent(out) :: dpx_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: dpy_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: dpz_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        integer iX, iY, iZ
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    dpx_th(iX, iY, iZ) = -2.0d0*sin(2.0d0*((X(iX,iY,iZ) + 0.5d0*dX) + Y(iX,iY,iZ) + Z(iX,iY,iZ)))
                    dpy_th(iX, iY, iZ) = -2.0d0*sin(2.0d0*(X(iX,iY,iZ) + (Y(iX,iY,iZ) + 0.5d0*dY) + Z(iX,iY,iZ)))
                    dpz_th(iX, iY, iZ) = -2.0d0*sin(2.0d0*(X(iX,iY,iZ) + Y(iX,iY,iZ) + (Z(iX,iY,iZ) + 0.5d0*dZ)))
                enddo
            enddo
        enddo
    end subroutine th_pressure
!********************************
!   誤差を計算するサブルーチン  *
!********************************
    subroutine cal_error(dpx, dpy, dpz, dpx_th, dpy_th, dpz_th, errorx, errory, errorz)
        double precision, intent(in) :: dpx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: dpy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: dpz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: dpx_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: dpy_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(in) :: dpz_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
        double precision, intent(out) :: errorx, errory, errorz
        integer iX, iY, iZ
        double precision errorx_norm, errory_norm, errorz_norm
        errorx_norm = 0.0d0
        errory_norm = 0.0d0
        errorz_norm = 0.0d0
        do iZ = NZmin-1, NZmax
            do iY = NYmin-1, NYmax
                do iX = NXmin-1, NXmax
                    errorx_norm = errorx_norm + (dpx(iX,iY,iZ) - dpx_th(iX,iY,iZ))**2
                    errory_norm = errory_norm + (dpy(iX,iY,iZ) - dpy_th(iX,iY,iZ))**2
                    errorz_norm = errorz_norm + (dpz(iX,iY,iZ) - dpz_th(iX,iY,iZ))**2
                enddo
            enddo
        enddo
        errorx = sqrt(errorx_norm / Ng)
        errory = sqrt(errory_norm / Ng)
        errorz = sqrt(errorz_norm / Ng)
    end subroutine cal_error
end module p_error

program main
    use p_error
    implicit none
    double precision X(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Y(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision Z(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision p(NXmin-1:NXmax+1, NYmin-1:NYmax+1, NZmin-1:NZmax+1)
    double precision dpx(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision dpy(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision dpz(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision dpx_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision dpy_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision dpz_th(NXmin-1:NXmax, NYmin-1:NYmax, NZmin-1:NZmax)
    double precision errorx, errory, errorz
    call set_grid(X, Y, Z)
    call set_pressure(p, X, Y, Z)
    call dif_pressure_x(p, dpx)
    call dif_pressure_y(p, dpy)
    call dif_pressure_z(p, dpz)
    call th_pressure(X, Y, Z, dpx_th, dpy_th, dpz_th)
    call cal_error(dpx, dpy, dpz, dpx_th, dpy_th, dpz_th, errorx, errory, errorz)
    open(11, file = 'chk_pressure_precision.dat', position = 'append')
    write(11, *) dX, dY, dZ, errorx, errory, errorz
    close(11)
end program main