module slvp
implicit none
real,parameter::R=287.04, &
                G=9.81,   &
                GAMMA=0.0065, &
                TC=273.16+17.5, &
                PCONST=10000

real,parameter::CEXP=GAMMA*R/G
contains

subroutine compute_seaprs_from_derived(nx,ny,nz,z,t,p,q,sp,mask)
implicit none

integer,intent(in)::nx,ny,nz
real,dimension(nx,ny,nz),intent(in)::z,t,p,q
real,dimension(nx,ny),intent(out)::sp
logical,dimension(nx,ny),intent(out)::mask

integer::i,j,k
integer,dimension(nx,ny)::level
integer::klo,khi
real::plo,phi,tlo,thi,zlo,zhi,p_at_pconst,t_at_pconst,z_at_pconst,t_surf,t_sea_level
logical::done

do j=1,ny
    do i=1,nx
        done=.false.
        level(i,j)=-1
        do k=1,nz
            if(.not.done.and.p(i,j,k) .lt. p(i,j,1) - PCONST) then
                level(i,j)=k
                mask(i,j)=.false.
                done=.true.
            endif
        enddo
    enddo
enddo

do j=1,ny
    do i=1,nx
        if(.not.mask(i,j))then
            klo=max(level(i,j)-1,1)
            khi=min(klo+1,nz-1)
            plo=p(i,j,klo)
            phi=p(i,j,khi)
            tlo=t(i,j,klo) * (1. + 0.608 * q(i,j,klo))
            thi=t(i,j,khi) * (1. + 0.608 * q(i,j,khi))
    
            zlo=z(i,j,klo)
            zhi=z(i,j,khi)
    
            p_at_pconst=p(i,j,1) - PCONST
            t_at_pconst=thi-(thi-tlo)*log(p_at_pconst/phi)*log(plo/phi)
            z_at_pconst=zhi-(zhi-zlo)*log(p_at_pconst/phi)*log(plo/phi)
    
            t_surf=t_at_pconst*(p(i,j,1)/p_at_pconst)**CEXP
            t_sea_level=t_at_pconst+GAMMA*z_at_pconst
    
            sp(i,j)=p(i,j,1)* &
                exp(2.*G*z(i,j,1)/(R*(t_sea_level+t_surf)))
        endif
    enddo
enddo

end subroutine compute_seaprs_from_derived

end module
