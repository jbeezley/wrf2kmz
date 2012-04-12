
subroutine reprojectionIdx(nx,ny,inx,iny,lon,lat,xi,yi,idx,ierr)
implicit none
integer,intent(in)::nx,ny,inx,iny
real,dimension(nx,ny),intent(in)::lon,lat
real,dimension(inx),intent(in)::xi
real,dimension(iny),intent(in)::yi
integer,dimension(2,inx,iny),intent(out)::idx
integer,intent(out)::ierr

integer::i,j,i0,j0,ic
integer,dimension(4),parameter::ip=(/-1,0,0,-1/),jp=(/-1,-1,0,0/)
logical,dimension(inx,iny)::m
real,dimension(4)::xc,yc
integer::inpolygon
external::inpolygon
do i=2,inx
    if(xi(i).le.xi(i-1)) then
        ierr=2
        return
    endif
enddo

do j=2,iny
    if(yi(j).le.yi(j-1)) then
        ierr=2
        return
    endif
enddo

m=.false.
idx=-1
do j=2,ny
    do i=2,nx
        if(lon(i,j).le.lon(i-1,j))then
            print*,'lon ',i,j,lon(i,j),'<=',lon(i-1,j)
            ierr=3
            return
        endif
        if(lat(i,j).le.lat(i,j-1))then
            print*,'lat ',i,j,lat(i,j),'<=',lat(i,j-1)
            ierr=3
            return
        endif
        do j0=1,iny
            if(yi(j0).lt.lat(i-1,j-1).and.yi(j0).lt.lat(i,j-1)) goto 10
            if(yi(j0).gt.lat(i-1,j).and.yi(j0).gt.lat(i,j)) goto 20
            do i0=1,inx
                if(m(i0,j0)) goto 30
                if(xi(i0).lt.lon(i-1,j-1).and.xi(i0).lt.lon(i-1,j)) goto 30
                if(xi(i0).gt.lon(i,j-1).and.xi(i0).gt.lon(i,j)) goto 20
                
                xc=(/lon(i-1,j-1),lon(i,j-1),lon(i,j),lon(i-1,j)/)
                yc=(/lat(i-1,j-1),lat(i,j-1),lat(i,j),lat(i-1,j)/)
                ic=inpolygon( xi(i0),yi(j0),xc,yc )
                if(ic.gt.0) then 
                    idx(1,i0,j0)=i+ip(ic)
                    idx(2,i0,j0)=j+jp(ic)
                    m(i0,j0)=.true.
                endif

                30 continue
            enddo
            10 continue
        enddo
        20 continue
    enddo
enddo

end subroutine reprojectionIdx

subroutine interpArray(nx,ny,inx,iny,idx,a,fill,b)
implicit none
integer,intent(in)::nx,ny,inx,iny
real,dimension(nx,ny),intent(in)::a
real,intent(in)::fill
integer,dimension(2,inx,iny),intent(in)::idx
real,dimension(inx,iny),intent(out)::b

integer::i,j
do j=1,iny
    do i=1,inx
        if(idx(1,i,j).gt.0)then
            b(i,j)=a(idx(1,i,j),idx(2,i,j))
        else
            b(i,j)=fill
        endif
    enddo
enddo

end subroutine interpArray

subroutine reprojectArray(nx,ny,inx,iny,lon,lat,a,xi,yi,fill,b,ierr)
implicit none
integer,intent(in)::nx,ny,inx,iny
real,dimension(nx,ny),intent(in)::lon,lat,a
real,dimension(inx),intent(in)::xi
real,dimension(iny),intent(in)::yi
real,intent(in)::fill
real,dimension(inx,iny),intent(out)::b
integer,intent(out)::ierr

integer::i,j,i0,j0,ic
integer,dimension(4),parameter::ip=(/-1,0,0,-1/),jp=(/-1,-1,0,0/)
logical,dimension(inx,iny)::m
real,dimension(4)::xc,yc
integer::inpolygon
external::inpolygon

do i=2,inx
    if(xi(i).le.xi(i-1)) then
        ierr=2
        return
    endif
enddo

do j=2,iny
    if(yi(j).le.yi(j-1)) then
        ierr=2
        return
    endif
enddo

b=fill
m=.false.
do j=2,ny
    do i=2,nx
        if(lon(i,j).le.lon(i-1,j))then
            print*,'lon ',i,j,lon(i,j),'<=',lon(i-1,j)
            ierr=3
            return
        endif
        if(lat(i,j).le.lat(i,j-1))then
            print*,'lat ',i,j,lat(i,j),'<=',lat(i,j-1)
            ierr=3
            return
        endif
        do j0=1,iny
            if(yi(j0).lt.lat(i-1,j-1).and.yi(j0).lt.lat(i,j-1)) goto 10
            if(yi(j0).gt.lat(i-1,j).and.yi(j0).gt.lat(i,j)) goto 20
            do i0=1,inx
                if(m(i0,j0))goto 20
                if(xi(i0).lt.lon(i-1,j-1).and.xi(i0).lt.lon(i-1,j)) goto 30
                if(xi(i0).gt.lon(i,j-1).and.xi(i0).gt.lon(i,j)) goto 20
                
                xc=(/lon(i-1,j-1),lon(i,j-1),lon(i,j),lon(i-1,j)/)
                yc=(/lat(i-1,j-1),lat(i,j-1),lat(i,j),lat(i-1,j)/)
                ic=inpolygon( xi(i0),yi(j0),xc,yc )
                !print*,i,j,i0,j0,ic
                if(ic.gt.0) then 
                    b(i0,j0)=a(i+ip(ic),j+jp(ic))
                    m(i0,j0)=.true.
                endif

                30 continue
            enddo
            10 continue
        enddo
        20 continue
    enddo
enddo

end subroutine reprojectArray

integer function inpolygon(x,y,xc,yc)
implicit none

integer,parameter::n=4
real,intent(in)::x,y
real,dimension(n),intent(in)::xc,yc

integer::i
real::a,b,c,d,m
logical::sgn
real::norm
external::sgn,norm

a=(xc(1)-x) * (yc(2)-y) - (xc(2)-x) * (yc(1)-y)
b=(xc(2)-x) * (yc(3)-y) - (xc(3)-x) * (yc(2)-y)
c=(xc(3)-x) * (yc(4)-y) - (xc(4)-x) * (yc(3)-y)
d=(xc(4)-x) * (yc(1)-y) - (xc(1)-x) * (yc(4)-y)

if(sgn(a,b).and.sgn(b,c).and.sgn(c,d))then
    a=norm(xc(1),yc(1),x,y)
    b=norm(xc(2),yc(2),x,y)
    c=norm(xc(3),yc(3),x,y)
    d=norm(xc(4),yc(4),x,y)

    i=1
    m=a
    if(b.lt.m)then
        i=2
        m=b
    endif
    if(c.lt.m)then
        i=3
        m=c
    endif
    if(d.lt.m)then
        i=4
        m=d
    endif
    inpolygon=i
else
    inpolygon=0
endif

end function inpolygon

pure logical function sgn(a,b)
implicit none
real,intent(in)::a,b
if((a.ge.0.and.b.ge.0).or.(a.lt.0.and.b.lt.0))then
    sgn=.true.
else
    sgn=.false.
endif
end function

pure real function norm(x1,y1,x2,y2)
real,intent(in)::x1,y1,x2,y2
real::x0,y0
x0=x1-x2
y0=y1-y2
norm=x0*x0+y0*y0
end function norm

