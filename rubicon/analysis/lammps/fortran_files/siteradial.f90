subroutine siteradial(x, y, z, aindex1, aindex2, mol, Lx, Ly, Lz,binsize, numbins, g1, maxr)

integer, dimension(:), intent(in)	:: aindex1
integer, dimension(:), intent(in)	:: aindex2
real, dimension(:), intent(in) :: x
real, dimension(:), intent(in) :: y
real, dimension(:), intent(in) :: z
integer, dimension(:), intent(in)	:: mol
real, intent(in)	:: Lx
real, intent(in)	:: Ly
real, intent(in)	:: Lz
real, intent(in)	:: binsize
integer, intent(in)	:: numbins
real, intent(in)	:: maxr
integer, dimension(0:numbins-1), intent(out) :: g1

real	:: dx,dy,dz
real	:: r
integer	:: ind,ind1,ind2,sa1,sa2
sa1 = size(aindex1)
sa2 = size(aindex2)
if (aindex1(1)==aindex2(1)) then	
	do i=1,sa1-1
		do j=i+1,sa2
			ind1 = aindex1(i)
			ind2 = aindex2(j)
			if (mol(ind1)/=mol(ind2)) then
				dx=x(ind1)-x(ind2)
				dy=y(ind1)-y(ind2)
				dz=z(ind1)-z(ind2)
				dx=dx-Lx*float(nint(dx/Lx))
				dy=dy-Ly*float(nint(dy/Ly))
				dz=dz-Lz*float(nint(dz/Lz))
				r=sqrt(dx**2+dy**2+dz**2)
				if (r<=maxr+3*binsize/2) then
					ind = nint(r/binsize)-1
					g1(ind)=g1(ind)+2
				end if
			end if
		end do
	end do
else
	do i=1,sa1
		do j=1,sa2
			ind1 = aindex1(i)
			ind2 = aindex2(j)
			if (mol(ind1)/=mol(ind2)) then
				dx=x(ind1)-x(ind2)
				dy=y(ind1)-y(ind2)
				dz=z(ind1)-z(ind2)
				dx=dx-Lx*float(nint(dx/Lx))
				dy=dy-Ly*float(nint(dy/Ly))
				dz=dz-Lz*float(nint(dz/Lz))
				r=sqrt(dx**2+dy**2+dz**2)
				if (r<=maxr+3*binsize/2) then
					ind = nint(r/binsize)-1
					g1(ind)=g1(ind)+1
				end if
			end if
		end do
	end do
end if
return 
end subroutine siteradial
