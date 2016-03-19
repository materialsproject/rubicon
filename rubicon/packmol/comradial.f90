subroutine comradial(mol1, mol2, nummol, moltype, comx, comy, comz, Lx, Ly, Lz,binsize, numbins, g1, maxr)

integer, intent(in)	:: mol1
integer, intent(in)	:: mol2
integer, intent(in)	:: nummol
integer, dimension(:), intent(in) :: moltype
real, dimension(:), intent(in) :: comx
real, dimension(:), intent(in) :: comy
real, dimension(:), intent(in) :: comz
real, intent(in)	:: Lx
real, intent(in)	:: Ly
real, intent(in)	:: Lz
real, intent(in)	:: binsize
integer, intent(in)	:: numbins
real, intent(in)	:: maxr
integer, dimension(0:numbins-1), intent(out)	:: g1
real	:: dx,dy,dz
real	:: r
integer	:: ind	
	if(mol1==mol2) then
		do i=1,nummol-1
			if (moltype(i)==mol1) then
				do j=i+1, nummol
					if (moltype(j)==mol2) then
						dx=comx(i)-comx(j)
						dy=comy(i)-comy(j)
						dz=comz(i)-comz(j)
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
			end if
		end do
	
	else
		do i=1,nummol
			if (moltype(i)==mol1) then
				do j=1,nummol
					if (moltype(j)==mol2) then
						dx=comx(i)-comx(j)
						dy=comy(i)-comy(j)
						dz=comz(i)-comz(j)
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
			end if
		end do 
	end if
return 
end subroutine comradial
