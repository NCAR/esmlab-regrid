C NCLFORTSTART
      subroutine fill_POP_core(nlon,nlat,var,fillmask,msv,tol,ltripole)
      implicit none
      integer nlon
      integer nlat
      real var(nlon,nlat)
      logical fillmask(nlon,nlat)
      real msv
      real tol
      logical ltripole
C NCLEND
Cf2py integer depend(var) :: nlon = shape(var,0), nlat = shape(var,1)
Cf2py intent(inout) var
Cf2py intent(in) fillmask
Cf2py intent(in) msv
Cf2py intent(in) tol
Cf2py intent(in) ltripole

      integer iter
      real work(nlon,nlat)
      real numer, denom, delta
      logical done
      integer i, ip1, im1
      integer j, jp1, jm1

      done = .false.
      iter = 0

      do while (.not. done)
        done = .true.
        iter = iter + 1

C assume bottom row is land, so skip it

        do j = 2, nlat

          jm1 = j-1
          jp1 = j+1

          do i = 1, nlon

            im1 = i-1
            if (i == 1) im1 = nlon
            ip1 = i+1
            if (i == nlon) ip1 = 1

            work(i,j) = var(i,j)

            if (fillmask(i,j)) then
              numer = 0.0
              denom = 0.0

              ! East
              if (var(ip1,j) .ne. msv) then
                numer = numer + var(ip1,j)
                denom = denom + 1.0
              end if

              ! North
              if (j .lt. nlat) then
                if (var(i,jp1) .ne. msv) then
                  numer = numer + var(i,jp1)
                  denom = denom + 1.0
                end if
              else
                ! assume only tripole has non-land top row
                if (ltripole) then
                  if (var(nlon-(i-1),j) .ne. msv) then
                    numer = numer + var(nlon-(i-1),j)
                    denom = denom + 1.0
                  end if
                end if
              end if

              ! West
              if (var(im1,j) .ne. msv) then
                numer = numer + var(im1,j)
                denom = denom + 1.0
              end if

              ! South
              if (var(i,jm1) .ne. msv) then
                numer = numer + var(i,jm1)
                denom = denom + 1.0
              end if

              ! self
              if (var(i,j) .ne. msv) then
                numer = numer + denom*var(i,j)
                denom = 2 * denom
              end if

              if (denom .gt. 0.0) then
                work(i,j) = numer / denom
                if (var(i,j) .eq. msv) then
                  done = .false.
                else
                  delta = abs(var(i,j)-work(i,j))
                  if (delta .gt. tol*abs(var(i,j))) then
                    done = .false.
                  end if
                end if
              end if
            end if

          end do

        end do

        var(:,2:nlat) = work(:,2:nlat)

      end do

      return
      end
