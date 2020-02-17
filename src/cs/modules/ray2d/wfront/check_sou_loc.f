c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

      subroutine check_sou_loc(x_int,z_int,d_int,c_int,b_int,npoints_int,
     &     x_source, z_source, dzdx_source,
     &     int_source, lay_source, error, f_out,
     &     n_int, MAXPOINTS_INT)

      implicit none

      integer n_int, MAXPOINTS_INT, MAX_CODESTEPS
      integer npoints_int(n_int)
      real x_int(MAXPOINTS_INT,n_int), z_int(MAXPOINTS_INT,n_int)
      real b_int(MAXPOINTS_INT,n_int), c_int(MAXPOINTS_INT,n_int)
      real d_int(MAXPOINTS_INT,n_int)
      real x_source, z_source, dzdx_source
      integer int_source, lay_source, error, f_out

      integer i, j, actp
      real ztmp, epsilon

      epsilon = 0.001
c-------------------------------
c Fall: source innerhalb Schicht
c
      if (int_source .eq. 0) then
c 1. Pruefe, ob Quelle wirklich innerhalb der Schicht liegt
         call get_z(x_source,x_int(1,lay_source),z_int(1,lay_source),d_int(1,lay_source),c_int(1,lay_source),
     &        b_int(1,lay_source),npoints_int(lay_source),
     &        ztmp,actp,
     &        MAXPOINTS_INT)
         if (z_source .lt. ztmp) then
            goto 910
         else if (abs(ztmp-z_source) .lt. epsilon) then
            goto 920
         end if
         call get_z(x_source,x_int(1,lay_source+1),z_int(1,lay_source+1),d_int(1,lay_source+1),c_int(1,lay_source+1),
     &        b_int(1,lay_source+1),npoints_int(lay_source+1),
     &        ztmp,actp,
     &        MAXPOINTS_INT)
         if (z_source .gt. ztmp) then
            goto 910
         else if (abs(ztmp-z_source) .lt. epsilon) then
            goto 920
         end if

         dzdx_source = 0.0
c------------------------------
c Fall: source auf Grenzflaeche
      else
         actp = 1
         call get_z1(x_source,x_int(1,int_source),z_int(1,int_source),d_int(1,int_source),c_int(1,int_source),
     &        b_int(1,int_source),npoints_int(int_source),
     &        ztmp,dzdx_source,actp,
     &        MAXPOINTS_INT)
         if (ztmp .ne. z_source) then
            z_source = ztmp
            write(*,1010) z_source
            write(f_out,1010) z_source
         end if
      end if

 1010 format(" WARNING: changed z_source to ",f10.5," ! Otherwise change INT_SOURCE in input file to '0'")

      goto 999
 910  write(*,*) "CHECK_SOU_LOC: source point not inside the right layer!"
      write(*,*) "e.g. change 'Z_SOURCE' in input file"
      write(*,*) "ztmp & z_source: ", ztmp, z_source
      error = 1
      goto 999
 920  write(*,*) "CHECK_SOU_LOC: source point too close to interface!"
      write(*,*) "e.g. change 'INT_SOURCE' from '0' to the right index of interface,"
      write(*,*) "or change 'Z_SOURCE' in such a way, that the distance between"
      write(*,*) "source point and interface is at least ",epsilon
      error = 1
      goto 999

 999   end

