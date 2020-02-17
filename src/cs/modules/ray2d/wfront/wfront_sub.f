c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

      subroutine wfsub(n_int,npoints_int,x_int,z_int,iii,rho1,rho2,ptos,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     b_int,c_int,d_int,
     &     left_border,right_border,iii2,
     &     nray_source, lay_source, int_source, x_source, z_source, t_source, dzdx_source, anglein, angleout,
     &     dtray,spread_dist,angle_max,flag_smooth,flag_compounds,flag_surface,flag_stack,flag_2point,
     &     nreclines,int_rec,lay_rec,nrec,xz_rec,nboreholes,
     &     step_spread,step_wfronts,step1_wfronts,last_wfronts,
     &     code, codesteps, ncodes, codestep_plot, hold,
     &     f_out,f_wfronts, f_xrays,
     &     f_rayout, f_rayouttmp, f_timeout,
     &     MAX_RECEIVERS, MAX_CODESTEPS,
     &     MAXP_XGRID, MAXP_ZGRID, MAX_2POINT_ANGLES,
     &     n_allint, MAXPOINTS_INT, MAX_ARRIVALS,
     &     time, amp_out, phase_out, angle_out, comp_out, ntimecodes, ccp )

      implicit none

c----------------------
c Parameter
c
      integer n_int, MAXPOINTS_INT, n_allint, MAX_RAYS, N_PARAM, N_RTPARAM
      integer MAXP_XGRID, MAXP_ZGRID, MAX_HALTS, MAX_CODESTEPS
      integer MAX_RECEIVERS, MAX_ARRIVALS, nreclines, N_ARRPARAM
      integer MAX_TIMEFILES

      parameter(MAX_RAYS = 40000, N_PARAM = 8, N_RTPARAM = 8)
      parameter(MAX_HALTS = 20)
      parameter(N_ARRPARAM = 6)
      parameter(MAX_TIMEFILES = 80)
c----------------------
c Das Modell, Interfaces, Geschwindigkeiten, Dichten...
c
      integer npoints_int(n_allint), iii(MAXPOINTS_INT,n_allint)
      integer iii2(MAXPOINTS_INT,n_int)
      real x_int(MAXPOINTS_INT,n_allint), z_int(MAXPOINTS_INT,n_allint)
      real b_int(MAXPOINTS_INT,n_allint), c_int(MAXPOINTS_INT,n_allint)
      real d_int(MAXPOINTS_INT,n_allint)
      integer nx_grid(n_int),nz_grid(n_int)
      real x_grid(MAXP_XGRID,n_int),z_grid(MAXP_ZGRID,n_int)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID,n_int)
      real rho1(n_int), rho2(n_int), ptos(n_int)
      real left_border, right_border
c----------------------
c Die Quelle  (isource ist neu fuer WFSTACK, zaehlt die berechneten Schuesse mit)
c
      integer lay_source, int_source, nray_source, isource, veltype_source
      real x_source, z_source, t_source, angle_source, dzdx_source, vp_source, vs_source, rho_source
      real v_source, d_angle, anglein, angleout
c----------------------
c Feste Werte, Genauigkeitswerte
c
      real pi, pihalf
      real spread_dist, spread_max, cos_min, angle_max, spread_max2, inf_dist, inf_angle
c Bedeutung: 
c FLAG_SMOOTH = 0: normal; != 0: Bei eckigen Punkten wird Wellenfront nicht getrennt
c FLAG_COMPOUNDS = 0: normal; != 0: keine compound elements
c FLAG_SURFACE = 0: normal; != 0: keine Reflexionen an Oberflaeche
      integer flag_smooth,flag_compounds,flag_surface,flag_stack
c----------------------
c Aktuelle und loslaufende Strahlen
c
      real ray(N_PARAM,MAX_RAYS),rayold(N_PARAM)
      real chkbound(MAX_RAYS), amplitude(MAX_RAYS), phase(MAX_RAYS), angle(MAX_RAYS)
      integer nextray(0:MAX_RAYS), kmah(MAX_RAYS)
      real ray0(N_PARAM,MAX_RAYS), tray0(MAX_RAYS)
      real reflpoint(2,MAX_RAYS)
      real boundtmp(N_PARAM),ttmp
      integer nextray0(0:MAX_RAYS), kmah0(MAX_RAYS)
      integer nray, nray0, nrayold, iray, irayold, iray0, iray0old, iraytmp, newcross
      integer halfsteps_count
c----------------------
c Momentane, aktuelle Werte
c
      integer actp,actint, updownflag, actlayer, oldlayer
      real zint, chkbound0,dzdx
      real acttime, dtray, dthalf, dtnew, dtnewhalf, vel(6)
      integer actncomp, actcompounds(MAX_CODESTEPS)
      integer veltype, oldveltype
      integer timesteps
c----------------------
c Intersection und Suche nach boundary ray
c
      real thesect(N_PARAM,MAX_RAYS),tsect(MAX_RAYS),rtsect(N_RTPARAM,MAX_RAYS)
      integer nextsect(0:MAX_RAYS),crossect(MAX_RAYS),kmahsect(MAX_RAYS), fictsect(MAX_RAYS)
      real raysect(4)
c----------------------
c Aufsammlung bei Interfaces
c
c TEMP
c      integer N_PARAM_HALT
c      parameter( N_PARAM_HALT=7 );
c END TEMP
      real thehalt(N_PARAM,MAX_RAYS,MAX_HALTS), thalt(MAX_RAYS,MAX_HALTS)
      real amplhalt(MAX_RAYS,MAX_HALTS), rthalt(N_RTPARAM,MAX_RAYS,MAX_HALTS)
      real phasehalt(MAX_RAYS,MAX_HALTS),anglehalt(MAX_RAYS,MAX_HALTS)
      integer nexthalt(0:MAX_RAYS,MAX_HALTS),kmahhalt(MAX_RAYS,MAX_HALTS),ficthalt(MAX_RAYS,MAX_HALTS)
      integer ncodehalt(MAX_HALTS),codestepshalt(MAX_HALTS)
      integer ncomphalt(MAX_HALTS),compoundshalt(MAX_CODESTEPS,MAX_HALTS)
      integer inthalt(MAX_HALTS), nrayhalt(MAX_HALTS)
      integer acthalt, actualhalt(MAX_CODESTEPS), ihaltup, ihaltdown, ihaltnorm, ihaltcomp
c----------------------
c Boreholes
c
      integer nboreholes, ibore
c----------------------
c Boreholes halts
c
      real leftboresect(N_PARAM,nboreholes,MAX_RAYS),tleftboresect(nboreholes,MAX_RAYS)
      integer kmahlboresect(nboreholes,MAX_RAYS)
      real rightboresect(N_PARAM,nboreholes,MAX_RAYS),trightboresect(nboreholes,MAX_RAYS)
      integer kmahrboresect(nboreholes,MAX_RAYS)
      integer nextboresect(0:MAX_RAYS)

      real borehalt(N_PARAM,MAX_RAYS), tborehalt(MAX_RAYS)
      real amplborehalt(MAX_RAYS), phaseborehalt(MAX_RAYS), angleborehalt(MAX_RAYS)
      integer nextborehalt(0:MAX_RAYS), kmahborehalt(MAX_RAYS)
      integer ncodeborehalt,codestborehalt
      integer ncompborehalt,compborehalt(MAX_CODESTEPS)
      integer nrayborehalt
c----------------------
c Codes...
c
      integer code(0:MAX_CODESTEPS,ncodes),codesteps(ncodes),hold(0:MAX_CODESTEPS,ncodes)
      integer ncodes
      integer actcode, oldcode, actcodestep(MAX_CODESTEPS), nactcodestep
c----------------------
c Reflection/Transmission
c
      integer signrt
c----------------------
c Receiver, Einsaetze
c
      real arrival(N_ARRPARAM,MAX_RECEIVERS,MAX_ARRIVALS,nreclines)
      real xz_rec(MAX_RECEIVERS,nreclines)
      integer int_rec(nreclines), lay_rec(nreclines), nrec(nreclines)
      integer narrivals(MAX_RECEIVERS,nreclines), xzindex, all_reclines
      integer iline
c----------------------
c Ausgabe...
c
      integer f_out, f_wfronts, f_xrays
      integer step_spread, step_wfronts, step1_wfronts, last_wfronts
      integer codestep_plot(ncodes)
      integer icompout(nreclines)
      integer f_rayout, f_rayouttmp
      character*30 fname_out,fname_wfronts,fname_rayout,fname_timeout
      character*30 fname_rayouttmp
c----------------------
c Laufvariablen
c
      integer n,i,j,k
c----------------------
c sonstige...
c
      integer idummy1,idummy2,idummy3,idummy4, lay_bore, int_bore, oldcross
      integer abbruch, ficttmp
      real ampltmp, phasetmp, arr2, arr3, arr4, arr5, x, z, ay, phy, tast
      integer flag, error, next
c
c----------------------
c New variables for kind of two-point ray-tracing
c Crude way of doing this...

c Allow first 3 arrivals
      integer MAX_2POINT_ANGLES
      integer flag_2point
      logical process_2point
      real angles_2point(MAX_2POINT_ANGLES,ncodes)
      integer nangles_2point(ncodes)
c
c TEMP
      integer f_timeout
c Arrival times for first receiver line:
c MAX_ARRIVALS: Maximum number of arrivals for each ray code
      integer ntimecodes, act_timecode
      real time(MAX_RECEIVERS,MAX_ARRIVALS,ntimecodes)
      real amp_out(MAX_RECEIVERS,MAX_ARRIVALS,ntimecodes)
      real phase_out(MAX_RECEIVERS,MAX_ARRIVALS,ntimecodes)
      real angle_out(MAX_RECEIVERS,MAX_ARRIVALS,ntimecodes)
      real ccp(2,MAX_RECEIVERS,ntimecodes)
      integer comp_out

c Ende Deklarationen
c================================================================================================
c================================================================================================

      abbruch = 0
      process_2point = .false.
c
      pi = 3.1415927
      pihalf = 1.5707963
      inf_angle = pi / 2000.
      inf_dist = 0.0001
c inf_dist = 0.0001 ausreichend, weniger bringt nichts!

      halfsteps_count = 0
      isource = 1

c Werte, die immer 0.0 bleiben (werden fuer Seismogramme gebraucht)
      ay = 0.0
      phy = 0.0
      tast = 0.0

c-------------------------------------------------
c Initialisierung
      do i = 1, MAX_HALTS
         codestepshalt(i) = 0
      enddo
c-------------------------------------------------
c Einlesen aller Eingabeparameter --> bereits erledigt
c
      error = 0

      if (error .ne. 0) then
         if (error .gt. 0) goto 990
         goto 995
      end if
c      angle_max = angle_max*pi/180.
      cos_min = cos(angle_max)
      spread_max  = spread_dist * spread_dist
      spread_max2 = spread_max*16.
c-------------------------------------------------
c Erstellen des Modells, Geschwindigkeitsmodell --> bereits erledigt

      all_reclines = nreclines + nboreholes
c-------------------------------------------------
c Sortieren der Codes --> bereits erledigt
c

c--------------------------------------------------
c Initialisierung, Setzen aller notwendigen Werte - bereits erledigt
c 
c Erst Werte, die nur einmal gesetzt werden muessen...
c
      dthalf = dtray / 2.0

      do i = 1, all_reclines
         do j = 1, MAX_RECEIVERS
            narrivals(j,i) = 0
         end do
      end do

c TEMP
c Try to avoid too many rays...
      do iray = 1, MAX_RAYS
         nextray0(iray) = 0
         nextray(iray) = 0
         nextsect(iray) = 0
      end do
c END TEMP
c
c ...dann Werte, die beim Stack neu gesetzt werden muessen
c
 111  nray = nray_source

      if( process_2point ) then
c         write(*,*)
c         write(*,*) ".....TWO-POINT RAYTRACING, FINAL RUN:"
         nray = nangles_2point(actcode)
c      else if( flag_2point .ne. 0 ) then
c         write(*,*)
c         write(*,*) ".....TWO-POINT RAYTRACING, INITIAL RUN:"
      endif

      timesteps = 0
      actcode = 1
      act_timecode = 1
      nactcodestep = 1
      actcodestep(1) = 1
      acthalt = 0
      actncomp = 0
      actlayer = abs(code(1,1))

c TIME.OUT, Ausgabedatei der Laufzeiten; ICOMPOUT: =1,2,3... fuer verschiedene Laufwege von ein und demselben
c Strahlencode.
      do iline = 1, all_reclines
         icompout(iline) = 1
      end do

c-----------------------------------------------------
c Indizes der naechsten Halts setzen
      ihaltnorm = 1
      if ((codesteps(1) .gt. 1 .and. code(2,1) .eq. code(1,1)) .or. hold(0,1) .ne. 0 .or.
     &     hold(1,1) .gt. 1) then
         if (flag_compounds .eq. 0) then
            ihaltcomp = 2
         else
            ihaltcomp = 0
         end if
      else
         ihaltcomp = 0
      end if
      if (actlayer .eq. int_source .or. code(0,1) .gt. 0) then
         ihaltdown = ihaltnorm
         ihaltup = ihaltcomp
      else
         ihaltup = ihaltnorm
         ihaltdown = ihaltcomp
      end if

      actualhalt(1) = ihaltnorm

c=======================================================================================
c Quelle, Winkel...
c
c------------------------------------------------------------------
c Quelle liegt innerhalb Schicht:
      if (int_source .eq. 0) then
         if (anglein .ne. 0.0 .or. angleout .ne. 0.0) then
            angle_source = anglein
            d_angle = (angleout - anglein) / float(nray-1)
         else
            angle_source = 0.0
            d_angle = 2.*pi / float(nray-1)
         end if
c------------------------------------------------------------------
c Quelle liegt auf Interface:
      else
         angle_source = sign(1.0,dzdx_source) * acos(1.0/sqrt(1.0 + dzdx_source*dzdx_source))
c----------------------
c kein beam:
         if (anglein .eq. 0.0 .and. angleout .eq. 0.0) then
            d_angle = (pi - 2.0*inf_angle) / float(nray-1)
            angle_source = angle_source + inf_angle
c Fall: Strahl laeuft nach oben:
            if (actlayer .lt. int_source) then
               angle_source = angle_source + pi
            end if
c----------------------
c beam, Einschraenkung durch ANGLEIN, ANGLEOUT
         else
c Fall: Strahl laeuft nach oben:
            if (actlayer .lt. int_source) then
               angle_source = angle_source + pi
            end if
            if (anglein .lt. angle_source) then
               if (angleout .gt. angle_source) then
                  goto 800
               end if
c Falls beam in die andere Schicht zeigt, hier kein Strahl!:
               nray = 0
            else if (anglein .lt. angle_source+pi) then
               if (angleout .gt. angle_source+pi) then
                  goto 800
               end if
               if (angleout .gt. angle_source+pi-inf_angle) then
                  angleout = angle_source+pi-inf_angle
                  write(f_out,*,err=915) "Changed angleout to ",angleout*180./pi
               end if
               if (anglein .lt. angle_source+inf_angle) then
                  anglein = angle_source + inf_angle
                  write(f_out,*,err=915) "Changed anglein to ",anglein*180./pi
               end if
               angle_source = anglein
               d_angle = (angleout - angle_source)/float(nray-1)
            else if (angleout .gt. angle_source+2.*pi) then
               goto 800
c Falls beam in die andere Schicht zeigt, hier kein Strahl!:
            else
               nray = 0
            end if
         end if
      end if
c=========================================================================================
c Geschwindigkeiten, Dichte setzen
c
      veltype = isign(1,code(1,1))

      call velocity(x_source, z_source, nx_grid(lay_source), nz_grid(lay_source),
     &     x_grid(1,lay_source), z_grid(1,lay_source), v_grid(1,1,1,lay_source),
     &     veltype, ptos(lay_source),lay_source,
     &     vel,
     &     MAXP_XGRID, MAXP_ZGRID)

      if (error .ne. 0) goto 990

      veltype_source = veltype
      v_source = vel(1)
      if (veltype_source .gt. 0) then
         vp_source = v_source
         vs_source = v_source/ptos(lay_source)
      else
         vs_source = v_source
         vp_source = v_source*ptos(lay_source)
      end if
      rho_source = rho1(lay_source) + rho2(lay_source)*vp_source

c-----------------------
c Ausgabe in F_SEISMOUT.
c
c ...geloescht...
c
c-----------------------
c Belegen der Initialisierungsstrahlen ray0(..)
c 
      do iray = 1, nray
         tray0(iray) = t_source
         ray0(1,iray) = x_source
         ray0(2,iray) = z_source

         nextray0(iray) = iray + 1
         nextray(iray) = 0
         nextsect(iray) = nextray0(iray)

         if( process_2point ) then
            angle(iray) = angles_2point(iray,actcode)
c Make next ray "unconnected" to this ray, to avoid interpolation of new rays:
            nextsect(iray) = -nextray0(iray)
         else
            angle(iray) = angle_source + d_angle*float(iray-1)
         endif

         ray0(3,iray) = cos(angle(iray))/vel(1)
         ray0(4,iray) = sin(angle(iray))/vel(1)
         ray0(5,iray) = 0.0
         ray0(6,iray) = 1.0
         kmah0(iray) = 0
         ray0(7,iray) = x_source
         ray0(8,iray) = z_source

         crossect(iray) = 0
         kmahsect(iray) = 0
         fictsect(iray) = 0
         amplitude(iray) = vel(1)*sqrt(rho_source)
         phase(iray) = 0.0
         nextboresect(iray) = nextray0(iray)
         do ibore = 1, nboreholes
            tleftboresect(ibore,iray) = 0.0
            trightboresect(ibore,iray) = 0.0
         end do
      end do
c---------
c Quelle liegt auf conciding interface:
c Noch ist nur Erdoberflaeche als Interface moeglich...
c
      if (int_source .eq. 1) then
         call get_actp(x_source,int_source,x_int,npoints_int,
     &                actp,
     &                n_int,MAXPOINTS_INT)
         if (iii2(actp,int_source) .ne. 0) then
            ficttmp = iii2(actp,int_source)*1000 - 1000
            do iray = 1, nray
               fictsect(iray) = ficttmp
            end do
         end if
      end if

c Zur Sicherheit, damit der erste und letzte Strahl in diesem Falle wirklich gleich sind:
      if (int_source .eq. 0 .and. anglein .eq. 0.0 .and. angleout .eq. 0.0 .and. .not. process_2point ) then
         do i = 3, 4
            ray0(i,nray) = ray0(i,1)
         end do
         angle(nray) = angle(1)
      end if

      nextray0(0) = 1
      nextray0(nray) = 0
      nextray(0) = 0
      nextsect(0) = -1
      nextsect(nray) = 0

      nextboresect(0) = -1
      nextboresect(nray) = 0

      nray0 = nray

      if (f_rayout .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
         if (actcodestep(nactcodestep) .ge. codestep_plot(actcode) .and. timesteps .le. last_wfronts) then
            do iray = 1, nray
               write(f_rayouttmp,err=960) iray,ray0(1,iray),ray0(2,iray)
            end do
         end if
      end if

c
c Initialisierung Ende
c*****************************************************************************
c*****************************************************************************
c Hier geht's los
c-------------------------------------------------------------------------
c Aeussere Schleife:
c Laeuft, bis alle Strahlen terminieren.
c 
 100  continue
      do while (nextray(0) .ne. 0 .or. nextray0(0) .ne. 0)
         iray = 0
         irayold = 0
         timesteps = timesteps + 1
         acttime = float(timesteps)*dtray
c==========================================================
c Innere Schleife:
c Alle aktuellen Strahlen werden um DTRAY weiterpropagiert.
         flag = 0
         next = abs(nextray(iray))
         do while(next .ne. 0)
            iray = next
            next = abs(nextray(iray))
            do i = 1, N_PARAM
               rayold(i) = ray(i,iray)
            end do
            oldcross = crossect(iray)

            call propagate(rayold,actlayer,acttime,dtray,dthalf,ihaltup,ihaltdown,veltype,flag,
     &           n_int,x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos(actlayer),left_border,right_border,
     &           nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),v_grid(1,1,1,actlayer),
     &           ray(1,iray),kmah(iray),chkbound(iray),irayold,iray,nextray(irayold),nextray(iray),
     &           thesect(1,iray),tsect(iray),rtsect(1,iray),kmahsect(iray),crossect(iray),fictsect(iray),
     &           leftboresect(1,1,iray),tleftboresect(1,iray),kmahlboresect(1,iray),
     &           rightboresect(1,1,iray),trightboresect(1,iray),kmahrboresect(1,iray),halfsteps_count,
     &           n_allint,MAXPOINTS_INT,N_PARAM,N_RTPARAM,nboreholes,MAXP_XGRID,MAXP_ZGRID)

            if (error .ne. 0) goto 990
            if (f_rayouttmp .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
               if (actcodestep(nactcodestep) .ge. codestep_plot(actcode) .and. timesteps.le.last_wfronts) then
                  if (crossect(iray) .eq. 0) then
                     write(f_rayouttmp,err=960) iray,ray(1,iray),ray(2,iray)
                  else if (oldcross .eq. 0) then
                     write(f_rayouttmp,err=960) iray,thesect(1,iray),thesect(2,iray)
c temporary
c                  else
c                     write(f_rayouttmp,err=960) iray,ray(1,iray),ray(2,iray)
                  end if
               end if
            end if

         end do
c Ende innere Schleife
c==========================================================


c==========================================================
c Dieselbe innere Schleife zur Herstellung der neuen Wellenfront
c
         flag = -1
         iray0 = 0
         iray0old = 0
         do while(nextray0(iray0) .ne. 0)
            iray0 = nextray0(iray0)
c Sicherheitsepsilon fuer den Fall, dass TRAY0  g e n a u  gleich ACTTIME ist
c (dann wuerde es naemlich einen unerwarteten Schnittpunkt geben!)
            if ((acttime - tray0(iray0)) .gt. 0.000001) then
c Nur Propagation, wenn man nicht innerhalb von coinciding interfaces ist
               if (fictsect(iray0) .eq. 0) then
                  dtnew = acttime - tray0(iray0)
                  dtnewhalf = 0.5*dtnew
                  idummy2 = 1
                  call propagate(ray0(1,iray0),actlayer,acttime,dtnew,dtnewhalf,ihaltup,ihaltdown,veltype,flag,
     &                 n_int,x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos(actlayer),left_border,right_border,
     &                 nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),v_grid(1,1,1,actlayer),
     &                 ray(1,iray0),kmah0(iray0),chkbound(iray0),idummy3,idummy1,idummy4,idummy2,
     &                 thesect(1,iray0),tsect(iray0),rtsect(1,iray0),kmahsect(iray0),crossect(iray0),fictsect(iray0),
     &                 leftboresect(1,1,iray0),tleftboresect(1,iray0),kmahlboresect(1,iray0),
     &                 rightboresect(1,1,iray0),trightboresect(1,iray0),kmahrboresect(1,iray0),halfsteps_count,
     &                 n_allint, MAXPOINTS_INT, N_PARAM, N_RTPARAM, nboreholes, MAXP_XGRID, MAXP_ZGRID )

                  if (error .ne. 0) goto 990
                  kmah(iray0) = kmah0(iray0)
                  if (f_rayouttmp .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
                     if (actcodestep(nactcodestep) .ge. codestep_plot(actcode).and.timesteps.le.last_wfronts) then
                        if (crossect(iray0) .eq. 0) then
                           write(f_rayouttmp,err=960) iray0,ray(1,iray0),ray(2,iray0)
                        else
                           write(f_rayouttmp,err=960) iray0,thesect(1,iray0),thesect(2,iray0)
                        end if
                     end if
                  end if
c----------------------------------------------------------
c Bestimmung von nextray():
c
                  if (idummy2 .ne. 0) then
c                  if (ray(1,iray0) .gt. left_border .and. ray(1,iray0) .lt. right_border) then
                     do i = iray0, 0, -1
                        iraytmp = i
                        next = nextray(i)
                        if (next .ne. 0) goto 22
                     end do
 22                  next = abs(next)
                     if (next .lt. iray0) then
                        if ((iray0 - next) .eq. 1) then
                           nextray(next) = nextsect(next)
                        else
                           nextray(next) = -iray0
                        end if
                        nextray(iray0) = 0
                     else if (next .gt. nray0) then
                        i = iraytmp
                        n = abs(nextsect(i))
                        do while((n .lt. iray0 .or. n .gt. nray0) .and. next .ne. 0)
                           i = n
                           n = abs(nextsect(i))
                           if (i .eq. next) then
                              iraytmp = i
                              next = abs(nextray(i))
                           end if
                        end do
                        if ((iray0-iraytmp) .eq. 1) then
                           nextray(iraytmp) = nextsect(iraytmp)
                        else
                           nextray(iraytmp) = -iray0
                        end if
                        if (next .eq. 0) then
                           nextray(iray0) = 0
                        else if ((next-iray0) .eq. 1) then
                           nextray(iray0) = nextsect(iray0)
                        else
                           nextray(iray0) = -next
                        end if
                     else
                        if ((iray0-iraytmp) .eq. 1) then
                           nextray(iraytmp) = nextsect(iraytmp)
                        else
                           nextray(iraytmp) = -iray0
                        end if
                        if ((next-iray0) .eq. 1) then
                           nextray(iray0) = nextsect(iray0)
                        else
                           nextray(iray0) = -next
                        end if
                     end if
                  end if
c--------------------------------
c innerhalb coinciding interfaces:
               else
                  do i = 1, N_PARAM
                     thesect(i,iray0) = ray0(i,iray0)
                  end do
c rtsect wird sowieso nicht gebraucht!
                  tsect(iray0) = tray0(iray0)
                  kmahsect(iray0) = kmah0(iray0)
                  tsect(iray0) = tray0(iray0)
                  if (updownflag .lt. 0) then
                     crossect(iray0) = -actlayer
                  else
                     crossect(iray0) = actlayer+1
                  end if
               end if
               nextray0(iray0old) = nextray0(iray0)
            else
               iray0old = iray0
            end if
         end do

c
c Ende Schleife: Aufbau der neuen Wellenfront
c==========================================================

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c Abfragerei ueber spez. Termination sowie Interpolation
c Dafuer werden alle aktuellen Strahlen nochmals durchlaufen
c Ausserdem wird nach boundary rays gesucht

c----------------------------------------------------------------------------------------
c spezielle Termination:
c Gilt fuer den Fall, dass ein Strahl eine Grenzflaeche ueberschritten hat.
c Haben ausserdem beide Nachbarstrahlen terminiert oder eine Grenzflaeche ueberschritten,
c terminiert der Strahl.
c
         irayold = 0
         iray = abs(nextray(0))
         if (iray .ne. 0) then
            next = abs(nextray(iray))
            do while (next .ne. 0)
               if (crossect(iray) .eq. 0) then
                  irayold = iray
               else if (crossect(next) .ne. 0) then
                  if (nextray(irayold) .lt. 0 .or. crossect(irayold) .ne. 0 .or.
     &                 kmah(irayold) .ne. kmah(iray)) then
                     nextray(irayold) = -next
                     nextray(iray) = 0
                  else
                     irayold = iray
                  end if
               else if (nextray(iray) .lt. 0 .or. kmah(iray) .ne. kmah(next)) then
                  if (nextray(irayold) .lt. 0 .or. crossect(irayold) .ne. 0 .or.
     &                 kmah(irayold) .ne. kmah(iray)) then
                     nextray(irayold) = -next
                     nextray(iray) = 0
                  else
                     irayold = iray
                  end if
               else
                  irayold = iray
               end if
               iray = next
               next = abs(nextray(iray))
            end do
            if (iray .ne. 0 .and. crossect(iray) .ne. 0) then
               if (nextray(irayold) .lt. 0 .or. crossect(irayold) .ne. 0 .or.
     &              kmah(irayold) .ne. kmah(iray)) then
                  nextray(irayold) = 0
               end if
            end if
         end if

c----------------------------------------------------------------------------------------------
c Interpolation und Boundray-Suche:
c (Termination nur evt. bei boundary ray moeglich)
c
         nrayold = nray
         irayold = 0
         iray = 0
         next = abs(nextray(iray))
         do while (next .ne. 0)
            iray = next
            next = abs(nextray(iray))
            if (next .ne. 0 .and. kmah(iray) .ne. kmah(next)) then
               irayold = iray
            else
               if (crossect(iray) .eq. 0) then
                  if (nextray(iray) .gt. 0) then
                     if (crossect(next) .eq. 0) then
c=======================
c Interpolation 1
c CROSSECT(IRAY) = 0; CROSSECT(NEXT) = 0

                        call do_interpol(actlayer,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &                       spread_max,cos_min,angle_max,next,
     &                       ray,amplitude,phase,angle,kmah,
     &                       crossect,fictsect,kmahsect,nextray,nextsect,irayold,iray,nray,
     &                       nextboresect,tleftboresect,trightboresect,error,
     &                       N_PARAM,n_int,MAXPOINTS_INT,nboreholes,MAX_RAYS)

                        if (error .ne. 0) goto 990
c Ende Interpolation 1
c===========================
                     else
c Faelle 5D,5C,9C,9D:
                        actint = abs(crossect(next))
                        call intersect_wfront(ray(1,next),ray(1,iray),actint,crossect(next),
     &                       npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos(actlayer),actlayer,
     &                       nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),
     &                       v_grid(1,1,1,actlayer),
     &                       dzdx,raysect,abbruch,
     &                       n_int,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID)                        
                        chkbound0 = dzdx*raysect(3) - raysect(4)
                        if (abbruch .eq. 0 .and. chkbound0*chkbound(next) .le. 0.0) then
                           nray = nray + 1
                           if (nray .eq. MAX_RAYS) goto 900
                           crossect(nray+1) = 0
                           call get_boundray(ray(1,next),ray(1,iray),amplitude(next),amplitude(iray),
     &                          phase(next),phase(iray),angle(next),angle(iray),
     &                          chkbound0,acttime,kmah(iray),dtray,dthalf,actint,crossect(next),actlayer,
     &                          npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos(actlayer),
     &                          nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),
     &                          v_grid(1,1,1,actlayer),
     &                          boundtmp,ttmp,ray(1,nray),amplitude(nray),phase(nray),angle(nray),kmah(nray),flag,
     &                          crossect(nray+1),thesect(1,nray+1),tsect(nray+1),rtsect(1,nray+1),fictsect(nray+1),
     &                          amplitude(nray+1),phase(nray+1),kmahsect(nray+1),angle(nray+1),
     &                          N_PARAM, n_allint, MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID,N_RTPARAM,inf_dist)
c-----------------------------------------------------
c mindestens ein boundary ray vorhanden:
c
                           if (flag .eq. 0) then
                              call get_z(boundtmp(1),x_int(1,actint),z_int(1,actint),d_int(1,actint),
     &                             c_int(1,actint),b_int(1,actint),npoints_int(actint),
     &                             zint,actp,
     &                             MAXPOINTS_INT)
                              if ((actint .gt. actlayer .and. boundtmp(2) .lt. zint) .or.
     &                             (actint .eq. actlayer .and. boundtmp(2) .gt. zint)) then
                                 newcross = 0
         call do_intersection(boundtmp,actlayer,ttmp,(acttime-ttmp),ihaltup,ihaltdown,veltype,0,
     &        x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos(actlayer),left_border,right_border,
     &        nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),v_grid(1,1,1,actlayer),
     &        ray(1,nray),kmah(nray),chkbound(nray),newcross,
     &        thesect(1,nray),tsect(nray),rtsect(1,nray),kmahsect(nray),crossect(nray),fictsect(nray),
     &        n_int,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

c Erster boundray hat nix geschnitten und wird normal eingebaut:
                                 if (newcross .eq. 0) then
                                    if (ray(1,nray) .gt. left_border .and. ray(1,nray) .lt. right_border) then
                                       if (f_xrays .ne. 0) write(f_xrays,*) boundtmp(1),boundtmp(2)," boundary ray"
c Ausgabe in RAY.TMP:
                                       if (f_rayouttmp .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
                                          if (actcodestep(nactcodestep) .ge. codestep_plot(actcode)) then
                                             write(f_rayouttmp,err=960) nray,boundtmp(1),boundtmp(2)
                                          end if
                                       end if
                                       crossect(nray) = 0
                                       nextray(nray) = -next
                                       nextray(iray) = nray
                                       nextsect(iray) = nray
                                       fictsect(nray) = 0
                                       nextboresect(nray) = next
                                       nextboresect(iray) = nray
                                       do ibore = 1, nboreholes
                                          if (tleftboresect(ibore,iray) .ne. 0.0 .or. tleftboresect(ibore,next) .ne. 0.0) then
                                             tleftboresect(ibore,nray) = -1.0
                                          else
                                             tleftboresect(ibore,nray) = 0.0
                                          end if
                                          if (trightboresect(ibore,iray) .ne. 0.0 .or. trightboresect(ibore,next) .ne. 0.0) then
                                             trightboresect(ibore,nray) = -1.0
                                          else
                                             trightboresect(ibore,nray) = 0.0
                                          end if
                                       end do
c Zweiter boundary ray vorhanden:
                                       if (crossect(nray+1) .ne. 0) then
                                          nextsect(nray) = -(nray+1)
                                          nextsect(nray+1) = next
                                          nextray(nray+1) = 0
                                          nray = nray + 1
                                          if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c     ...nicht vorhanden:
                                       else
                                          nextsect(nray) = -next
                                       end if
                                       irayold = nray
c------------
c Erster boundray wird gekillt (weil er Rand ueberschritten hat):
c Zweiter Boundary vorhanden:
                                    else if (crossect(nray+1) .ne. 0) then
                                       do i = 1, N_PARAM
                                          thesect(i,nray) = thesect(i,nray+1)
                                       end do
                                       do i = 1, N_RTPARAM
                                          rtsect(i,nray) = rtsect(i,nray+1)
                                       end do
                                       amplitude(nray) = amplitude(nray+1)
                                       crossect(nray) = crossect(nray+1)
                                       phase(nray) = phase(nray+1)
                                       kmahsect(nray) = kmahsect(nray+1)
                                       angle(nray) = angle(nray+1)
                                       nray = nray - 1
                                       nextray(iray) = -next
                                       nextray(nray) = 0
                                       nextsect(iray) = -nray
                                       nextsect(nray) = next
                                       if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c Kein boundary ray vorhanden:
                                    else
                                       nray = nray - 1
                                       nextray(iray) = -next
                                       nextsect(iray) = -next
                                    end if
c Erster boundary ray hat spaeter Interface geschnitten, wird als SECT uebernommen:
                                 else
                                    if (f_xrays .ne. 0) write(f_xrays,*) boundtmp(1),boundtmp(2)," boundary ray"
                                    nextray(iray) = -next
                                    nextray(nray) = 0
                                    nextsect(iray) = nray
c Zweiter boundary ray vorhanden:
                                    if (crossect(nray+1) .ne. 0) then
                                       nextray(nray+1) = 0
                                       nextsect(nray) = -(nray+1)
                                       nextsect(nray+1) = next
                                       nray = nray + 1
                                       if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c ...nicht vorhanden:
                                    else
                                       nextsect(nray) = -next
                                       nextray(nray) = 0
                                    end if
                                 end if
c------------------------------------------
c Erster boundray wird gekillt (weil er direkt schon in der gegenueberliegenden Schicht liegt):
c Zweiter Boundary vorhanden:
                              else if (crossect(nray+1) .ne. 0) then
                                 do i = 1, N_PARAM
                                    thesect(i,nray) = thesect(i,nray+1)
                                 end do
                                 do i = 1, N_RTPARAM
                                    rtsect(i,nray) = rtsect(1,nray+1)
                                 end do
                                 amplitude(nray) = amplitude(nray+1)
                                 crossect(nray) = crossect(nray+1)
                                 phase(nray) = phase(nray+1)
                                 kmahsect(nray) = kmahsect(nray+1)
                                 angle(nray) = angle(nray+1)
                                 fictsect(nray) = fictsect(nray+1)
                                 nextray(iray) = -next
                                 nextray(nray) = 0
                                 nextsect(iray) = -nray
                                 nextsect(nray) = next
                                 if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c Kein boundary ray vorhanden:
                              else
                                 nray = nray - 1
                                 nextray(iray) = -next
                                 nextsect(iray) = -next
                              end if
c-----------------------------------------
c Strahl 'iray' schnitt zweimal die Grenzflaeche und wird daher geloescht:
                           else if (flag .lt. 0) then
                              nray = nray -1
                              nextray(irayold) = next*isign(1,nextray(irayold))
                              nextray(iray) = 0
                              if (nextray(irayold) .gt. 0) then
                                 nextsect(irayold) = next
                              end if
                              iray0 = 0
                              do while(abs(nextray(iray0)) .ne. irayold)
                                 iray0 = abs(nextray(iray0))
                              end do
                              irayold = iray0
                              next = abs(nextray(irayold))
c---------------------------------------------
c Abbruch der boundray-Suche wegen caustic zone, bzw. Sztrahl IRAY ist bereits gut genug:
                           else
                              nray = nray - 1
                              nextray(iray) = -next
                              nextsect(iray) = -next
                           end if
                        else
c=======================
c Interpolation 2
c CROSSECT(IRAY) = 0; CROSSECT(NEXT) != 0

                           call do_interpol(actlayer,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &                          spread_max,cos_min,angle_max,next,
     &                          ray,amplitude,phase,angle,kmah,
     &                          crossect,fictsect,kmahsect,nextray,nextsect,irayold,iray,nray,
     &                          nextboresect,tleftboresect,trightboresect,error,
     &                          N_PARAM,n_int,MAXPOINTS_INT,nboreholes,MAX_RAYS)
                           if (error .ne. 0) goto 990
c Ende Interpolation 2
c=======================
                        end if
                     end if
                  else
                     irayold = iray
                  end if
               else if (nextray(iray) .gt. 0) then
                  if (crossect(next) .eq. 0) then
                     actint = abs(crossect(iray))
                     call intersect_wfront(ray(1,iray),ray(1,next),actint,crossect(iray),
     &                    npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos(actlayer),actlayer,
     &                    nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),
     &                    v_grid(1,1,1,actlayer),
     &                    dzdx,raysect,abbruch,
     &                    n_int,MAXPOINTS_INT,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
                     chkbound0 = dzdx*raysect(3) - raysect(4)
c                     if (chkbound0*chkbound(iray).le.0.0 .and. tsect(iray)-tray0().gt.dtray) then
c==========================================
c boundary ray Fall ...

                     if (abbruch .eq. 0 .and. chkbound0*chkbound(iray) .le. 0.0) then
                        nray = nray + 1
                        if ((nray+1) .eq. MAX_RAYS) goto 900
                        crossect(nray+1) = 0
                        call get_boundray(ray(1,iray),ray(1,next),amplitude(iray),amplitude(next),
     &                       phase(iray),phase(next),angle(iray),angle(next),
     &                       chkbound0,acttime,kmah(iray),dtray,dthalf,actint,crossect(iray),actlayer,
     &                       npoints_int,x_int,z_int,b_int,c_int,d_int,iii,iii2,veltype,ptos(actlayer),
     &                       nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),
     &                       v_grid(1,1,1,actlayer),
     &                       boundtmp,ttmp,ray(1,nray),amplitude(nray),phase(nray),angle(nray),kmah(nray),flag,
     &                       crossect(nray+1),thesect(1,nray+1),tsect(nray+1),rtsect(1,nray+1),fictsect(nray+1),
     &                       amplitude(nray+1),phase(nray+1),kmahsect(nray+1),angle(nray+1),
     &                       N_PARAM, n_allint, MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID,N_RTPARAM,inf_dist)
c-------------------------
c mindestens ein boundary ray wurde berechnet:
c
                        if (flag .eq. 0) then
                           call get_z(boundtmp(1),x_int(1,actint),z_int(1,actint),d_int(1,actint),
     &                          c_int(1,actint),b_int(1,actint),npoints_int(actint),
     &                          zint,actp,
     &                          MAXPOINTS_INT)
c-------------------------------
c pruefe auf Schnittpunkt mit oberer und unterer Grenzflaeche:
c
                           if ((actint .gt. actlayer .and. boundtmp(2) .lt. zint) .or.
     &                          (actint .eq. actlayer .and. boundtmp(2) .gt. zint)) then
                              newcross = 0
         call do_intersection(boundtmp,actlayer,ttmp,(acttime-ttmp),ihaltup,ihaltdown,veltype,0,
     &        x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,ptos(actlayer),left_border,right_border,
     &        nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),v_grid(1,1,1,actlayer),
     &        ray(1,nray),kmah(nray),chkbound(nray),newcross,
     &        thesect(1,nray),tsect(nray),rtsect(1,nray),kmahsect(nray),crossect(nray),fictsect(nray),
     &        n_int,MAXPOINTS_INT,N_PARAM,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)

c--------------------
c Erster boundray hat nichts geschnitten und wird normal eingebaut:
                              if (newcross .eq. 0) then
                                 if (ray(1,nray) .gt. left_border .and. ray(1,nray) .lt. right_border) then
                                    if (f_xrays .ne. 0) write(f_xrays,*) boundtmp(1),boundtmp(2)," boundary ray"
c Ausgabe in RAY.TMP:
                                    if (f_rayouttmp .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
                                       if (actcodestep(nactcodestep) .ge. codestep_plot(actcode)) then
                                          write(f_rayouttmp,err=960) nray,boundtmp(1),boundtmp(2)
                                       end if
                                    end if
                                    crossect(nray) = 0
                                    nextray(nray) = next
                                    nextsect(nray) = next
                                    nextray(iray) = 0
                                    nextray(irayold) = -nray
                                    fictsect(nray) = 0
                                    nextboresect(nray) = next
                                    nextboresect(iray) = nray
                                    do ibore = 1, nboreholes
                                       if (tleftboresect(ibore,iray) .ne. 0.0 .or. tleftboresect(ibore,next) .ne. 0.0) then
                                          tleftboresect(ibore,nray) = -1.0
                                       else
                                          tleftboresect(ibore,nray) = 0.0
                                       end if
                                       if (trightboresect(ibore,iray) .ne. 0.0 .or. trightboresect(ibore,next) .ne. 0.0) then
                                          trightboresect(ibore,nray) = -1.0
                                       else
                                          trightboresect(ibore,nray) = 0.0
                                       end if
c                                       tleftboresect(ibore,nray) = 0.0
c                                       trightboresect(ibore,nray) = 0.0
                                    end do
c Zweiter boundary ray vorhanden:
                                    if (crossect(nray+1) .ne. 0) then
                                       nextsect(iray) = nray+1
                                       nextsect(nray+1) = -nray
                                       nextray(nray+1) = 0
                                       nray = nray + 1
                                       if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c ...nicht vorhanden:
                                    else
                                       nextsect(iray) = -nray
                                    end if
                                    irayold = nray
c--------------------------------------
c Erster boundray wird gekillt (weil er Rand ueberschritten hat):
c Zweiter Boundary vorhanden:
                                 else if (crossect(nray+1) .ne. 0) then
                                    do i = 1, N_PARAM
                                       thesect(i,nray) = thesect(i,nray+1)
                                    end do
                                    do i = 1, N_RTPARAM
                                       rtsect(i,nray) = rtsect(i,nray+1)
                                    end do
                                    amplitude(nray) = amplitude(nray+1)
                                    crossect(nray) = crossect(nray+1)
                                    phase(nray) = phase(nray+1)
                                    kmahsect(nray) = kmahsect(nray+1)
                                    angle(nray) = angle(nray+1)
                                    fictsect(nray) = fictsect(nray+1)
                                    nray = nray - 1
                                    nextsect(iray) = nray
                                    nextsect(nray) = -next
                                    nextray(irayold) = -next
                                    nextray(iray) = 0
                                    nextray(nray) = 0
                                    if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c Kein boundary ray vorhanden:
                                 else
                                    nray = nray - 1
                                    nextsect(iray) = -next
                                    nextray(iray) = 0
                                    nextray(irayold) = -next
                                 end if
                                 iray = irayold
c------------------------
c Erster boundary ray hat spaeter Interface geschnitten, wird als SECT uebernommen:
                              else
                                 if (f_xrays .ne. 0) write(f_xrays,*) boundtmp(1),boundtmp(2)," boundary ray"
                                 nextray(irayold) = -next
                                 nextray(iray) = 0
                                 nextray(nray) = 0
                                 nextsect(nray) = next
c Zweiter boundary ray vorhanden:
                                 if (crossect(nray+1) .ne. 0) then
                                    nextsect(iray) = nray+1
                                    nextsect(nray+1) = -nray
                                    nextray(nray+1) = 0
                                    nray = nray + 1
                                    if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c ...nicht vorhanden:
                                 else
                                    nextsect(iray) = -nray
                                 end if
                                 iray = irayold
                              end if
c-------------------------------------------
c Erster boundray wird gekillt (weil er direkt schon in der gegenueberliegenden Schicht liegt):
c Zweiter Boundary vorhanden:
                           else if (crossect(nray+1) .ne. 0) then
                              do i = 1, N_PARAM
                                 thesect(i,nray) = thesect(i,nray+1)
                              end do
                              do i = 1, N_RTPARAM
                                 rtsect(i,nray) = rtsect(1,nray+1)
                              end do
                              amplitude(nray) = amplitude(nray+1)
                              crossect(nray) = crossect(nray+1)
                              phase(nray) = phase(nray+1)
                              kmahsect(nray) = kmahsect(nray+1)
                              angle(nray) = angle(nray+1)
                              nextsect(iray) = nray
                              nextsect(nray) = -next
                              nextray(irayold) = -next
                              nextray(iray) = 0
                              nextray(nray) = 0
                              if (f_xrays .ne. 0) write(f_xrays,*) thesect(1,nray),thesect(2,nray)," boundary ray"
c Kein boundary ray vorhanden:
                           else
                              nray = nray - 1
                              nextsect(iray) = -next
                              nextray(iray) = 0
                              nextray(irayold) = -next
                           end if
                           iray = irayold
c------------------------------------
c Strahl 'next' schnitt zweimal die Grenzflaeche und wird daher geloescht:
                        else if (flag .lt. 0) then
                           nray = nray -1
                           nextray(iray) = nextray(next)
                           nextray(next) = 0
                           if (nextray(iray) .gt. 0) then
                              nextsect(iray) = nextray(iray)
                           end if
                           next = iray
                           iray = irayold
c---------------------------------------------------
c Abbruch der boundray-Suche wegen caustic zone, bzw. Strahl NEXT ist bereits OK
                        else
                           nray = nray - 1
                           nextsect(iray) = -next
                           nextray(iray) = 0
                           nextray(irayold) = -next
                        end if
c boundary ray Fall... Ende
c==========================================
c=======================
c Interpolation 3
c CROSSECT(IRAY) != 0; CROSSECT(NEXT) = 0
                     else
                        call do_interpol(actlayer,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &                       spread_max,cos_min,angle_max,next,
     &                       ray,amplitude,phase,angle,kmah,
     &                       crossect,fictsect,kmahsect,nextray,nextsect,irayold,iray,nray,
     &                       nextboresect,tleftboresect,trightboresect,error,
     &                       N_PARAM,n_int,MAXPOINTS_INT,nboreholes,MAX_RAYS)

                        if (error .ne. 0) goto 990
c Ende Interpolation 3
c=======================
                     end if
                  else
                     irayold = iray
                  end if
               else
                  irayold = iray
               end if
            end if
c End if kmah
         end do

c-------------------------------
c Ausgabe in RAY.TMP:
c
         if (f_rayouttmp .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
            if (actcodestep(nactcodestep) .ge. codestep_plot(actcode).and.timesteps .le.last_wfronts) then
               if (nray .gt. nrayold) then
                  do iray = nrayold+1, nray
                     write(f_rayouttmp,err=960) iray,ray(1,iray),ray(2,iray)
                  end do
               end if
            end if
         end if
c----------------------------------------
c Ausgabe der Wellenfronten in F_WFRONTS:
c
         if (f_wfronts .ne. 0) then
            if (step_wfronts .ne. 0 .and.
     &           timesteps .ge. step1_wfronts .and. actcodestep(nactcodestep) .ge. codestep_plot(actcode)) then
               if (mod(timesteps,step_wfronts) .eq. 0 .and. timesteps .le. last_wfronts) then
                  iray = 0
                  i = 0
                  do while  (nextray(iray) .ne. 0)
                     i = i + 1
                     irayold = iray
                     iray = abs(nextray(iray))
                     if (crossect(iray) .eq. 0) then
                        if (mod(i,step_spread) .eq. 0) then
                           write(f_wfronts,*,err=970) ray(1,iray),ray(2,iray)
                           if (nextray(iray) .lt. 0) then
                              write(f_wfronts,'()',err=970)
                           end if
                        else if (nextray(iray) .lt. 1 .or. nextray(irayold) .lt. 1) then
                           write(f_wfronts,*,err=970) ray(1,iray),ray(2,iray)
                           i = 0
                           if (nextray(iray) .lt. 0) then
                              write(f_wfronts,'()',err=970)
                           end if
                        else if (crossect(nextray(iray)) .ne. 0 .or. crossect(irayold) .ne. 0) then
                           write(f_wfronts,*,err=970) ray(1,iray),ray(2,iray)
                           i = 0
                        else if (kmah(iray) .ne. kmah(irayold)) then
                           write(f_wfronts,*,err=970) ray(1,iray),ray(2,iray)
                           i = 0
                        else if (nextray(iray) .gt. 0 .and. kmah(iray) .ne. kmah(nextray(iray))) then
                           write(f_wfronts,*,err=970) ray(1,iray),ray(2,iray)
                           i = 0
                        end if
                     else
                        write(f_wfronts,'()',err=970)
                     end if
                  end do
                  write(f_wfronts,'()',err=970)
               end if
            end if
         end if

c Ende Interpolation
c++++++++++++++++++++++++++++++++++++++++++++

      end do
c Ende aeussere Schleife... naechster Zeitschritt

c--------------------------------------------------------------------------------
c
c Eine Schicht ist damit durchlaufen...
c
c==================================================================================================
c Schreiben der Strahlen in Ausgabedatei:
c
      if (f_rayout .ne. 0 .and. actcodestep(nactcodestep) .ge. codestep_plot(actcode) .and.
     &     (flag_2point .eq. 0 .or. process_2point ) ) then
c         if( .not. process_2point ) write(*,"('Writing rays in output file ',a30)") fname_rayout
         do i = 1, nray
            rewind(f_rayouttmp)
            iray = 0
 210        read(f_rayouttmp,end = 211,err=962) n,x,z
            if (n .eq. i) then
               if (iray .lt. MAX_RAYS) then
                  iray = iray + 1
                  ray(1,iray) = x
                  ray(2,iray) = z
               else if (iray .lt. 2*MAX_RAYS) then
                  iray = iray + 1
                  ray(3,iray-MAX_RAYS) = x
                  ray(4,iray-MAX_RAYS) = z
               else
                  write(*,*) "WARNING: Too many timesteps! Couldn't write all rays to output file!"
                  goto 211
               end if
            end if
            goto 210
 211        if (iray .ne. 0) then
               if (iray .le. MAX_RAYS) then
                  do iray0 = 1, iray
                     write(f_rayout,*,err=965) ray(1,iray0),ray(2,iray0)
                  end do
               else
                  do iray0 = 1, MAX_RAYS
                     write(f_rayout,*,err=965) ray(1,iray0),ray(2,iray0)
                  end do
                  do iray0 = 1, iray-MAX_RAYS
                     write(f_rayout,*,err=965) ray(3,iray0),ray(4,iray0)
                  end do
               end if
               write(f_rayout,'()',err=965)
            end if
         end do
         rewind(f_rayouttmp)
         write(f_rayouttmp,err=960) n,x,z
         rewind(f_rayouttmp)
      end if
c
c Schreibe Strahlen Ende
c--------------------------------------------------------------------------------

c      iray = 0
c      do while(nextsect(iray) .ne. 0)
c         iray = abs(nextsect(iray))
c         if (crossect(iray) .ne. 0 .and. iray .le. nray) then
c            write(89,*) sect(1,iray),sect(2,iray),nextsect(iray),fictsect(iray)
c         end if
c      end do
c--------------------------------------------------------------------------------
c Aufsammeln der Schnittpunkte SECT in entsprechenden HALTS (sowie viele neue boundary rays interpolieren)
c
      call sortsect(thesect,tsect,rtsect,amplitude,phase,angle,kmahsect,crossect,ihaltup,ihaltdown,
     &     dtray,actlayer,left_border,right_border,nray,spread_max2,inf_dist,
     &     x_int,z_int,d_int,c_int,b_int,npoints_int,iii,iii2,f_out,f_xrays,flag_smooth,
     &     nx_grid(actlayer),nz_grid(actlayer),x_grid(1,actlayer),z_grid(1,actlayer),v_grid(1,1,1,actlayer),veltype,
     &     ptos(actlayer),
     &     nextsect,fictsect,error,
     &     N_PARAM,N_RTPARAM,MAX_RAYS,n_allint,MAXPOINTS_INT,MAXP_XGRID,MAXP_ZGRID)

      if (error .ne. 0) goto 990

 222  call collect_halt(thesect,tsect,rtsect,amplitude,phase,angle,kmahsect,nextsect,crossect,fictsect,actlayer,
     &     code,actcode,actcodestep,nactcodestep,actcompounds,actncomp,hold,
     &     thehalt,thalt,rthalt,amplhalt,phasehalt,anglehalt,kmahhalt,nexthalt,ficthalt,nrayhalt,ncodehalt,
     &     codestepshalt,inthalt,compoundshalt,ncomphalt,ihaltup,ihaltdown,ihaltcomp,error,
     &     N_PARAM,MAX_RAYS,MAX_HALTS,ncodes,MAX_CODESTEPS,N_RTPARAM)
      if (error .ne. 0) goto 990
c=====================================================================================================
c Aufnahme der Einsaetze bei den receivern
c
      do iline = 1, nreclines
         acthalt = 0
         if (ihaltup .ne. 0) then
            if (inthalt(ihaltup) .eq. int_rec(iline)) then
               acthalt = ihaltup
            end if
         end if
         if (ihaltdown .ne. 0) then
            if (inthalt(ihaltdown) .eq. int_rec(iline)) then
               acthalt = ihaltdown
            end if
         end if

         if (acthalt .ne. 0) then
            xzindex = 1

c            print*,"-------------record---------------"
            if (f_timeout .ne. 0 .and. flag_stack .eq. 0) then
               icompout(iline) = icompout(iline) + 1
            end if

            if (flag_stack .ne. 0) then
               iray = 0
               do while (nexthalt(iray,acthalt) .ne. 0)
                  iray = abs(nexthalt(iray,acthalt))
                  narrivals(isource,1) = narrivals(isource,1) + 1
                  arrival(1,isource,narrivals(isource,1),1) = thalt(iray,acthalt)
                  arrival(6,isource,narrivals(isource,1),1) = anglehalt(iray,acthalt)
c OBS! Nur zum Test gedacht, sollte spaeter geaendert werden!
c                  arrival(1,isource,narrivals(isource,1),1) = anglehalt(iray,acthalt)
c                  arrival(6,isource,narrivals(isource,1),1) = 0.0
                  ampltmp = amplhalt(iray,acthalt) * sqrt( v_source / (abs(thehalt(5,iray,acthalt)) * rho_source))
                  phasetmp = phasehalt(iray,acthalt) - 0.5 * pi * kmahhalt(iray,acthalt)
c---------------------- Amplitude ---
c     p-wave (V_REC steckt bereits in AMPLTMP drin)
                  if (veltype_source .gt. 0) then
                     arr2 = ampltmp * thehalt(3,iray,acthalt)
                     arr3 = ampltmp * thehalt(4,iray,acthalt)
c s-wave
                  else
                     arr2 = ampltmp * thehalt(4,iray,acthalt)
                     arr3 = -ampltmp * thehalt(3,iray,acthalt)
                  end if
c---------------------- Phase --------
c
                  if (arr2 .lt. 0.0) then
                     arr2 = -arr2
                     arr4 = phasetmp - pi
                  else
                     arr4 = phasetmp
                  end if
                  if (arr3 .lt. 0.0) then
                     arr3 = -arr3
                     arr5 = phasetmp - pi
                  else
                     arr5 = phasetmp
                  end if
                  arrival(2,isource,narrivals(isource,1),1) = arr2
                  arrival(3,isource,narrivals(isource,1),1) = arr3
                  arrival(4,isource,narrivals(isource,1),1) = arr4
                  arrival(5,isource,narrivals(isource,1),1) = arr5

               end do
               
               if (f_timeout .ne. 0) then
                  do i = 1, narrivals(isource,1)
                     write(f_timeout,1080,err=999) x_source,arrival(1,isource,i,1),arrival(2,isource,i,1),arrival(3,isource,i,1)
                     time(isource,i,1) = arrival(1,isource,i,1)
                     amp_out(isource,i,1) = arrival(comp_out,isource,i,1)
                     phase_out(isource,i,1) = arrival(comp_out+2,isource,i,1)
                  end do
               end if
            else if ( .not. process_2point ) then
c               write(*,*) "-----------------------------------------------"
c               write(*,*) " Start record "
               if( act_timecode .gt. ntimecodes ) then
                  goto 975
               end if
               call record(thehalt(1,1,acthalt),thalt(1,acthalt),amplhalt(1,acthalt),phasehalt(1,acthalt),
     &              anglehalt(1,acthalt),kmahhalt(1,acthalt),nexthalt(0,acthalt),ficthalt(1,acthalt),xzindex,
     &              xz_rec(1,iline),nrec(iline),int_rec(iline),lay_rec(iline),iii2,f_timeout,
     &              time(1,1,act_timecode), amp_out(1,1,act_timecode),
     &              phase_out(1,1,act_timecode), angle_out(1,1,act_timecode), comp_out, ccp(1,1,act_timecode),
     &              actlayer, npoints_int,x_int,z_int,b_int,c_int,d_int,veltype, ptos,
     &              nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,flag_smooth,
     &              arrival(1,1,1,iline),narrivals(1,iline),error,
     &              n_allint,MAXPOINTS_INT,MAX_RECEIVERS,MAX_ARRIVALS,MAX_RAYS,N_PARAM,N_RTPARAM,N_ARRPARAM,
     &              MAXP_XGRID,MAXP_ZGRID)
c              write(*,*) "TIME code ",actcode,act_timecode
               act_timecode = act_timecode + 1
               if( flag_2point .ne. 0 .and. .not. process_2point ) then
                  nangles_2point(actcode) = 0
                  do k = 1, nrec(iline)
c Only save first 3 arrivals (or less if not available)
                     do i = 1, min( 3, narrivals(k,iline) )
                        nangles_2point(actcode) = nangles_2point(actcode)+1
                        angles_2point( nangles_2point(actcode) ,actcode) = arrival(6,k,i,iline)
                     enddo
                  enddo
c                  call sort_array( angles_2point(1,actcode), nangles_2point(actcode) )
               endif
            end if
c-----------------------------
c Ausgabe in f_seismout
c
c            if (f_seismout(1) .ne. 0) then
c               icode = ncodehalt(acthalt)
c               if (code(1,icode) .lt. 0) then
c                  icode = -icode
c               end if
c               do k = 1, nrec(iline)
c                  do j = 1, narrivals(k,iline)
c                     write(f_seismout(iline),1050,err=940)
c     &                    icode,k,arrival(1,k,j,iline),arrival(2,k,j,iline),ay,arrival(3,k,j,iline),
c     &                    arrival(4,k,j,iline),phy,arrival(5,k,j,iline),arrival(6,k,j,iline),tast
c                  end do
c               end do
c            end if
c ARRIVAL(1,irec,a)   : Laufzeit an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(2,irec,a)   : Amplitude(x) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(3,irec,a)   : Amplitude(z) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(4,irec,a)   : Phase(x) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(5,irec,a)   : Phase(z) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(6,irec,a)   : Ausgangswinkel von der Quelle an receiver 'irec', Anzahl der Einsaetze 'a'
c---------------------------------------
c Ausgabe in f_codeout, Vorbereitung von f_timeout
c
c ..veraendert..
c

            if (f_timeout .ne. 0) then
               k = 0
               flag = 0
               do while (k .lt. nrec(iline) .and. flag .eq. 0)
                  k = k + 1
                  if (narrivals(k,iline) .ne. 0) flag = -1
               end do
c Es sind Einsaetze aufgenommen worden:
               if (flag .ne. 0) then
                  write(f_timeout,'()',err=925)
               else
                  icompout(iline) = icompout(iline) - 1
               end if
            end if

            if (error .ne. 0) goto 990
         end if
c----------------------------------------------------
c Ausgabe in f_amplout
c
c ...geloescht...
c
         do k = 1, nrec(iline)
            narrivals(k,iline) = 0
         end do
      end do
c
c (Ende Aufnahme)
c===================================================================
c Aufnahme der Einsaetze an den Bohrloechern
c
      lay_bore = actlayer
      xzindex = 2
      do ibore = 1, nboreholes
         call collect_bore(leftboresect,tleftboresect,kmahlboresect,rightboresect,trightboresect,kmahrboresect,
     &        nextboresect,ibore,amplitude,phase,angle,
     &        actcode,actcodestep,nactcodestep,actcompounds,actncomp,
     &        borehalt,tborehalt,amplborehalt,phaseborehalt,angleborehalt,kmahborehalt,nextborehalt,
     &        nrayborehalt,ncodeborehalt,codestborehalt,compborehalt,ncompborehalt,flag,
     &        N_PARAM,MAX_RAYS,MAX_CODESTEPS,nboreholes)
         if (flag .eq. 0) then
            int_bore = ibore + n_int
            iline = ibore + nreclines
c fictsect ist hier ein dummy!(?)
            call record(borehalt,tborehalt,amplborehalt,phaseborehalt,
     &           angleborehalt,kmahborehalt,nextborehalt,fictsect,xzindex,
     &           xz_rec(1,iline),nrec(iline),int_bore,lay_bore,iii2,0,
     &           time(1,1,1), amp_out(1,1,1), phase_out(1,1,1), angle_out(1,1,1), comp_out, ccp(1,1,1),
     &           actlayer, npoints_int,x_int,z_int,b_int,c_int,d_int,veltype, ptos,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,flag_smooth,
     &           arrival(1,1,1,iline),narrivals(1,iline),error,
     &           n_allint,MAXPOINTS_INT,MAX_RECEIVERS,MAX_ARRIVALS,MAX_RAYS,N_PARAM,N_RTPARAM,N_ARRPARAM,
     &           MAXP_XGRID,MAXP_ZGRID)
            if (error .ne. 0) goto 990
c-----------------------------
c Ausgabe in f_seismout
c
c            if (f_seismout(1) .ne. 0) then
c               icode = actcode
c               if (code(1,actcode) .lt. 0) then
c                  icode = -icode
c               end if
c               do k = 1, nrec(iline)
c                  do j = 1, narrivals(k,iline)
c                     write(f_seismout(iline),1050,err=940)
c     &                    icode,k,arrival(1,k,j,iline),arrival(2,k,j,iline),ay,arrival(3,k,j,iline),
c     &                    arrival(4,k,j,iline),phy,arrival(5,k,j,iline),arrival(6,k,j,iline),tast
c                  end do
c               end do
c            end if
c NEW: Just dump all times into f_timeout:
            if (f_timeout .ne. 0) then
               do k = 1, nrec(iline)
                  do j = 1, narrivals(k,iline)
                     write(f_timeout,1080,err=999)
     &                    xz_rec(k,iline), arrival(1,k,j,iline), arrival(2,k,j,iline), arrival(3,k,j,iline)
                     time(k,j,1) = arrival(1,k,j,iline)
                     amp_out(k,j,1) = arrival(comp_out,k,j,iline)
                     phase_out(k,j,1) = arrival(comp_out+2,k,j,iline)
                  end do
               end do
            end if

c-----------------------------
c Loeschen aller Einsaetze
            do k = 1, nrec(iline)
               narrivals(k,iline) = 0
            end do
c-----------------------------
         else
            write(*,*) "no borehole hits, layer ",actlayer
         end if

      end do

c=====================================================================================================
c Freistellen der halts, falls code ganz durchgelaufen ist:
c
      if (ihaltup .ne. 0) then
         if (ihaltup .eq. ihaltcomp .and. hold(actcodestep(nactcodestep),actcode) .gt. 1) then
            if (codestepshalt(ihaltup) .eq. codesteps(hold(actcodestep(nactcodestep),actcode))) then
               codestepshalt(ihaltup) = 0
            else
               ihaltup = 0
            end if
         else if (codestepshalt(ihaltup) .eq. codesteps(actcode) .and.
     &        hold(codestepshalt(ihaltup),actcode) .eq. 0) then
            codestepshalt(ihaltup) = 0
         end if
      end if
      if (ihaltdown .ne. 0) then
         if (ihaltdown .eq. ihaltcomp .and. hold(actcodestep(nactcodestep),actcode) .gt. 1) then
            if (codestepshalt(ihaltdown) .eq. codesteps(hold(actcodestep(nactcodestep),actcode))) then
               codestepshalt(ihaltdown) = 0
            else
               ihaltdown = 0
            end if
         else if (codestepshalt(ihaltdown) .eq. codesteps(actcode) .and.
     &        hold(codestepshalt(ihaltdown),actcode) .eq. 0) then
            codestepshalt(ihaltdown) = 0
         end if
      end if
c---------------------------------------------------------------------
c Bestimmung der naechsten Halts sowie weiterer entsprechender Werte
c Bestimmung des naechsten Codeschritts oder sogar Codes
c
      oldcode = actcode
      call get_acthalt(code,codesteps,ncodehalt,ncodes,ihaltup,ihaltdown,
     &     ncomphalt,compoundshalt,hold,codestepshalt,flag_compounds,all_reclines,
     &     actcode,actcodestep,nactcodestep,acthalt,actualhalt,ihaltnorm,ihaltcomp,actncomp,
     &     actcompounds,icompout,flag,
     &     MAX_HALTS, ncodes, MAX_CODESTEPS, nreclines)

c FLAG != 0: Letzter Code erreicht, gehe an's Ende des Programms
      if (flag .ne. 0) then
         if( flag_2point .eq. 0 .or. process_2point ) then
            goto 200            
         else
            process_2point = .true.
            goto 111
         endif
      endif

c Trick to reset halts so that new code is started from beginning (no previous halts possible)
      if( oldcode .ne. actcode .and. process_2point ) then
         acthalt = 0
         actcodestep(1) = 1
         actncomp = 0
      endif

      if (actcodestep(nactcodestep) .gt. 1) then
         oldlayer = code(actcodestep(nactcodestep)-1,actcode)
         oldveltype = isign(1,oldlayer)
         oldlayer = abs(oldlayer)
      end if
      actlayer = abs(code(actcodestep(nactcodestep),actcode))
      veltype = isign(1,code(actcodestep(nactcodestep),actcode))
c===========================================================================================
c Fuer acthalt=0 neuer Beginn; sonst Reflexion/Transmission
c
c Fall: Es wird von der Quelle aus propagiert:
c
      if (acthalt .eq. 0) then

         if (actlayer .eq. int_source .or. code(0,actcode) .gt. 0) then
            ihaltdown = ihaltnorm
            ihaltup = ihaltcomp
         else
            ihaltup = ihaltnorm
            ihaltdown = ihaltcomp
         end if
         nray = nray_source
         if( process_2point ) then
            nray = nangles_2point(actcode)
         endif
c------------------------------------------------------------------
c Quelle liegt innerhalb Schicht:
c (nichts aendert sich)
c------------------------------------------------------------------
c Quelle liegt auf Interface:
c (d_angle aendert sich nicht)
         if (int_source .ne. 0) then
            angle_source = sign(1.0,dzdx_source) * acos(1.0/sqrt(1.0 + dzdx_source*dzdx_source))
c----------------------
c kein beam:
            if (anglein .eq. 0.0 .and. angleout .eq. 0.0) then
               angle_source = angle_source + inf_angle
c Fall: Strahl laeuft nach oben:
               if (actlayer .lt. int_source) then
                  angle_source = angle_source + pi
               end if
c----------------------
c beam, Einschraenkung durch ANGLEIN, ANGLEOUT
            else
c     Fall: Strahl laeuft nach oben:
               if (actlayer .lt. int_source) then
                  angle_source = angle_source + pi
               end if
               if (anglein .lt. angle_source .or. anglein .gt. angle_source+pi) then
c Falls beam in die andere Schicht zeigt, hier kein Strahl!:
                  nray = 0
               else
                  angle_source = anglein
               end if
            end if
         end if
c--------------------------
         call velocity(x_source, z_source, nx_grid(lay_source), nz_grid(lay_source),
     &        x_grid(1,lay_source), z_grid(1,lay_source), v_grid(1,1,1,lay_source),
     &        veltype, ptos(lay_source),lay_source,
     &        vel,
     &        MAXP_XGRID, MAXP_ZGRID)
         v_source = vel(1)
         veltype_source = veltype

         acthalt = 0
         actualhalt(1) = ihaltnorm
         do i = 1, nray
            tray0(i) = t_source
            ray0(1,i) = x_source
            ray0(2,i) = z_source

            nextray0(i) = i + 1
            nextray(i) = 0
            nextsect(i) = nextray0(i)

            if( process_2point ) then
               angle(i) = angles_2point(i,actcode)
c Make next ray "unconnected" to this ray, to avoid interpolation of new rays:
               nextsect(i) = -nextray0(i)
            else
               angle(i) = angle_source + d_angle*float(i-1)
            endif

            ray0(3,i) = cos(angle(i))/vel(1)
            ray0(4,i) = sin(angle(i))/vel(1)
            ray0(5,i) = 0.0
            ray0(6,i) = 1.0
            ray0(7,i) = x_source
            ray0(8,i) = z_source

            kmah0(i) = 0
            crossect(i) = 0
            fictsect(i) = 0
            amplitude(i) = vel(1)*sqrt(rho_source)
            phase(i) = 0.0
            nextboresect(i) = nextray0(i)
            do ibore = 1, nboreholes
               tleftboresect(ibore,i) = 0.0
               trightboresect(ibore,i) = 0.0
            end do
         end do

c Zur Sicherheit, damit der erste und letzte Strahl in diesem Falle wirklich gleich sind:
         if (int_source .eq. 0 .and. anglein .eq. 0.0 .and. angleout .eq. 0.0 .and. .not. process_2point ) then
            do i = 3, 4
               ray0(i,nray) = ray0(i,1)
            end do
         end if

c---------
c Quelle liegt auf conciding interface:
c Noch ist nur Erdoberflaeche als Interface moeglich...
c
         if (int_source .eq. 1) then
            call get_actp(x_source,int_source,x_int,npoints_int,
     &           actp,
     &           n_allint,MAXPOINTS_INT)
            if (iii2(actp,int_source) .ne. 0) then
               ficttmp = iii2(actp,int_source)*1000 - 1000
               do iray = 1, nray
                  fictsect(iray) = ficttmp
               end do
            end if
         end if

         nextray0(0) = 1
         nextray0(nray) = 0
         nextray(0) = 0
         nextsect(0) = -1
         nextsect(nray) = 0
         actncomp = 0
         nray0 = nray
         nextboresect(0) = -1
         nextboresect(nray) = 0
         acttime = t_source
c===========================================================================================
c Transmission / Reflexion
      else
         if (actlayer .eq. inthalt(acthalt)) then
            ihaltup = ihaltcomp
            ihaltdown = ihaltnorm
            updownflag = 1
         else
            ihaltdown = ihaltcomp
            ihaltup = ihaltnorm
            updownflag = -1
         end if
c---------
c Fall: Reflexion
         if (abs(code(actcodestep(nactcodestep)-1,actcode)) .eq. actlayer) then
            signrt = -1
c Fall: Transmission
         else
            signrt = +1
         end if
         call refstack(thehalt(1,1,acthalt),thalt(1,acthalt),rthalt(1,1,acthalt),amplhalt(1,acthalt),
     &        phasehalt(1,acthalt),anglehalt(1,acthalt),inthalt(acthalt),nrayhalt(acthalt),kmahhalt(1,acthalt),
     &        ficthalt(1,acthalt),nexthalt(0,acthalt),signrt,x_int,z_int,d_int,c_int,b_int,npoints_int,
     &        ptos,actlayer,oldlayer,veltype,oldveltype,inf_angle,inf_dist,updownflag,flag_surface,
     &        nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,flag_stack,
     &        ray0,tray0,amplitude,phase,angle,kmah0,nextsect,nray0,error,f_xrays,
     &        crossect,fictsect,process_2point,ccp(1,1,actcode),MAX_RECEIVERS,
     &        n_allint,MAXPOINTS_INT,N_PARAM,MAX_RAYS,N_RTPARAM,MAXP_XGRID,MAXP_ZGRID,reflpoint)

         if (error .ne. 0) goto 990

c R/T fand nicht statt aufgrund ueberkritischen Winkeln... zurueck
         if (nray0 .eq. 0) then
            write(*,*) "Code versiegt, R/T findet nicht statt"
            nextsect(0) = 0
            goto 222
         end if
         acttime = 1000000.0
         irayold = 0
         iray = 0
         do while (nextsect(iray) .ne. 0)
            iray = abs(nextsect(iray))
            nextray0(irayold) = iray
            irayold = iray
            nextray(iray) = 0
            nextboresect(iray) = nextsect(iray)
            do ibore = 1, nboreholes
               tleftboresect(ibore,iray) = 0.0
               trightboresect(ibore,iray) = 0.0
            end do
            if (tray0(iray) .lt. acttime) acttime = tray0(iray)
         end do
         nextray0(irayold) = 0
         nextboresect(irayold) = 0
      end if
c
c-------------------------------------------------------------------

      nray = nray0
      acttime = float(int(acttime/dtray))*dtray
      timesteps = int(acttime/dtray)

      if (f_rayout .ne. 0 .and. (flag_2point .eq. 0 .or. process_2point ) ) then
         if (actcodestep(nactcodestep) .ge. codestep_plot(actcode)) then
            iray = 0
            do while (nextray0(iray) .ne. 0)
               iray = abs(nextray0(iray))
               if (fictsect(iray) .eq. 0) then
                  write(f_rayouttmp,err=960) iray,ray0(1,iray),ray0(2,iray)
               end if
            end do
         end if
      end if

      goto 100
c
c Travel times and amplitudes computation end
c=============================================================================================
c=============================================================================================
c Optional seismogram plot
c

 200  continue

      if (flag_stack .ne. 0) then
         if (isource .lt. nrec(1)) then
            isource = isource + 1
            x_source = xz_rec(isource,1)
            call get_z(x_source,x_int(1,int_source),z_int(1,int_source),d_int(1,int_source),
     &           c_int(1,int_source),b_int(1,int_source),npoints_int(int_source),
     &           z_source,actp,
     &           MAXPOINTS_INT)
            goto 111
         end if
      end if

c---------------------------------------------------------------------------------------------
c Screen output
c

c Halbierung des Zeitschritts, Ausgabe der Anzahl der Halbierungen:
      if (halfsteps_count .gt. 0) then
         write(f_out,1060,err=915) halfsteps_count
      end if

      if (f_timeout .ne. 0) then
         close(f_timeout)
         if (j .eq. MAX_RAYS) then
            write(f_out,*,err=915) "Travel time output file doesn't really contain all travel times due to array overflow"
         end if
      end if

      goto 999


c 1001 format(" actcodestep,ihaltup,ihaltdown,nray: ",i3,i3,i3,i5)
 1060 format(/" WARNING: Halfing of the time interval took place ",i10," times during the computation"/,
     &     "to get an acceptable accuracy with the Runge Kutta method.")
 1080 format(1x, f10.5, 3xf12.8)

 3000 format(" WFRONT: Open error, file no. ",i2,", name ",a30)
 3010 format(" WFRONT: Write error, file no. ",i2,", name ",a30)
 3020 format(" WFRONT: Read error, file no. ",i2,", name ",a30)
 3100 format(/," WFRONT: Problems with ray beam:",/,
     &     " A beam may only point into ONE layer. In THIS case, you can choose ANGLEIN and ANGLEOUT between",/,
     &     1x,f10.7," and ",f10.7,"       or between ",f10.7," and ",f10.7,".",/)

 800  if (angle_source .gt. pi) then
         angle_source = angle_source - pi
      end if
      write(*,3100) angle_source,angle_source+pi,angle_source+pi,angle_source+2.*pi
      goto 990
 900  write(*,*) "WFRONT : Too many rays!"
      write(*,*) "Enlarge the PARAMETER 'MAX_RAYS' in WFRONT.F"
      goto 990
 915  write(*,3010) f_out,fname_out
      goto 990
 920  write(*,3000) f_timeout,fname_timeout
      goto 990
 925  write(*,3010) f_timeout,fname_timeout
      goto 990
 960  write(*,3010) f_rayouttmp,fname_rayouttmp
      goto 990
 962  write(*,3020) f_rayouttmp,fname_rayouttmp
      goto 990
 965  write(*,3010) f_rayout,fname_rayout
      goto 990
 970  write(*,3000) f_wfronts,fname_wfronts
      goto 990
 975  write(*,*) "Too many codes recorded > ", ntimecodes
      goto 990
 990  write(*,*) "...program terminated on error"
      goto 999
 995  write(*,*) "...program terminated"
      goto 999

 999  end


