c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm RECORD / Hauptprogramm WFRONT
c letzte Aenderung 6.02.97 / nicht getestet. Es fehlt: Ueberpruefung der Abfrage des KMAH-Index
c Wird nur von WFRONT aufgerufen
c Aufruf des Unterprogramms GET_RECORD
c Zeichnet ankommende Signale eines durchlaufenen Strahlencodes auf.
c Dazu wird zwischen den auf der Grenzflaeche angekommenen Strahlen
c paraxial approximiert. Ergebnis:
c ARRIVAL(1,irec,a)   : Laufzeit an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(2,irec,a)   : Amplitude(x) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(3,irec,a)   : Amplitude(z) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(4,irec,a)   : Phase(x) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(5,irec,a)   : Phase(z) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(6,irec,a)   : Ausgangswinkel von der Quelle an receiver 'irec', Anzahl der Einsaetze 'a'
c
      subroutine record(ray,tray,amplitude,phase,angle,kmah,nextray,fict,xzindex,
     &     xz_rec,nrec,int_rec,lay_rec,iii2,f_timeout,
     &     time,amp_out,phase_out,angle_out,comp_out,ccp,
     &     layer, npoints_int,x_int,z_int,b_int,c_int,d_int,veltype, ptos,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,flag_smooth,
     &     arrival,narrivals,error,
     &     MAX_INT,MAXPOINTS_INT,MAX_RECEIVERS,MAX_ARRIVALS,MAX_RAYS,N_PARAM,N_RTPARAM,N_ARRPARAM,
     &     MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAX_INT,MAXPOINTS_INT,MAX_RAYS,N_PARAM,MAX_RECEIVERS,MAX_ARRIVALS,N_RTPARAM,N_ARRPARAM
      integer MAXP_XGRID,MAXP_ZGRID
      integer npoints_int(MAX_INT), iii2(MAXPOINTS_INT,MAX_INT),layer, veltype
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      real ptos(MAX_INT),rho1(MAX_INT),rho2(MAX_INT)
      integer nx_grid(MAX_INT),nz_grid(MAX_INT)
      real x_grid(MAXP_XGRID,MAX_INT),z_grid(MAXP_ZGRID,MAX_INT)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID,MAX_INT)
      real ray(N_PARAM,MAX_RAYS),tray(MAX_RAYS), amplitude(MAX_RAYS), phase(MAX_RAYS), angle(MAX_RAYS)
      integer nextray(0:MAX_RAYS),kmah(MAX_RAYS),fict(MAX_RAYS)
      real xz_rec(MAX_RECEIVERS),dtray
      integer xzindex,f_timeout,nrec,int_rec,lay_rec,flag_smooth

      real arrival(N_ARRPARAM,MAX_RECEIVERS,MAX_ARRIVALS)
      real time(MAX_RECEIVERS,MAX_ARRIVALS)
      real amp_out(MAX_RECEIVERS,MAX_ARRIVALS)
      real phase_out(MAX_RECEIVERS,MAX_ARRIVALS)
      real angle_out(MAX_RECEIVERS,MAX_ARRIVALS)
      integer comp_out
      real ccp(2,MAX_RECEIVERS)
      integer narrivals(MAX_RECEIVERS), error
      
      integer iray, irec, next, flag, iray1, iray2, narrtmp, updownflag

      iray = abs(nextray(0))
      next = abs(nextray(iray))

c      if (next .eq. 0) then
c         write(*,*) "No records available"
c      end if

c-----------------------------------
c Richtung des Strahls setzen
c UPDOWNFLAG = +1 : nach unten
c UPDOWNFLAG = -1 : nach oben
c
      if (layer .lt. int_rec) then
         updownflag = 1
      else
         updownflag = -1
      end if

c===================================================================================
c Ebene (0)
c Schleife ueber alle Strahlen, aufgezeichnet wird nur zwischen zwei eingefallenen Strahlen IRAY und NEXT
c
      do while (next .ne. 0)
c----------------------------
c Ebene (1)
c Nur Aufzeichnung, wenn Strahlen benachbart und KMAH gleich
c
         if (nextray(iray) .gt. 0 .and. kmah(iray) .eq. kmah(next)) then
c---------------
c Ebene (2)
c Strahlen sind zuerst rechts eingefallen
c
            if (ray(xzindex,iray) .lt. ray(xzindex,next)) then
c----------
c Ebene (3)
c Strahl IRAY und NEXT ausserhalb der Empfaengerlinie
c
               if (ray(xzindex,iray) .gt. xz_rec(nrec) .or. ray(xzindex,next) .lt. xz_rec(1)) then
                  iray = next
                  next = abs(nextray(iray))
c----------
c Ebene (3)
c Von rechts nach links werden nun die Strahlen durchlaufen.
c Sobald ein Empfaenger dazwischen ist, wird aufgezeichnet.
c Sobald KMAH-Index ungleich oder nicht mehr benachbart oder ausserhalb Empfaengerlinie
c  -> Sprung zurueck auf Ebene (1)
c
               else
                  irec = 1
                  do while (ray(xzindex,iray) .gt. xz_rec(irec))
                     irec = irec + 1
                  end do
                  do while (nextray(iray) .gt. 0 .and. kmah(iray) .eq. kmah(next) .and.
     &                 irec .le. nrec)
c-----------
c Normal: Geophon ist von Treffer links und rechts umschlossen
                     if (ray(xzindex,next) .gt. xz_rec(irec)) then
                        if ((ray(xzindex,iray)+ray(xzindex,next)) .lt. 2.*xz_rec(irec)) then
                           iray1 = next
                           iray2 = iray
                        else
                           iray1 = iray
                           iray2 = next
                        end if
                        if (narrivals(irec) .eq. MAX_ARRIVALS) goto 900
                        narrtmp = narrivals(irec) + 1
                        call get_record(ray(1,iray1),tray(iray1),amplitude(iray1),phase(iray1),angle(iray1),
     &                       ray(1,iray2),tray(iray2),amplitude(iray2),phase(iray2),angle(iray2),kmah(iray1),
     &                       xz_rec(irec),xzindex,int_rec,lay_rec,f_timeout,updownflag,fict(iray1),
     &                       layer,npoints_int,x_int,z_int,b_int,c_int,d_int,iii2,veltype,ptos,
     &                       nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &                       arrival(1,irec,narrtmp),time(irec,narrtmp),amp_out(irec,narrtmp),
     &                       phase_out(irec,narrtmp),angle_out(irec,narrtmp),comp_out,ccp(1,irec),flag,error,
     &                       MAX_INT,MAXPOINTS_INT,N_PARAM,N_RTPARAM,N_ARRPARAM,MAXP_XGRID,MAXP_ZGRID)
                        if (error .ne. 0) goto 999
                        if (flag .eq. 0) narrivals(irec) = narrtmp
                        irec = irec + 1
c-------
c Sonderfall: Geophon wurde exakt getroffen:
                     else if (ray(xzindex,next) .eq. xz_rec(irec)) then
                        if (nextray(next) .le. 0) then
                           if (narrivals(irec) .eq. MAX_ARRIVALS) goto 900
                           narrtmp = narrivals(irec) + 1
                           call get_record(ray(1,next),tray(next),amplitude(next),phase(next),angle(next),
     &                          ray(1,iray),tray(iray),amplitude(iray),phase(iray),angle(iray),kmah(next),
     &                          xz_rec(irec),xzindex,int_rec,lay_rec,f_timeout,updownflag,fict(iray1),
     &                          layer, npoints_int, x_int, z_int,b_int, c_int, d_int,iii2,veltype, ptos,
     &                          nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &                          arrival(1,irec,narrtmp),time(irec,narrtmp),amp_out(irec,narrtmp),
     &                          phase_out(irec,narrtmp),angle_out(irec,narrtmp),comp_out,ccp(1,irec),flag,error,
     &                          MAX_INT, MAXPOINTS_INT, N_PARAM, N_RTPARAM, N_ARRPARAM,MAXP_XGRID,MAXP_ZGRID)
                           if (error .ne. 0) goto 999
                           if (flag .eq. 0) narrivals(irec) = narrtmp
                           irec = irec + 1
                        end if
                        iray = next
                        next = abs(nextray(iray))
c-----
c Dieser Einsatz ist bereits berechnet (mehrere Treffer zwischen zwei Geophonen), naechster Strahl
                     else
                        iray = next
                        next = abs(nextray(iray))
                     end if
c Wieso geht das sonst nicht? (Abbruch wegen kmah(0)!)
c Zweite Abfrage ist noetig, da sich hier Strahlen schneiden koennen, ohne dass der KMAH-Index ungleich ist!
                     if (next .eq. 0 .or. (ray(1,iray) .gt. ray(1,next) .and. flag_smooth .ne. 0)) goto 1
                  end do
 1                if (irec .gt. nrec) then
                     iray = next
                     next = abs(nextray(iray))
                  end if
               end if
c
c Ende Ebene (3)
c---------------
c Ebene (2)
c Strahl IRAY und NEXT ausserhalb der Empfaengerlinie
c
            else if (ray(xzindex,next) .gt. xz_rec(nrec) .or. ray(xzindex,iray) .lt. xz_rec(1)) then
               iray = next
               next = abs(nextray(iray))
c---------------
c Ebene (2)
c Strahlen sind zuerst links eingefallen
c Von links nach rechts werden nun die Strahlen durchlaufen.
c Sobald ein Empfaenger dazwischen ist, wird aufgezeichnet.
c Sobald KMAH-Index ungleich oder nicht mehr benachbart oder ausserhalb Empfaengerlinie
c  -> Sprung zurueck auf Ebene (1)
c
            else
               irec = nrec
               do while (ray(xzindex,iray) .lt. xz_rec(irec))
                  irec = irec - 1
               end do
               do while (nextray(iray) .gt. 0 .and. (kmah(iray) .eq. kmah(next) .and.
     &              irec .ge. 1))
c-----------
c Normal: Geophon ist von Treffer links und rechts umschlossen
c                  write(*,*) " Treffer ",irec,narrivals(irec)
                  if (ray(xzindex,next) .lt. xz_rec(irec)) then
                     if ((ray(xzindex,iray)+ray(xzindex,next)) .gt. 2.*xz_rec(irec)) then
                        iray1 = next
                        iray2 = iray
                     else
                        iray1 = iray
                        iray2 = next
                     end if
                     if (narrivals(irec) .eq. MAX_ARRIVALS) goto 900
                     narrtmp = narrivals(irec) + 1
                     call get_record(ray(1,iray1),tray(iray1),amplitude(iray1),phase(iray1),angle(iray1),
     &                    ray(1,iray2),tray(iray2),amplitude(iray2),phase(iray2),angle(iray2),kmah(iray1),
     &                    xz_rec(irec),xzindex,int_rec,lay_rec,f_timeout,updownflag,fict(iray1),
     &                    layer, npoints_int, x_int, z_int,b_int, c_int, d_int,iii2,veltype, ptos,
     &                    nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &                    arrival(1,irec,narrtmp),time(irec,narrtmp),amp_out(irec,narrtmp),
     &                    phase_out(irec,narrtmp),angle_out(irec,narrtmp),comp_out,ccp(1,irec),flag,error,
     &                    MAX_INT, MAXPOINTS_INT, N_PARAM, N_RTPARAM, N_ARRPARAM,MAXP_XGRID,MAXP_ZGRID)
                     if (error .ne. 0) goto 999
                     if (flag .eq. 0) narrivals(irec) = narrtmp
                     irec = irec - 1
c-------
c Sonderfall: Geophon wurde exakt getroffen:
                  else if (ray(xzindex,next) .eq. xz_rec(irec)) then
                     if (nextray(next) .le. 0) then
                        if (narrivals(irec) .eq. MAX_ARRIVALS) goto 900
                        narrtmp = narrivals(irec) + 1
                        call get_record(ray(1,next),tray(next),amplitude(next),phase(next),angle(next),
     &                       ray(1,iray),tray(iray),amplitude(iray),phase(iray),angle(iray),kmah(next),
     &                       xz_rec(irec),xzindex,int_rec,lay_rec,f_timeout,updownflag,fict(iray1),
     &                       layer, npoints_int, x_int, z_int,b_int, c_int, d_int,iii2, veltype, ptos,
     &                       nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &                       arrival(1,irec,narrtmp),time(irec,narrtmp),amp_out(irec,narrtmp),
     &                       phase_out(irec,narrtmp),angle_out(irec,narrtmp),comp_out,ccp(1,irec),flag,error,
     &                       MAX_INT, MAXPOINTS_INT, N_PARAM, N_RTPARAM, N_ARRPARAM,MAXP_XGRID,MAXP_ZGRID)
                        if (error .ne. 0) goto 999
                        if (flag .eq. 0) narrivals(irec) = narrtmp
                        irec = irec - 1
                     end if
                     iray = next
                     next = abs(nextray(iray))
c-----
c Dieser Einsatz ist bereits berechnet (mehrere Treffer zwischen zwei Geophonen), naechster Strahl
                  else
                     iray = next
                     next = abs(nextray(iray))
                  end if
c Wieso geht das sonst nicht? (Abbruch wegen kmah(0)!)
c Zweite Abfrage ist noetig, da sich hier Strahlen schneiden koennen, ohne dass der KMAH-Index ungleich ist!
                  if (next .eq. 0 .or. (ray(1,iray) .lt. ray(1,next) .and. flag_smooth .ne. 0)) goto 2
               end do
 2             if (irec .lt. 1) then
                  iray = next
                  next = abs(nextray(iray))
               end if
            end if
c
c Ende Ebene (2)
c---------------------------------
c Ebene (1)
c Strahlen nicht benachbart oder KMAH-Index ungleich: naechster Strahl ist dran
c
         else
            iray = next
            next = abs(nextray(iray))
         end if
c
c Ende Ebene (1)
c----------------------------------------------------------------
      end do
c
c Ende Ebene (0)
c====================================================================================

      goto 999
 900  write(*,*) "----------------WARNING--------------------"
      write(*,*) "RECORD: Too many arrivals at receiver no.",irec
      write(*,*) "Enlarge parameter MAX_ARRIVALS in WFRONT"
      write(*,*) "-------------------------------------------"

 999  end

c****************************************************************************
c Unterprogramm GET_RECORD / Hauptprogramm WFRONT
c letzte Aenderung: 6.02.97  / nicht getestet
c Wird nur von RECORD aufgerufen
c RAY1 ist naeher an XZ_REC (receiver position) als RAY2
c FLAG =  0 : Alles in Butter
c FLAG = -1 : GET_RAYONINT hatte Probleme; kein Einsatz berechenbar
c
      subroutine get_record(ray1,tray1,amplitude1,phase1,angle1,ray2,tray2,amplitude2,phase2,angle2,kmah1,
     &     xz_rec,xzindex,int_rec,lay_rec,f_timeout,updownflag,fict,
     &     layer, npoints_int, x_int, z_int,b_int, c_int, d_int,iii2,veltype, ptos,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,dtray,
     &     arrival,time,amp_out,phase_out,angle_out,comp_out,ccp,flag,error,
     &     MAX_INT, MAXPOINTS_INT, N_PARAM, N_RTPARAM, N_ARRPARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer MAX_INT, MAXPOINTS_INT, N_PARAM, N_RTPARAM, N_ARRPARAM,MAXP_XGRID,MAXP_ZGRID
      real ray1(N_PARAM),tray1,amplitude1,phase1,angle1,ray2(N_PARAM),tray2,amplitude2,phase2,angle2,xz_rec,dtray
      integer kmah1,xzindex,int_rec,lay_rec, f_timeout, updownflag, fict
      integer npoints_int(MAX_INT), iii2(MAXPOINTS_INT,MAX_INT), layer, veltype
      real x_int(MAXPOINTS_INT,MAX_INT), z_int(MAXPOINTS_INT,MAX_INT)
      real b_int(MAXPOINTS_INT,MAX_INT), c_int(MAXPOINTS_INT,MAX_INT)
      real d_int(MAXPOINTS_INT,MAX_INT)
      real ptos(MAX_INT),rho1(MAX_INT),rho2(MAX_INT)
      integer nx_grid(MAX_INT),nz_grid(MAX_INT)
      real x_grid(MAXP_XGRID,MAX_INT),z_grid(MAXP_ZGRID,MAX_INT)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID,MAX_INT)

      real arrival(N_ARRPARAM)
      real time, amp_out, phase_out, angle_out, ccp(2)
      integer flag,error,comp_out

      integer NRTPARAM_LOCAL, actp, oldveltype, otherlayer
      parameter (NRTPARAM_LOCAL = 8)
      real raynew(8), rtray(NRTPARAM_LOCAL), dzdx, d2zdx2, ztmp
      real raynew0(8),amplnew, phasenew, ampltmp,phasetmp,pi,rho_rec,v_rec
      integer i, ficttmp, layer2, int2

      pi = 3.1415927
c ARRIVAL(1,irec,a)   : Laufzeit an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(2,irec,a)   : Amplitude(x) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(3,irec,a)   : Amplitude(z) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(4,irec,a)   : Phase(x) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(5,irec,a)   : Phase(z) an receiver 'irec', Anzahl der Einsaetze 'a'
c ARRIVAL(6,irec,a)   : Ausgangswinkel von der Quelle an receiver 'irec', Anzahl der Einsaetze 'a'

      rho_rec = 1.0

      raynew(xzindex) = xz_rec

c=============================================================================
c Setzen der zweiten Koordinate des Einsatzes... und layer2
c LAYER2 ist die Schicht, deren Groessen fuer die genaue Berechnung des Einsatzes eingehen soll.
c Im Normalfall die zuletzt durchlaufene Schicht, bei coinciding interfaces allerdings
c die Schicht, die zuletzt durchlaufen wurde  u n d  n i c h t  die Dicke '0' hatte
c----------------------------
c Grenzflaeche
c
      if (xzindex .eq. 1) then
         call get_z(raynew(1),x_int(1,int_rec),z_int(1,int_rec),d_int(1,int_rec),c_int(1,int_rec),
     &        b_int(1,int_rec),npoints_int(int_rec),
     &        raynew(2),actp,
     &        MAXPOINTS_INT)
c----------------------------
c Erst Setzen von LAYER2. Ist dieses Stueck ACTP coinciding, so muss
c eventuell das Geschwindigkeitsfeld einer anderen Schicht benutzt werden
c
c ACTP muss extra f"ur den Punkt des Strahls RAY1 (weil FICT vom Strahl 1 kommt) gefunden werden, da
c speziell im Fall FLAG_SMOOTH != 0 ueber einen corner point hinweg approximiert werden darf...
c
         call get_actp(ray1(1), int_rec, x_int, npoints_int,
     &        actp,
     &        MAX_INT, MAXPOINTS_INT)
         if (abs(fict) .gt. 999) then
c     wir befinden uns mitten innerhalb von einer unbekannten Anzahl an coinciding interfaces
c III2 am Rand wird gesucht... (III2 enthaelt die Gesamtanzahl zusammenfallender Grenzflaechen)
            if (abs(fict) .gt. 1999) then
               int2 = int_rec + updownflag * (abs(fict/1000) - 1)
               call get_actp(raynew(1),int2,x_int,npoints_int,
     &              actp,
     &              MAX_INT,MAXPOINTS_INT)
               if (updownflag .gt. 0) then
                  layer2 = int2 - abs(iii2(actp,int2))
               else
                  layer2 = int2 + abs(iii2(actp,int2)) - 1
               end if
c Wir befinden uns am letzten coinciding interface.
c Hier koennen wir direkt III2 erfahren
            else
               layer2 = layer + updownflag * (1 - abs(iii2(actp,int_rec)))
            end if
c Normalfall:
         else
            layer2 = layer
         end if
c-----------------------------
c Bohrloch
c
      else
         call get_z(raynew(2),x_int(1,int_rec),z_int(1,int_rec),d_int(1,int_rec),c_int(1,int_rec),
     &        b_int(1,int_rec),npoints_int(int_rec),
     &        raynew(1),actp,
     &        MAXPOINTS_INT)
         layer2 = layer
      end if
c====================================================================
c Slowness, Laufzeit, Spreading, Amplitude, Phase, Angle am Receiverpunkt:
c

      call get_rayonint(ray1,amplitude1,phase1,tray1,angle1,ray2,amplitude2,phase2,tray2,angle2,
     &     arrival(1),veltype,ptos(layer2),layer2,dtray,
     &     nx_grid(layer2),nz_grid(layer2),x_grid(1,layer2),z_grid(1,layer2),
     &     v_grid(1,1,1,layer2),
     &     raynew,amplnew,phasenew,arrival(6),flag,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      if (flag .ne. 0) goto 999

c==========================================================================
c Fall: Receiver liegen in der aktuellen Schicht
c
      if (lay_rec .eq. layer) then
         flag = 0

         layer2 = lay_rec
c==========================================================================
c Fall: Receiver liegen in der naechsten Schicht,
c       Strahlen muessen erst transmittiert werden:
      else
c-----------------------------
c in gewissen Faellen muss noch transmittiert werden...
c Faelle : normal(0) oder coinciding interfaces-Paket gerade getroffen
         if (fict .eq. 0 .or. mod(fict,1000) .ne. 0) then
            if (fict .gt. 999) then
               ficttmp = fict / 1000
               layer2 = layer + updownflag * ficttmp
            else
               layer2 = lay_rec
            end if
            if (NRTPARAM_LOCAL .ne. N_RTPARAM) goto 900
            call get_z2(raynew(1),x_int(1,int_rec),z_int(1,int_rec),d_int(1,int_rec),c_int(1,int_rec),
     &           b_int(1,int_rec),npoints_int(int_rec),
     &           ztmp,dzdx,d2zdx2,actp,
     &           MAXPOINTS_INT)
            call get_rtsect(raynew,dzdx,d2zdx2,
     &           nx_grid(layer),nz_grid(layer),x_grid(1,layer),z_grid(1,layer),
     &           v_grid(1,1,1,layer),ptos(layer),veltype,layer,
     &           rtray,
     &           N_RTPARAM,MAXP_XGRID,MAXP_ZGRID)
            oldveltype = veltype
            otherlayer = layer
            call rtcalc(raynew,rtray,amplnew,phasenew,layer,layer2,
     &           otherlayer,oldveltype,veltype,int_rec,1,ptos,0,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,rho1,rho2,
     &           raynew0,ampltmp,phasetmp,flag,
     &           MAX_INT,N_RTPARAM,N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            if (flag .ne. 0) goto 999

            do i = 1, N_PARAM
               raynew(i) = raynew0(i)
            end do
            phasenew = phasetmp
            amplnew = ampltmp
c---------------------
c ...in anderen laeuft der Strahl direkt durch
c Faelle fictious interface/s(< 0) oder innerhalb coinciding interfaces (> 999)
         else
            ficttmp = abs(fict / 1000)
            layer2 = lay_rec + updownflag * ficttmp
            oldveltype = veltype
            otherlayer = layer
         end if
      end if
c==========================================================================


c Sicherheit: In der Naehe der Quelle sowie caustics...
      if (abs(raynew(5)) .lt. 0.000001) raynew(5) = 0.000001

c Geschwindigkeit am receiver:
      call velonly(raynew(1),raynew(2), nx_grid(layer2), nz_grid(layer2),
     &     x_grid(1,layer2), z_grid(1,layer2), v_grid(1,1,1,layer2),
     &     veltype, ptos(layer2),layer2,
     &     v_rec,
     &     MAXP_XGRID, MAXP_ZGRID)
      
c Dichte am receiver:
      if (veltype .gt. 0) then
         rho_rec = rho1(layer2) + rho2(layer2)*v_rec
      else
         rho_rec = rho1(layer2) + rho2(layer2)*v_rec*ptos(layer2)
      end if

c Amplitude = Amplitude0 * sqrt( 1.0 / | Q * V_REC * RHO_REC | )
c dann: AX = Amplitude * (PX * V_REC)
c       AZ = Amplitude * (PZ * V_REC)
c daher direkt: Amplitude = Amplitude0 * sqrt( 1.0 / | Q * RHO_REC * V_REC | ) * V_REC

      ampltmp = amplnew * sqrt( v_rec / (abs(raynew(5)) * rho_rec))

      phasetmp = phasenew - 0.5 * pi * kmah1
c---------------------- Amplitude ---
c p-wave (V_REC steckt bereits in AMPLTMP drin)
      if (veltype .gt. 0) then
         arrival(2) = ampltmp * raynew(3)
         arrival(3) = ampltmp * raynew(4)
c s-wave
      else
         arrival(2) = ampltmp * raynew(4)
         arrival(3) = -ampltmp * raynew(3)
      end if
c---------------------- Phase --------
c
      if (arrival(2) .lt. 0.0) then
         arrival(2) = -arrival(2)
         arrival(4) = phasetmp - pi
      else
         arrival(4) = phasetmp
      end if
      if (arrival(3) .lt. 0.0) then
         arrival(3) = -arrival(3)
         arrival(5) = phasetmp - pi
      else
         arrival(5) = phasetmp
      end if
c
c Ende Setzen von Slowness, Amplitude und Phase
c-------------------------------------------------------------------------
c Setzen der Phase auf Wert zwischen -PI und PI
c
      do i = 4, 5
         do while (arrival(i) .gt. pi)
            arrival(i) = arrival(i) - 2.*pi
         end do
         do while (arrival(i) .lt. (-pi))
            arrival(i) = arrival(i) + 2.*pi
         end do
      end do

      ccp(1) = raynew(7)
      ccp(2) = raynew(8)
      time = arrival(1)
      amp_out = arrival(comp_out)
      phase_out = arrival(comp_out+2)
      angle_out = arrival(6)

      if (f_timeout .ne. 0) then
c 7&8 not fully initialized here:
         write(f_timeout,*,err=999) xz_rec,arrival(1),raynew(7),raynew(8),arrival(2),arrival(3)
      end if

 1000 format(1x, f10.5, f12.8)
 1001 format(1x, f10.5, 4f12.8)
      goto 999
 900  write(*,*) "GET_RECORD: Parameter NRTPARAM_LOCAL not equal N_RTPARAM from WFRONT!"
      write(*,*) "Change NRTPARAM_LOCAL in GET_RECORD to ",N_RTPARAM
      error = 21

 999  end

c*************************************************************************************************
c Unterprogramm GET_RAYONINT / Hauptprogramm WFRONT
c letzte Aenderung 6.02.97 und getestet / Es fehlt: Ueberpruefung, ob neuer Strahl in die
c richtige Richtung verlaeuft (koennte ganz eventuell auftreten)
c Wird aufgerufen von: GET_RECORD, GET_CRITRAY, SORTSECT
c Aufruf der Unterprogramme PARAXIAL, RUNKUTTA, VELOCITY
c Berechnet alle neuen Strahl-Parameter eines bekannten Punktes (RAYNEW(1),RAYNEW(2)).
c Der Punkt liegt auf einer Grenzflaeche, mit Nachbarstrahlen RAY1 und RAY2.
c Die Slowness wird (ungefaehr entlang einer Wellenfront) paraxial bestimmt.
c Alle anderen Parameter werden entlang einer Wellenfront linear interpoliert!
c
      subroutine get_rayonint(ray1,amplitude1,phase1,tray1,angle1,ray2,amplitude2,phase2,tray2,angle2,
     &     tnew,veltype,ptos,actlayer,dtray,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     raynew,amplnew,phasenew,anglenew,flag,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

      implicit none

      integer N_PARAM, MAXP_XGRID,MAXP_ZGRID
      real ray1(N_PARAM),tray1,amplitude1,phase1,angle1,ray2(N_PARAM),tray2,amplitude2,phase2,angle2
      integer veltype, actlayer
      real tnew, ptos, dtray
      integer nx_grid, nz_grid, flag
      real x_grid(MAXP_XGRID), z_grid(MAXP_ZGRID)
      real v_grid(4,MAXP_ZGRID,MAXP_XGRID)

      real raynew(N_PARAM),amplnew,phasenew,anglenew

      real rayleft(N_PARAM), rayleft0(N_PARAM)
      real rayright(N_PARAM), rayright0(N_PARAM)
      real vel1(6), vel2(6), pxtmp, pztmp, ttmp, tmid, dtnew, dtnewhalf, dtnew0
      real totaldist, dist1, dist10, dist2, distall, part
      integer i, nloops, iloop

      if (abs(ray1(5)) .lt. 0.00001 .or. abs(ray2(5)) .lt. 0.00001) then
         write(*,*) 'GET_RAYONINT: CAUSTIC! New ray could not be detected.'
         flag = -1
         goto 999
      end if

      flag = 0

c======================================================================
c 1. Verhaeltnis der Abstaende der Strahlen zum neuen Punkt feststellen
c
      totaldist = (ray2(1) - ray1(1))**2 + (ray2(2) - ray1(2))**2
      dist1 = (ray1(1) - raynew(1))**2 + (ray1(2) - raynew(2))**2
      if (totaldist .lt. 0.000001) then
         part = 0.5
      else
         part = sqrt(dist1 / totaldist)
      end if
c======================================================================
c 2. Laufzeit mit gleichem Verhaeltnis mitteln.
c    Dies ist die vermutete neue Laufzeit.
c    ...sind die Zeiten gleich, kann man sich den Anfang ersparen:
c
      tmid = tray1 + part * (tray2 - tray1)
c======================================================================
c 3. Beide Strahlen bis zur Laufzeit TMID propagieren.
c    Sinn: Dann ist der neue Punkt in einem abschaetzbaren Abstand, naemlich hoechstens
c          SPREAD_MAX von den Strahlen entfernt -> paraxiale Approximation gueltig
c
      do i = 1, N_PARAM
         rayleft(i) = ray1(i)
         rayright(i) = ray2(i)
      end do
c----------------------------------
c Strahl 1 auf RAYLEFT0 propagieren:
c
c Zeitschritt DTRAY darf nicht ueberschritten werden, deshalb diese muehsame Prozedur:
      dtnew = tmid - tray1
      if (abs(dtnew) .gt. dtray) then
         nloops = abs(int(dtnew/dtray))
         dtnew = mod(dtnew,dtray)

         dtnew0 = sign(1.0,dtnew) * dtray
         dtnewhalf = 0.5 * dtnew0
         do iloop = 1, nloops
            call runkutta(rayleft,dtnew0,dtnewhalf,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           rayleft0,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            do i = 1, N_PARAM
               rayleft(i) = rayleft0(i)
            end do
         end do
      end if

      dtnewhalf = 0.5*dtnew
      call runkutta(rayleft,dtnew,dtnewhalf,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     rayleft0,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

c-----------------------------------
c Strahl 2 auf RAYRIGHT0 propagieren:
c
      dtnew = tmid - tray2
      if (abs(dtnew) .gt. dtray) then
         nloops = abs(int(dtnew/dtray))
         dtnew = mod(dtnew,dtray)

         dtnew0 = sign(1.0,dtnew) * dtray
         dtnewhalf = 0.5*dtnew0
         do iloop = 1, nloops
            call runkutta(rayright,dtnew0,dtnewhalf,
     &           nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &           ptos,veltype,actlayer,
     &           rayright0,
     &           N_PARAM,MAXP_XGRID,MAXP_ZGRID)
            do i = 1, N_PARAM
               rayright(i) = rayright0(i)
            end do
         end do
      end if

      dtnewhalf = 0.5*dtnew
      call runkutta(rayright,dtnew,dtnewhalf,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     rayright0,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

c======================================================================
c 4. (quadr.) Verhaeltnis der neuen Abstaende zur quadratischen Mittelung der beiden paraxialen
c    Approximationen errechnen (evt totaldist und dist1 nicht noch einmal extra ausrechnen)
c    "quadratische" Mittelung:
c
      totaldist = (rayleft0(1) - rayright0(1))**2 + (rayleft0(2) - rayright0(2))**2
      dist1 = (rayleft0(1) - raynew(1))**2 + (rayleft0(2) - raynew(2))**2
      dist10 = sqrt(dist1)
      dist2 = sqrt(totaldist) - dist10
      distall = totaldist - 2.*dist10*dist2
      if (distall .lt. 0.000001) then
         part = 0.5
      else
         part = dist1/distall
      end if
c======================================================================
c 5. paraxiale Approximation der neuen Slowness sowie der Laufzeit
c    von beiden Strahlen aus, danach quadratische Mittelung
c
      call velocity(rayleft0(1), rayleft0(2),nx_grid,nz_grid,
     &     x_grid,z_grid,v_grid,veltype,ptos,actlayer,
     &     vel1,
     &     MAXP_XGRID, MAXP_ZGRID)
      call velocity(rayright0(1), rayright0(2), nx_grid, nz_grid,
     &     x_grid,z_grid,v_grid,veltype,ptos,actlayer,
     &     vel2,
     &     MAXP_XGRID, MAXP_ZGRID)

      call paraxial(rayleft0,tmid,vel1,
     &     raynew,ttmp,N_PARAM)
      pxtmp = raynew(3)
      pztmp = raynew(4)
      call paraxial(rayright0,tmid,vel2,
     &     raynew,tnew,N_PARAM)
      raynew(3) = pxtmp + part * (raynew(3) - pxtmp)
      raynew(4) = pztmp + part * (raynew(4) - pztmp)
      tnew = ttmp + part * (tnew - ttmp)

c---------------------------
c KEIN neuer Strahl, wenn RAYNEW in die falsche Richtung laeuft!
c eventuell eine Ueberpruefung starten, etwa wie:
c  (dzdx bereits importieren)
c
c      chkbound = raynew(3)*dzdx - raynew(4)
c      if (chkbound*float(cross) .gt. 0.0) then
c         flag = -1
c         goto 999
c      end if

c======================================================================
c 6. Propagation der Elternstrahlen zur exakten Laufzeit TNEW
c    Sinn: Bestimmung der restlichen neuen Parameter.
c    Diese Extrapropagation ist nur fuer Q und P notwendig (und fast unnoetig).
c---------------------------------------------
c Strahl 1 (RAYLEFT0) auf RAYLEFT propagieren:
c
      dtnew = tnew - tmid
      dtnewhalf = 0.5*dtnew

      call runkutta(rayleft0,dtnew,dtnewhalf,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     rayleft,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

c-----------------------------------------------
c Strahl 2 (RAYRIGHT) auf RAYRIGHT0 propagieren:

      call runkutta(rayright0,dtnew,dtnewhalf,
     &     nx_grid,nz_grid,x_grid,z_grid,v_grid,
     &     ptos,veltype,actlayer,
     &     rayright,
     &     N_PARAM,MAXP_XGRID,MAXP_ZGRID)

c======================================================================
c 7. Bestimmung des neuen (linearen) Verhaeltnisses der Abstaende, oder einfach altes nehmen
c
c      totaldist = (rayleft(1)-rayright(1))**2 + (rayleft(2)-rayright(2))**2
c      dist1 = (rayleft(1)-raynew(1))**2 + (rayleft(2)-raynew(2))**2

      if (totaldist .lt. 0.000001) then
         part = 0.5
      else
         part = sqrt(dist1 / totaldist)
      end if

c======================================================================
c 8. Bestimmung der restlichen Strahlparameter durch lineare Interpolation
c
      do i = 5, 8
         raynew(i) = rayleft(i) + part*(rayright(i) - rayleft(i))
      end do
      amplnew = amplitude1 + part*(amplitude2 - amplitude1)
      phasenew = phase1 + part*(phase2 - phase1)
      anglenew = angle1 + part*(angle2 - angle1)

 999  end

