c Copyright (c) Colorado School of Mines, 2013.
c All rights reserved.

c Unterprogramm SORTRAYCODE / Hauptprogramm WFRONT
c Aufruf des Unterprogramms CODE_GEN1, GET_Z, GET_Z1
c
c Sortiert die eingegebenen Strahlen-Codes und setzt das Feld HOLD(,) folgendermassen:
c - erster Code ist der kuerzeste
c - Codes mit gleichem Anfang stehen hintereinander
c - Reihenfolge bei ungleichem Anfang: CODE(nstep,ncode) <= CODE(nstep,ncode+1)...
c - Strahlencodestuecke, die gespeichert werden sollen, erhalten ein HOLD(nstep,ncode) = 1
c   d.h.: der entsprechende HALT, der die Schnittpunkte der Strahlen enthaelt, die den
c   Strahlencode 'ncode' durchlaufen haben, soll erhalten bleiben
c
      subroutine sortraycode( codesteps, ncodes, codestep_plot, f_out,
     &     int_source, lay_source, code, hold, flag_dup, error,
     &     MAX_CODESTEPS)

      implicit none

      integer MAX_CODESTEPS
      integer ncodes, int_source, lay_source, error, f_out, flag_dup
      integer code(0:MAX_CODESTEPS,ncodes),codesteps(ncodes),hold(0:MAX_CODESTEPS,ncodes)
      integer codestep_plot(ncodes)

      integer CODESTEPS_LOCAL
      parameter(CODESTEPS_LOCAL = 100)
      integer ncode1(CODESTEPS_LOCAL),ncode2(CODESTEPS_LOCAL),actstep,act_recursion
      integer codetmp(0:CODESTEPS_LOCAL), codestepstmp, codeplottmp
      integer icode, jcode

      integer i, j
      integer minus_flags, plus_flags

c-------------------------
c Initialisierung von HOLD
c
      do i = 1, ncodes
         do j = 0, MAX_CODESTEPS
            hold(j,i) = 0
         end do
      end do
c===================================================================================
c Pruefung auf Konsistenz
c Danach evt. Vorsortierung nach code(0,...) (int_source .eq. 0)
c
c-------------------------------
c Fall: source innerhalb Schicht
c
      if (int_source .eq. 0) then

c 2. Sortieren nach Up/Down flag code(0,...)

c
c Code flag 0 funktioniert noch nicht...
c
         plus_flags = 0
         minus_flags = 0
         do i = 1, ncodes
            if( code(0,i) .gt. 0 ) then
               plus_flags = plus_flags + 1
            else if( code(0,i) .lt. 0 ) then
               minus_flags = minus_flags + 1
            else
               
            endif
         enddo
         if( minus_flags .gt. 0 .and. plus_flags .gt. 0 ) then
            do icode = 1, minus_flags
               if( code(0,icode) .gt. 0 ) then
                  do jcode = icode, ncodes
                     if( code(0,jcode) .lt. 0 ) then
c     Vertausche die beiden codes:
                        if (codesteps(jcode) .gt. CODESTEPS_LOCAL) goto 960
                        do i = 0, codesteps(jcode)
                           codetmp(i) = code(i,jcode)
                        end do
                        codestepstmp = codesteps(jcode)
                        codeplottmp = codestep_plot(jcode)
                        
                        do i = 0, codesteps(icode)
                           code(i,jcode) = code(i,icode)
                        end do
                        do i = codesteps(icode)+1, codesteps(jcode)
                           code(i,jcode) = 0
                        end do
                        codesteps(jcode) = codesteps(icode)
                        codestep_plot(jcode) = codestep_plot(icode)
                        
                        do i = 0, codestepstmp
                           code(i,icode) = codetmp(i)
                        end do
                        do i = codestepstmp+1, codesteps(icode)
                           code(i,icode) = 0
                        end do
                        codesteps(icode) = codestepstmp
                        codestep_plot(icode) = codeplottmp

                     endif
                  enddo
               endif
            enddo
            icode = minus_flags
         else
            icode = ncodes
         endif

c Setzen von hold(0,..) im Falle, dass derselbe Anfangscodeschritt fuer oben und unten vorliegt:
         if (icode .lt. ncodes) then
            do 10 i = 1, icode
               do 20 j = icode+1, ncodes
                  if (code(1,j) .eq. code(1,i)) then
                     hold(0,i) = j
                     goto 10
                  end if
 20            continue
 10         continue
         end if
c------------------------------
c Fall: source auf Grenzflaeche
      else
c Loeschen des auf/ab-Werts code(0,.), da er hier nicht gebraucht wird
         do i = 1, ncodes
            code(0,i) = 0
         end do
         icode = ncodes
      end if

c===========================================================================================
c Sortieren der Codes
c
c Fall up/down-flag = 0, -1 (oder +1, falls kein -1 vorhanden)
      ncode1(1) = 1
      ncode2(1) = icode
      actstep = 1
      act_recursion = 1
      call code_gen1(code,codesteps,codestep_plot,ncodes,codetmp(1),
     &     ncode1,ncode2,actstep,act_recursion,hold,error,
     &     MAX_CODESTEPS,CODESTEPS_LOCAL)

      if (error .ne. 0) goto 940

c Fall up/down-flag = +1 (falls vorhanden)
      ncode1(1) = icode + 1
      ncode2(1) = ncodes
      actstep = 1
      act_recursion = 1
      call code_gen1(code,codesteps,codestep_plot,ncodes,codetmp(1),
     &     ncode1,ncode2,actstep,act_recursion,hold,error,
     &     MAX_CODESTEPS,CODESTEPS_LOCAL)

      if (error .ne. 0) goto 940


c-------------------------
c Ausgabe der sortierten Codes in output file
c
c      write(f_out,*) "sorted raycodes:"
c      do i = 1, ncodes
c         write(f_out,1000) code(0,i),codesteps(i),(code(j,i),j = 1, codesteps(i))
c      end do

c Ausgabe des Feldes HOLD(,)
c      write(f_out,*) "hold:"
c      do i = 1, ncodes
c         write(f_out,1005) (hold(j,i),j = 0, codesteps(i))
c      end do

c-------------------------
c Allow duplicate ray codes
      if ( flag_dup .eq. 1 ) goto 999
c Do NOT allow duplicate ray codes
      do i = 1, ncodes - 1

c Pruefe auf genau gleiche Codes, oder, ob ein Code ganz in einem anderen vorhanden ist:
         j = 0
         do while (j .le. codesteps(i))
            if (code(j,i) .ne. code(j,i+1)) goto 13
            j = j + 1
         end do
 13      if (j .gt. codesteps(i)) goto 950
      end do
      
 1000 format(60i3)

      goto 999
 940  write(*,*) "SORTRAYCODE: Too many recursions."
      write(*,*) "Enlarge the PARAMETER 'CODESTEPS_LOCAL' in SORTCODE.F"
      error = 2
      goto 999
 950  write(*,*) "SORTRAYCODE: This code is fully contained in another one"
      write(*,*) " Delete this code:"
      write(*,1000) (code(j,i),j = 1, codesteps(i))
      error = 3
      goto 999
 960  write(*,*) "SORRAYTCODE: Too many codesteps!"
      write(*,*) "          Enlarge PARAMETER 'MAX_CODESTEPS_LOCAL' in SORTCODE.F"
      error = 4

 999  end

c************************************************************************
c CODE_GEN1 und CODE_GEN2 sind fast exakt die gleichen Unterprogramme. Ueber
c sie wird eine Rekursion simuliert, die in Fortran eigentlich nicht moeglich ist.
c Der rekursive Aufruf geschieht immer an das jeweils andere Unterprogramm.
c Bei jedem rekursiven Aufruf wird die Variable ACT_RECURSION um eins erhoeht; ACT_RECURSION
c gibt also die Tiefe der Rekursion an.
c CODE_GEN* ordnet die Codes NCODE1 bis NCODE2, und zwar lediglich den Schritt ACTSTEP.
c
      subroutine code_gen1(code,codesteps,codestep_plot,ncodes,codetmp,
     &     ncode1,ncode2,actstep,act_recursion,hold,error,
     &     MAX_CODESTEPS,CODESTEPS_LOCAL)
      implicit none

      integer MAX_CODESTEPS,CODESTEPS_LOCAL, ncodes
      integer code(0:MAX_CODESTEPS,ncodes),codesteps(ncodes),hold(0:MAX_CODESTEPS,ncodes)
      integer codestep_plot(ncodes)
      integer ncode1(CODESTEPS_LOCAL),ncode2(CODESTEPS_LOCAL),actstep,act_recursion,error
      integer codetmp(CODESTEPS_LOCAL),codestepstmp,codeplottmp,actcode,holdtmp
      integer i,j,i1,k

      if (act_recursion .eq. CODESTEPS_LOCAL) error = 12
      if (error .ne. 0) goto 999

      do 1 i = ncode1(act_recursion)+1, ncode2(act_recursion)
         do 2 j = ncode1(act_recursion), i-1
c----------------------
c Einsetzen von code 'I' in code 'J', sowie Aufruecken aller codes dazwischen
c
c ...wenn : ('Codeschritt' ist immer der aktuelle)
c    Code 'J' darf nicht zu Ende sein!
c 1. Codeschritt von 'I' ist 'kleiner' als Codeschritt von 'J'  .AND.
c    Code 'J' ist nicht zu Ende  .AND.
c    ( es ist nicht der erste Schritt  .AND.  vorheriger Codeschritt von 'J' ist nicht gleich diesem )
c 2. es ist der letzte Schritt von 'I'
c 3. 'I' > 1 (ist wahrscheinlich ueberfluessig)  .AND.
c    vorheriger Codeschritt von 'I' ist gleich diesem  .AND.
c
            if (actstep .le. codesteps(j)) then
               if ((abs(code(actstep,i)) .lt. abs(code(actstep,j)) .and. code(actstep,j) .ne. 0 .and.
     &              (actstep .gt. 1 .and. code(actstep,j) .ne. code(actstep-1,j))) .or. 
     &              code(actstep,i) .eq. 0 .or. (i .gt. 1 .and. code(actstep,i) .eq. code(actstep-1,i))) then
                  do k = actstep, codesteps(i)
                     codetmp(k) = code(k,i)
                  end do
                  codestepstmp = codesteps(i)
                  codeplottmp = codestep_plot(i)
                  holdtmp = hold(0,i)
                  do i1 = i-1, j, -1
                     do k = actstep, codesteps(i1)
                        code(k,i1+1) = code(k,i1)
                     end do
                     codesteps(i1+1) = codesteps(i1)
                     codestep_plot(i1+1) = codestep_plot(i1)
                     hold(0,i1+1) = hold(0,i1)
                     do k = codesteps(i1)+1, MAX_CODESTEPS
                        code(k,i1+1) = 0
                     end do
                  end do

                  codesteps(j) = codestepstmp
                  codestep_plot(j) = codeplottmp
                  hold(0,j) = holdtmp
                  do k = actstep, codestepstmp
                     code(k,j) = codetmp(k)
                  end do
                  do k = codestepstmp+1, MAX_CODESTEPS
                     code(k,j) = 0
                  end do
                  goto 1
               end if
            end if
c----------
 2       continue
 1    continue

      actcode = ncode1(act_recursion)
      if (actcode .gt. ncodes ) goto 999
      if (code(actstep,actcode) .eq. 0) then
         actcode = actcode + 1
      end if
      if (actcode .lt. ncode2(act_recursion)) then
         do while(actcode .lt. ncode2(act_recursion))
            do i = actcode+1, ncode2(act_recursion)
               if (code(actstep,i) .ne. code(actstep,actcode)) then
                  ncode2(act_recursion+1) = i-1
                  goto 3
               end if
            end do
            ncode2(act_recursion+1) = ncode2(act_recursion)
 3          continue
            if (actcode .ne. ncode2(act_recursion+1)) then
               do i = actcode, ncode2(act_recursion+1)-1
                  hold(actstep,i) = 1
                  if (actstep .gt. 1) then
                     if (ncode2(act_recursion) .eq. ncode2(act_recursion+1)) then
                        hold(actstep-1,i) = 0
                     end if
                  end if
               end do
               act_recursion = act_recursion + 1
               actstep = actstep + 1
               ncode1(act_recursion) = actcode
               call code_gen2(code,codesteps,codestep_plot,ncodes,codetmp,
     &              ncode1,ncode2,actstep,act_recursion,hold,error,
     &              MAX_CODESTEPS,CODESTEPS_LOCAL)
               act_recursion = act_recursion - 1
               actstep = actstep - 1
            end if
            actcode = ncode2(act_recursion+1) + 1
         end do
      end if

 999  end

c***********************************************************************************

      subroutine code_gen2(code,codesteps,codestep_plot,ncodes,codetmp,
     &     ncode1,ncode2,actstep,act_recursion,hold,error,
     &     MAX_CODESTEPS,CODESTEPS_LOCAL)
      implicit none

      integer MAX_CODESTEPS, CODESTEPS_LOCAL, ncodes
      integer code(0:MAX_CODESTEPS,ncodes),codesteps(ncodes),hold(0:MAX_CODESTEPS,ncodes)
      integer codestep_plot(ncodes)
      integer ncode1(MAX_CODESTEPS),ncode2(MAX_CODESTEPS),actstep,act_recursion,error
      integer codetmp(CODESTEPS_LOCAL),codestepstmp,codeplottmp,actcode,holdtmp
      integer i,j,i1,k

      if (act_recursion .eq. CODESTEPS_LOCAL) error = 12
      if (error .ne. 0) goto 999
      do 1 i = ncode1(act_recursion)+1, ncode2(act_recursion)
         do 2 j = ncode1(act_recursion), i-1
c----------------------
c Einsetzen von code 'I' in code 'J', sowie Aufruecken aller codes dazwischen
c
c ...wenn : ('Codeschritt' ist immer der aktuelle)
c    Code 'J' darf nicht zu Ende sein!
c 1. Codeschritt von 'I' ist 'kleiner' als Codeschritt von 'J'  .AND.
c    Code 'J' ist nicht zu Ende  .AND.
c    ( es ist nicht der erste Schritt  .AND.  vorheriger Codeschritt von 'J' ist nicht gleich diesem )
c 2. es ist der letzte Schritt von 'I'
c 3. 'I' > 1 (ist wahrscheinlich ueberfluessig)  .AND.
c    vorheriger Codeschritt von 'I' ist gleich diesem  .AND.
c
            if (actstep .le. codesteps(j)) then
               if ((abs(code(actstep,i)) .lt. abs(code(actstep,j)) .and. code(actstep,j) .ne. 0 .and.
     &              (actstep .gt. 1 .and. code(actstep,j) .ne. code(actstep-1,j))) .or. 
     &              code(actstep,i) .eq. 0 .or. (i .gt. 1 .and. code(actstep,i) .eq. code(actstep-1,i))) then
                  do k = actstep, codesteps(i)
                     codetmp(k) = code(k,i)
                  end do
                  codestepstmp = codesteps(i)
                  codeplottmp = codestep_plot(i)
                  holdtmp = hold(0,i)
                  do i1 = i-1, j, -1
                     do k = actstep, codesteps(i1)
                        code(k,i1+1) = code(k,i1)
                     end do
                     codesteps(i1+1) = codesteps(i1)
                     codestep_plot(i1+1) = codestep_plot(i1)
                     hold(0,i1+1) = hold(0,i1)
                     do k = codesteps(i1)+1, MAX_CODESTEPS
                        code(k,i1+1) = 0
                     end do
                  end do
                  
                  codesteps(j) = codestepstmp
                  codestep_plot(j) = codeplottmp
                  hold(0,j) = holdtmp
                  do k = actstep, codestepstmp
                     code(k,j) = codetmp(k)
                  end do
                  do k = codestepstmp+1, MAX_CODESTEPS
                     code(k,j) = 0
                  end do
                  goto 1
               end if
            end if
c----------
 2       continue
 1    continue

      actcode = ncode1(act_recursion)

      if (actcode .gt. ncodes ) goto 999
      if (code(actstep,actcode) .eq. 0) then
         actcode = actcode + 1
      end if
      if (actcode .lt. ncode2(act_recursion)) then
         do while(actcode .lt. ncode2(act_recursion))
            do i = actcode+1, ncode2(act_recursion)
               if (code(actstep,i) .ne. code(actstep,actcode)) then
                  ncode2(act_recursion+1) = i-1
                  goto 3
               end if
            end do
            ncode2(act_recursion+1) = ncode2(act_recursion)
 3          continue
            if (actcode .ne. ncode2(act_recursion+1)) then
               do i = actcode, ncode2(act_recursion+1)-1
                  hold(actstep,i) = 1
                  if (actstep .gt. 1) then
                     if (ncode2(act_recursion) .eq. ncode2(act_recursion+1)) then
                        hold(actstep-1,i) = 0
                     end if
                  end if
               end do
               act_recursion = act_recursion + 1
               actstep = actstep + 1
               ncode1(act_recursion) = actcode
               call code_gen1(code,codesteps,codestep_plot,ncodes,codetmp,
     &              ncode1,ncode2,actstep,act_recursion,hold,error,
     &              MAX_CODESTEPS,CODESTEPS_LOCAL)
               act_recursion = act_recursion - 1
               actstep = actstep - 1
            end if
            actcode = ncode2(act_recursion+1) + 1
         end do
      end if
      
 999  end

