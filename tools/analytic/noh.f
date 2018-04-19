      program noh_test
      implicit none

c..tests the noh solver
c..taken from http://cococubed.asu.edu/research_pages/noh.shtml

c..declare
      character*80     outfile,string
      integer          i,nstep,iargc
      double precision time,zpos,
     1                 rho0,vel0,gamma,xgeom,
     2                 den,ener,pres,vel,
     3                 zlo,zhi,zstep,value


c..popular formats
 01   format(1x,t4,a,t8,a,t22,a,t36,a,t50,a,t64,a,t78,a,t92,a)
 02   format(1x,i4,1p8e12.4)



c..if your compiler handles command line arguments
c..get the input arguments from the command line
c      i = iargc()
c      if (i. lt. 2) stop 'too few arguments'
c      call getarg(1,string)
c      nstep = int(value(string))
c      call getarg(2,outfile)


c..input parameters in cgs
      time   =  0.2
      rho0   =  1
      vel0   = -0.1d0
      gamma  = 5.d0/3.d0
      xgeom  = 2.0d0


c..number of grid points, spatial domain, spatial step size
      nstep = 1000
      zlo   = 0.0d0
      zhi   = 0.5d0
      zstep = (zhi - zlo)/float(nstep)


c..output file 
      outfile = 'noh.dat'
      open(unit=2,file=outfile,status='unknown')
      write(2,02) nstep,time
c      write(6,01) 'i','x','den','ener','pres','vel'
      write(2,01) 'i','x','den','ener','pres','vel'


c..to match rage output, use the mid-cell points
      do i=1,nstep
       zpos   = zlo + 0.5d0*zstep + float(i-1)*zstep

       call noh_1d(time,zpos,
     1             rho0,vel0,gamma,xgeom,
     2             den,ener,pres,vel)

c       write(6,40) i,zpos,den,ener,pres,vel
       write(2,40) i,zpos,den,ener,pres,vel
 40    format(1x,i4,1p8e14.6)

      enddo

c..close up stop
      close(unit=2)
      end







      subroutine noh_1d(time,xpos,
     1                 rho1,u1,gamma,xgeom,
     2                 den,ener,pres,vel)
      implicit none
      save


c..solves the standard case, (as opposed to the singular or vacuum case), 
c..constant density (omega = 0) sedov problem in one-dimension.


c..input: 
c..time     = temporal point where solution is desired seconds
c..xpos     = spatial point where solution is desired cm


c..output:
c..den  = density  g/cm**3
c..ener = specific internal energy erg/g
c..pres = presssure erg/cm**3
c..vel  = velocity cm/sh


c..declare the pass
      double precision time,xpos,
     1                 rho1,u1,gamma,xgeom,
     3                 den,ener,pres,vel

c..local variables
      double precision gamm1,gamp1,gpogm,xgm1,us,r2,rhop,rho2,u2,e2,p2


c..some parameters
      gamm1 = gamma - 1.0d0
      gamp1 = gamma + 1.0d0
      gpogm = gamp1 / gamm1
      xgm1  = xgeom - 1.0d0


c..immediate post-chock values using strong shock relations
c..shock velocity, position, pre- and post-shock density,
c..flow velocity, internal energy, and pressure

      us   = 0.5d0 * gamm1 * abs(u1)                  
      r2   = us * time                                
      rhop = rho1 * (1.0d0 - (u1*time/r2))**xgm1
      rho2 = rho1 * gpogm**xgeom      
      u2   = 0.0d0
      e2   = 0.5d0 * u1**2
      p2   = gamm1 * rho2 * e2


c..if we are farther out than the shock front
      if (xpos .gt. r2) then
       den  = rho1 * (1.0d0 - (u1*time/xpos))**xgm1
       vel  = u1
       ener = 0.0d0
       pres = 0.0d0

c..if we are between the origin and the shock front
      else  
       den  = rho2
       vel  = u2
       ener = e2
       pres = p2
      end if

      return
      end





      double precision function value(string)
      implicit none
      save


c..this routine takes a character string and converts it to a real number. 
c..on error during the conversion, a fortran stop is issued

c..declare
      logical          pflag
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,se1,sd1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      double precision x,sign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   ,
     1                  blank = ' ' , se = 'e'    , sd = 'd'       ,
     2                  se1 = 'E'   , sd1 = 'D'   , rten =  10.0,
     3                  iten = 10                                   )

c..initialize
      x      =  0.0d0
      sign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)


c..remove any leading blanks and get the sign of the number
      do z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto  30
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         sign =  -1.0d0
        end if
        goto 10
       end if
      enddo


c..main number conversion loop
 10   continue
      do i = noblnk,long
       ipoint = i + 1


c..if a blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x     = x * sign
        value = x 
        return


c..if an exponent character, process the whole exponent, and return
       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or.
     1          string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0 .and. ipoint.eq.2)     x = 1.0d0
        if (sign .eq. -1.0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = sign * x * rten**(power*psign)
          value = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) goto 30
          power= (power * iten)  + j
         end if
        enddo


c..if an ascii number character, process ie
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) goto 30
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         goto 20
        end if

c..must be a decimal point if none of the above
c..check that there are not two decimal points!
       else
        if (pflag) goto 30
        pflag = .true.
       end if
 20   continue
      end do

c..if we got through the do loop ok, then we must be done
      x     = x * sign
      value = x 
      return
      

c..error processing the number
 30   write(6,40) long,string(1:long)
 40   format(' error converting the ',i4,' characters ',/,
     1       ' >',a,'< ',/,
     2       ' into a real number in function value')
      stop ' error in routine value'
      end



