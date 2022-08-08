        subroutine azdist(stalat, stalon, evtlat, evtlon, delta, az, baz)
!
! Subroutine to calculate the Great Circle Arc distance
!    between two sets of geographic coordinates
!
! Given:  stalat => Latitude of first point (+N, -S) in degrees
!         stalon => Longitude of first point (+E, -W) in degrees
!         evtlat => Latitude of second point
!         evtlon => Longitude of second point
!
! Returns:  delta => Great Circle Arc distance in degrees
!           az    => Azimuth from pt. 1 to pt. 2 in degrees
!           baz   => Back Azimuth from pt. 2 to pt. 1 in degrees
!
! If you are calculating station-epicenter pairs, pt. 1 is the station
!
! Equations take from Bullen, pages 154, 155
!
! T. Owens, September 19, 1991
!           Sept. 25 -- fixed az and baz calculations
!           Dec. 2006, changed for fortran95 
!           May, 2007 -- added predel to get around OSX acos round-off NaN issue
!
      !double precision scolat, slon, ecolat, elon
      !double precision a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk
      !double precision rhs1,rhs2,sph,rad,del,daz,dbaz,pi
      
      real*8 :: scolat, slon, ecolat, elon
      real*8 :: a,b,c,d,e,aa,bb,cc,dd,ee,g,gg,h,hh,k,kk
      real*8 :: rhs1,rhs2,sph,rad,del,daz,dbaz,pi
!
      pi=4.0*atan(1.0)
      piby2=pi/2.
      rad=2.*pi/360.
!
! scolat and ecolat are the geocentric colatitudes
! as defined by Richter (pg. 318)
!
! Earth Flattening of 1/298.257 take from Bott (pg. 3)
!
      sph=1.0/298.257
!
      scolat=piby2 - atan((1.-sph)*(1.-sph)*tan(dble(stalat)*rad))
      ecolat=piby2 - atan((1.-sph)*(1.-sph)*tan(dble(evtlat)*rad))
      slon=dble(stalon)*rad
      elon=dble(evtlon)*rad
!
!  a - e are as defined by Bullen (pg. 154, Sec 10.2)
!     These are defined for the pt. 1
!
      a=sin(scolat)*cos(slon)
      b=sin(scolat)*sin(slon)
      c=cos(scolat)
      d=sin(slon)
      e=-cos(slon)
      g=-c*e
      h=c*d
      k=-sin(scolat)
!
!  aa - ee are the same as a - e, except for pt. 2
!
      aa=sin(ecolat)*cos(elon)
      bb=sin(ecolat)*sin(elon)
      cc=cos(ecolat)
      dd=sin(elon)
      ee=-cos(elon)
      gg=-cc*ee
      hh=cc*dd
      kk=-sin(ecolat)
!
!  Bullen, Sec 10.2, eqn. 4
!
      predel=a*aa + b*bb + c*cc
      if(abs(predel+1.).lt..000001) then
        predel=-1.
      endif
      if(abs(predel-1.).lt..000001) then
        predel=1.
      endif
      del=acos(predel)
      delta=del/rad
!
!  Bullen, Sec 10.2, eqn 7 / eqn 8
!
!    pt. 1 is unprimed, so this is technically the baz
!
!  Calculate baz this way to avoid quadrant problems
!
      rhs1=(aa-d)*(aa-d)+(bb-e)*(bb-e)+cc*cc - 2.
      rhs2=(aa-g)*(aa-g)+(bb-h)*(bb-h)+(cc-k)*(cc-k) - 2.
      dbaz=atan2(rhs1,rhs2)
      if(dbaz.lt.0.0d0) dbaz=dbaz+2*pi
      baz=dbaz/rad
!
!  Bullen, Sec 10.2, eqn 7 / eqn 8
!
!    pt. 2 is unprimed, so this is technically the az
!
      rhs1=(a-dd)*(a-dd)+(b-ee)*(b-ee)+c*c - 2.
      rhs2=(a-gg)*(a-gg)+(b-hh)*(b-hh)+(c-kk)*(c-kk) - 2.
      daz=atan2(rhs1,rhs2)
      if(daz.lt.0.0d0) daz=daz+2*pi
      az=daz/rad
!
!   Make sure 0.0 is always 0.0, not 360.
!
      if(abs(baz-360.).lt..00001) baz=0.0
      if(abs(az-360.).lt..00001) az=0.0
      return
      end
