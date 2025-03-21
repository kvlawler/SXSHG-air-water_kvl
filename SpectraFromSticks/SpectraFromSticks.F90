Program SpectraFromSticks

!0) variables
integer :: ngrid,ierr,igrid,ipeak,which,BroadType
real*8 :: PeakE,f,pi,fwhm,sigma,hwhm,a,c,Elo,Ehi,dE,thisE,E2,AbsUnit,cfa,maxI
real*8, allocatable :: spectra(:)
integer, parameter :: stdin = 5
real*8, parameter :: hc = 1239.8419300924
logical :: found,tonm,normalize
character(len=80) :: tmp,line
character*10 :: Etag,Efwhm,label
NAMELIST /runtype/ which
NAMELIST /params/ BroadType,Elo,Ehi,dE,fwhm,fwhm_hi,tonm,normalize

! *   which = spectrum type = 0 Electronic, 1 IR, 2 Raman        *  
! *   BroadType = 1 Gaussian, 2 Lorentzian, Default off          * 
! *     + 10 converts into experimental units L mol^-1 cm^-1     *
! *   fwhm = Full-width at half max                              *
! *     >0 is a static broadening where FWHM=read value          *
! *     <0 is a dynamic broadening where FWHM=read value*PeakE   *

!1) initialize
!1.1)  read the spectrum type and initialize defaults 
! AbsUnit is for a normalized Gaussian in L mol^-1 cm^-2
!   which converts oscillator strength to absorptivity
!   and it comes from N_a e^2 / (2 m_e c^2 epsilon_0) *ln 10)
! If electronic spectra convert AbsUnit to L mol^-1 cm^-1 eV
!   & 1 eV/8065.5 cm^-1
AbsUnit = 2.1751279267d8
pi = 4d0 * atan(1d0) !4 times pi/4 to machine precision
which = 0 !the run default is electronic it is electronic
read(stdin,nml=runtype,iostat=ierr)
if(ierr.ne.0) then
  write(*,*) "runtype namelist not found, and this is needed for operation"
  write(*,runtype) ! dump the namelist for debugging purposes
  stop
endif
select case(which)
case(0)
  fwhm = 1d-1
  dE = 1d-2
  Elo = 0d0
  Ehi = 1d3
  Etag = "eV"
  label = "ELECTRONIC" 
  AbsUnit = AbsUnit / 8065.54429d0  
case(1) 
  fwhm = 5d1
  Elo = 4d2
  Ehi = 4d3
  dE = 1d1
  Etag = "cm-1"
  label = "IR"
case(2) 
  fwhm = 5d1
  Elo = 4d2
  Ehi = 4d3
  dE = 1d0
  Etag = "cm-1"
  label = "RAMAN"
case default
  write(*,*) "Unknown spectrum type selected"
  stop
end select
tonm = .false.
Efwhm = Etag
normalize = .false.

!1.2) read the working parameters and error check
read(stdin,nml=params,iostat=ierr)
if(ierr.ne.0) then
  write(*,*) "params namelist not found, and this is needed for operation"
  write(*,params) ! dump the namelist for debugging purposes
  stop
endif
if(tonm) then
  if(which.eq.0) then
    Etag = 'nm'
  else
    tonm = .false.
  endif
endif

!1.3) Error check grid size
if(Elo.ge.Ehi) then
  write(*,*) "Energy range input backwards, aborting"
  stop
endif
ngrid=ceiling((Ehi-Elo)/dE) + 1
if(ngrid.le.1) then
  write(*,*) "Spectra will have less than 2 points, please fix inputs"
  stop
endif
!write(*,*) "Ehi is ",Ehi
!write(*,*) "Elo is ",Elo
!write(*,*) "dE  is ",dE
!write(*,*) "this gives a grid size: ",ngrid
allocate(spectra(ngrid),stat=ierr)
if(ierr.ne.0) then
  write(*,*) "Error allocating spectra array"
  stop
endif
spectra(:)=0d0

!1.4) Error check FWHM details ***OLD WAY OF INCREASING FWHM as a function of E
!if(fwhm.le.0d0) then
!  write(*,*) "ERROR: fwhm <= 0, aborting" 
!  stop
!elseif(fwhm.gt.Ehi) then
!  write(*,*) "ERROR: fwhm > Ehi, aborting" 
!  stop
!endif
!dFW=0d0
!if(fwhm_hi.gt.fwhm) then
!  dFW = (fwhm_hi-fwhm) / (ngrid-1)
!  write(*,*) "The broadening will expand linearly from",fwhm," to ",fwhm_hi
!endif

!2) read the peaks and build them into the spectra array
!2.1) find the PEAKS section in stdin
found=.false.
rewind(stdin)
do while(.not.found)
  read(stdin,'(A)',END=10) tmp
  line=trim(adjustl(tmp))
  if(line(:5).eq.'PEAKS') found=.true.
enddo
10 continue
if(.not.found) then
  write(*,*) "Did not find the header PEAKS that indicates where to start reading comuted peaks, aborting"
  stop
endif 

!2.2) Read peaks and build them into spectra
select case(BroadType)
!2.2.1) Gaussian: f/sigma/sqrt(2pi) * exp(-0.5*(E-E0)^2/sigma^2)
  case(1)
    if(fwhm.lt.0) then
      write(*,1010) "Gaussian",abs(fwhm)*1d2
    else
      write(*,1000) "Gaussian",fwhm,trim(adjustl(Efwhm))
    endif
    write(*,1020) trim(adjustl(Etag))
    ipeak=0
    do
      ipeak=ipeak+1
      read(stdin,*,END=20) PeakE,f
      !write(*,*) "peak: ",ipeak,", Energy=",PeakE," eV f=",f
      if(PeakE.lt.Elo) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is lower than the energy window"
        stop
      elseif(PeakE.gt.Ehi) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is greater than the energy window"
        stop
      endif
      if(fwhm.lt.0) then
        sigma = (abs(fwhm)*PeakE) / (2d0 * sqrt(2d0*log(2d0)))
      else 
        sigma = fwhm / (2d0 * sqrt(2d0*log(2d0))) 
      endif
      a = 1d0 / sigma / sqrt(2d0*pi)
      c = -5d-1 / sigma / sigma
      do igrid = 1,ngrid
        thisE = Elo + ((igrid-1) * dE)
        E2 = c * ((thisE - PeakE)**2) 
        spectra(igrid) = spectra(igrid) + (f * a * exp(E2))
      enddo
    enddo
    20 continue

!2.2.2) Lorentzian: f/pi/hwhm * [ hwhm^2 / ((E-E0)^2 + hwhm^2)]
!                 = f/pi*hwhm / ((E-E0)^2 + hwhm^2)
  case(2)
    if(fwhm.lt.0) then
      write(*,1010) "Lorentzian",abs(fwhm)*1d2
    else
      write(*,1000) "Lorentzian",fwhm,trim(adjustl(Efwhm))
    endif
    write(*,1020) trim(adjustl(Etag))
    do
      ipeak=ipeak+1
      read(stdin,*,END=30) PeakE,f
      !write(*,*) "peak: ",ipeak,", Energy=",PeakE," eV f=",f
      if(PeakE.lt.Elo) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is lower than the energy window"
        stop
      elseif(PeakE.gt.Ehi) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is greater than the energy window"
        stop
      endif
      if(fwhm.lt.0) then
        hwhm = (abs(fwhm)*PeakE) * 5d-1 
      else
        hwhm = fwhm * 5d-1 
      endif
      a = hwhm / pi
      c = hwhm * hwhm
      do igrid = 1,ngrid
        thisE = Elo + ((igrid-1) * dE)
        E2 = c + ((thisE - PeakE)**2) 
        spectra(igrid) = spectra(igrid) + (f * a / E2)
      enddo
    enddo
    30 continue

!2.2.11) Gaussian in exptl units : AbsUnit*f/FWHM * exp(-0.5*(E-E0)^2/sigma^2)
  case(11)
    if(fwhm.lt.0) then
      write(*,1010) "Gaussian",abs(fwhm)*1d2
    else
      write(*,1000) "Gaussian",fwhm,trim(adjustl(Efwhm))
    endif
    write(*,1030) trim(adjustl(Etag))
    do
      ipeak=ipeak+1
      read(stdin,*,END=40) PeakE,f
      !write(*,*) "peak: ",ipeak,", Energy=",PeakE," eV f=",f
      if(PeakE.lt.Elo) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is lower than the energy window"
        stop
      elseif(PeakE.gt.Ehi) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is greater than the energy window"
        stop
      endif
      if(fwhm.lt.0) then
        a = abs(fwhm) * PeakE
      else 
        a = fwhm
      endif
      cfa = AbsUnit * (f / a)
      c = -4d0 * log(2d0) / a**2
      do igrid = 1,ngrid
        thisE = Elo + ((igrid-1) * dE)
        E2 = c * ((thisE - PeakE)**2) 
        spectra(igrid) = spectra(igrid) + (cfa * exp (E2)) 
      enddo
    enddo
    40 continue

!2.2.12) Lorentzian in exptl units: AbsUnit*f/fwhm * [ hwhm^2 / ((E-E0)^2 + hwhm^2)]
!          AbsUnit for Lorentzian = SQRT(pi*ln(2)) * AbsUnit for Gaussian
  case(12)
    if(fwhm.lt.0) then
      write(*,1010) "Lorentzian",abs(fwhm)*1d2
    else
      write(*,1000) "Lorentzian",fwhm,trim(adjustl(Efwhm))
    endif
    write(*,1030) trim(adjustl(Etag))
    do
      ipeak=ipeak+1
      read(stdin,*,END=50) PeakE,f
      !write(*,*) "peak: ",ipeak,", Energy=",PeakE," eV f=",f
      if(PeakE.lt.Elo) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is lower than the energy window"
        stop
      elseif(PeakE.gt.Ehi) then
        write(*,*) "ERROR: read excitation energy ",ipeak," is greater than the energy window"
        stop
      endif
      if(fwhm.lt.0) then
        a = abs(fwhm) * PeakE
      else 
        a = fwhm
      endif
      cfa = AbsUnit * SQRT(pi * log(2d0)) * (f / a)
      hwhm = a * 5d-1 
      do igrid = 1,ngrid
        thisE = Elo + ((igrid-1) * dE)
        E2 = 1d0 + (((thisE - PeakE) / hwhm)**2) 
        spectra(igrid) = spectra(igrid) + (cfa / E2)
      enddo
    enddo
    50 continue

  case default
    write(*,*) "Unknown functional form for spectral broadening"
    stop
end select

!3) Dump the computed spectra 
!3.1) find max value if needed
maxI=1d0
if(normalize) then
  do igrid=1,ngrid
    if(spectra(igrid).gt.maxI) then
      maxI=spectra(igrid)
    endif
  enddo
endif
!3.2) actually write the spectra
if(tonm) then
  do igrid=ngrid,1,-1
    Egrid = Elo + ((igrid-1)*dE)
    write(*,'(2F20.10)') hc/Egrid,spectra(igrid)/maxI
  enddo
else
  do igrid=1,ngrid
    thisE = Elo + ((igrid-1)*dE)
    write(*,*) thisE,spectra(igrid)/maxI
  enddo
endif

!4) clean up and leave
deallocate(spectra,stat=ierr)
if(ierr.ne.0) then
  write(*,*) "Error deallocating spectra array"
  stop
endif
1000  FORMAT(A,' broadening with full width at half max of ',F7.3,x,A4)
1010  FORMAT('Dynamic ',A,' broadening with full width at half max of ',F7.3,'% peak Energy')
1020  FORMAT(8x,'Energy(',A,')',14x,'Intensity (arb. units)')
1030  FORMAT(8x,'Energy(',A,')',14x,'Intensity (L mol^-1 cm^-1)')
End Program SpectraFromSticks
