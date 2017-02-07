PRO overlay
fitsname='/nfshome/lqian/clumps/t13_new.fits';default file name
fitsname2='/nfshome/lqian/clumps/t12_new.fits'
fitsname3='/nfshome/lqian/clumps/HI.fits'



G = 6.67e-8; gravitational constant
c=3.0d10; speed of light, units in cm/s
ckm=3.0d5; speed of light, units in km/s
mu=0.122d0; dipole of CO molecular, unit in Debye, 1 Debye = 10^-18 e.s.u.
;mu=mu*1e-18; dipole of CO molecular, e.s.u. units
pi=3.1415926d0
h=6.63d-27; Planck constant in cgs units
nu12=115.271203d0; Frequency of 12CO J 1->0 transition, units in GHz
nu13=110.201370d0; Frequency of 13CO J 1->0 transition, units in GHz
A12=64.0*pi^4/3.0/h*(nu12/c)^3*mu^2*(1.0/3.0); Transition probability rate of 12CO J 1->0
A13=64.0*pi^4/3.0/h*(nu13/c)^3*mu^2*(1.0/3.0); Transition probability rate of 13CO J 1->0
A12=A12*(1.0d9)^3*(1.0d-18)^2; Frequency in GHz, dipole in Debye
A13=A13*(1.0d9)^3*(1.0d-18)^2; Frequency in GHz, dipole in Debye

;print,nu13^2/A13

;deal,fitsname,fitsname2,nu13,A13,28
deal,fitsname,nu13,A13,29,fitsname2,nu12,A12,28,fitsname3

END

PRO deal,fitsname,nu13,Aul13,weight13,fitsname2,nu12,Aul12,weight12,fitsname3
print,Aul13

;PRO deal,fitsname,fitsname2,nu,Aul13,weight13
G = 6.67e-8; gravitational constant
h=6.63d-27; Planck constant in cgs units
c=3.0d10; speed of light, units in cm/s
ck=1.38e-16; Boltzmann constant
ckm=3.0d5; speed of light, units in km/s
mu=0.122d0; dipole of CO molecular, unit in Debye, 1 Debye = 10^-18 e.s.u.
;mu=mu*1e-18; dipole of CO molecular, e.s.u. units
pi=3.1415926d0

head=headfits(fitsname);read the header of the fits file to a vector

bw      =   fxpar(head,'BW'); band width (MHz)
freq    =   fxpar(head,'LINEFREQ'); central frequency (GHz)

nx      =   fxpar(head,'NAXIS1'); number of elements in the first dimension
ny      =   fxpar(head,'NAXIS2');
nz      =   fxpar(head,'NAXIS3');
crvalx  =   fxpar(head,'CRVAL1'); reference value of the first dimension
cdeltax =   fxpar(head,'CDELT1'); increasement of the first dimension
                                ; in units of degree, when calculate physical
				; scale, must changed to arcdegree
crpixx  =   fxpar(head,'CRPIX1'); reference position of the first dimension
crvaly  =   fxpar(head,'CRVAL2');
cdeltay =   fxpar(head,'CDELT2');
crpixy  =   fxpar(head,'CRPIX2');
crvalz  =   fxpar(head,'CRVAL3');
cdeltaz =   fxpar(head,'CDELT3');
crpixz  =   fxpar(head,'CRPIX3');
ctype1  =   fxpar(head,'CTYPE1');
ctype2  =   fxpar(head,'CTYPE2');
print,ctype1,ctype2

x=(dindgen(nx)+0.5-crpixx)*cdeltax+crvalx;
y=(dindgen(ny)+0.5-crpixy)*cdeltay+crvaly;
z=(dindgen(nz)+0.5-crpixz)*cdeltaz+crvalz;


Distance = 140.0; distance in units of parcec
Distance = Distance*3.086d18; distance in units of cm
vchannel=abs(cdeltaz)/1.0e3; km/s?
area = Distance^2*(abs(cdeltax)/57.3)*(abs(cdeltay)/57.3)
;13CO

a=mrdfits(fitsname)

head2=headfits(fitsname2);read the header of the fits file to a vector
bw2      =   fxpar(head2,'BW'); band width (MHz)
freq2    =   fxpar(head2,'LINEFREQ'); central frequency (GHz)


nx2      =   fxpar(head2,'NAXIS1'); number of elements in the first dimension
ny2      =   fxpar(head2,'NAXIS2');
nz2      =   fxpar(head2,'NAXIS3');
crvalx2  =   fxpar(head2,'CRVAL1'); reference value of the first dimension
cdeltax2 =   fxpar(head2,'CDELT1'); increasement of the first dimension
                                ; in units of degree, when calculate physical
				; scale, must changed to arcdegree
crpixx2  =   fxpar(head2,'CRPIX1'); reference position of the first dimension
crvaly2  =   fxpar(head2,'CRVAL2');
cdeltay2 =   fxpar(head2,'CDELT2');
crpixy2  =   fxpar(head2,'CRPIX2');
crvalz2  =   fxpar(head2,'CRVAL3');
cdeltaz2 =   fxpar(head2,'CDELT3');
crpixz2  =   fxpar(head2,'CRPIX3');
x2=(dindgen(nx2)+0.5-crpixx2)*cdeltax2+crvalx2;
y2=(dindgen(ny2)+0.5-crpixy2)*cdeltay2+crvaly2;
z2=(dindgen(nz2)+0.5-crpixz2)*cdeltaz2+crvalz2;
a2=mrdfits(fitsname2)
vchannel2=abs(cdeltaz2)/1.0e3
;12CO


head3=headfits(fitsname3);read the header of the fits file to a vector

bw3      =   fxpar(head3,'BW'); band width
freq3    =   fxpar(head3,'LINEFREQ'); central frequency 

nx3      =   fxpar(head3,'NAXIS1'); number of elements in the first dimension
ny3      =   fxpar(head3,'NAXIS2');
nz3      =   fxpar(head3,'NAXIS3');
crvalx3  =   fxpar(head3,'CRVAL1'); reference value of the first dimension
cdeltax3 =   fxpar(head3,'CDELT1'); increasement of the first dimension
                                ; in units of degree, when calculate physical
				; scale, must changed to arcdegree
crpixx3  =   fxpar(head3,'CRPIX1'); reference position of the first dimension
crvaly3  =   fxpar(head3,'CRVAL2');
cdeltay3 =   fxpar(head3,'CDELT2');
crpixy3  =   fxpar(head3,'CRPIX2');
crvalz3  =   fxpar(head3,'CRVAL3');
cdeltaz3 =   fxpar(head3,'CDELT3');
crpixz3  =   fxpar(head3,'CRPIX3');

x3=(dindgen(nx3)+0.5-crpixx3)*cdeltax3+crvalx3;
y3=(dindgen(ny3)+0.5-crpixy3)*cdeltay3+crvaly3;
z3=(dindgen(nz3)+0.5-crpixz3)*cdeltaz3+crvalz3;


b=mrdfits(fitsname3)

t12=total(a2,3)
bad=where(finite(t12) eq 0, count)
if(count gt 0) then t12(bad)=0.0 ; here just define t12 as a 2D array
                                 ; value not important
t13=total(a,3)

for i=0,nx-1 do begin
  for j=0,ny-1 do begin
      t12[i,j]=max(a2(i,j,*))/0.47;correct for the beam efficiency
      t13[i,j]=max(a(i,j,*))/0.47
  end
end

;tk=t12+2.7
tbg=2.7
tk=1.0/Alog(1.0/(ck*t12/h/(nu12*1.0e9)+1.0/(exp(h*(nu12*1.0e9)/ck/tbg)-1.0)) $
   +1.0)*h*nu12*1.0e9/ck
;print,max(tk)
f1=t13/t12
tau13=-Alog(1.0-f1)
ftau=tau13/f1
bad=where(finite(ftau) eq 0, count)
if(count gt 0) then ftau(bad)=1.0
;print,ftau
fb=1.0/(1.0-(exp(h*nu13/ck/tk)-1.0)/(exp(h*nu13/ck/2.7)-1.0))
;print,fb*ftau

QT=tk/(h*nu13*1.0e9/ck/2.0)
fu=QT/(3.0*exp(-h*nu13*1.0e9/ck/tk))

nt=total(a,3)*vchannel*1.93e3*nu13^2/Aul13*(1+1.0/3.0*exp(h*nu13/ck/tk))*fb*ftau
;*****Note: correct is an array!!!************** 
;correct=fb*ftau*(1+1.0/3.0*exp(h*nu13/ck/tk))
correct=fb*ftau*fu
;***********************************************
image=total(a,3)
;image=total(a[*,*,37:38],3)
image=(abs(image)+image)/2.0
;image1=image*vchannel*1.93e3*nu13^2/Aul13*correct/1.7e-6
;image1=image1/1.0d19
image1=image*0.266
image2=image*0.266
print,max(image1)
plotscale=3.0
bad=where(image2 gt max(image2)/plotscale, count)
if(count gt 0) then image2(bad)=max(image2)/plotscale


ncolor=250
LoadCT, 0, NColors=ncolor, Bottom=1, /Silent
tvlct,rn,gn,bn,/get 
rr=reverse(rn) 
gg=reverse(gn) 
bb=reverse(bn) 
tvlct,rr,gg,bb 

device,decompose=0

xrange=[x[0],x[nx-1]]
yrange=[y[0],y[ny-1]]
position = [0.15, 0.0, 0.85, 1.0]
xrange2=[x2[0],x2[nx2-1]]
yrange2=[y2[0],y2[ny2-1]]

margin=0.12
wall=0.03
xsize=8.8
aa=xsize/8.8-(margin+wall)
bb=aa*2d/(1+sqrt(5))
ysize=(margin+bb+wall+bb+wall)*8.8

xsize=10.0
ysize=8.8

openr,lun,'txt/catalog.txt',/get_lun

;================================================
set_plot,'PS'
filename='overlay.eps' ; set the file name of the output ps file
device,file=filename,/ENCAPSULATED,/COLOR, BITS=8,xsize=xsize,ysize=ysize,set_font='Times-Roman'

!p.font=-1
!x.style=1
!y.style=1
!p.charsize=1.0
!p.thick=1.8
!p.charthick=1.8
!x.thick=1.8
!y.thick=1.8
!z.thick=1.8

TVimage,BytScl(image2, Top=ncolor),x,y,Position=position, $
/Keep_Aspect,/Erase, /NoInterpolate
rastr=TeXtoIDL('RA(J2000)')
decstr=TeXtoIDL('Dec(J2000)')
statstr=TeXtoIDL('!1704^h48^m')
ystr1=TeXtoIDL('!1722^\circ')
ystr2=TeXtoIDL('!1724^\circ')
ystr3=TeXtoIDL('!1726^\circ')
ystr4=TeXtoIDL('!1728^\circ')
ystr5=TeXtoIDL('!1730^\circ')
!X.TICKNAME=[statstr,'!17','!1732','!1724','!1716']
!Y.TICKNAME=[ystr1,ystr2,ystr3,ystr4,ystr5]
unitstr=TeXtoIDL('!17K\cdot km/s')
Plot, xrange, yrange, XRANGE=xrange, YRANGE=yrange, $
Position=position, XStyle=1, YStyle=1, /NoData, /NoErase, $
xtitle='!17RA', $
ytitle='!17Dec',color=ncolor,ticklen=0.0
Colorbar, Range=[Min(image1), Max(image1)],Charsize=1,charthick=2, $
Divisions=3, Minor=5, NColors=ncolor, Bottom=1,color=ncolor, $
Position=[0.9, 0.20, 0.92, 0.80],/VERTICAL
;!p.charsize=0.1
xyouts,62.4,30.4,unitstr,Charsize=0.7,charthick=2,color=ncolor
;!p.charsize=1.0

ncolor=250
LoadCT, 39, NColors=ncolor, Bottom=1, /Silent

liney=y
linex=(y*0.0+64.0-crvalx)*cos(y*!pi/180.0)+crvalx
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=y
linex=(y*0.0+66.0-crvalx)*cos(y*!pi/180.0)+crvalx
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=y
linex=(y*0.0+68.0-crvalx)*cos(y*!pi/180.0)+crvalx
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=y
linex=(y*0.0+70.0-crvalx)*cos(y*!pi/180.0)+crvalx
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=y
linex=(y*0.0+72.0-crvalx)*cos(y*!pi/180.0)+crvalx
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=x*0.0+24.0
linex=x
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=x*0.0+26.0
linex=x
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors
liney=x*0.0+28.0
linex=x
;oplot,linex,liney,linestyle=1,thick=3,color=0.90*!d.n_colors

liney=(dindgen(1529-400)+400+0.5-crpixy)*cdeltay+crvaly
linex=liney*0.0+(300+0.5-crpixx)*cdeltax+crvalx
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
liney=(dindgen(1150)+0.5-crpixy)*cdeltay+crvaly
linex=liney*0.0+(600+0.5-crpixx)*cdeltax+crvalx
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
liney=(dindgen(1150-400)+400+0.5-crpixy)*cdeltay+crvaly
linex=liney*0.0+(950+0.5-crpixx)*cdeltax+crvalx
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
liney=(dindgen(1529-400)+400+0.5-crpixy)*cdeltay+crvaly
linex=liney*0.0+(1250+0.5-crpixx)*cdeltax+crvalx
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
;liney=(dindgen(1529)+0.5-crpixy)*cdeltay+crvaly
;linex=liney*0.0+(1500+0.5-crpixx)*cdeltax+crvalx
;oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
liney=(dindgen(1529)+0.5-crpixy)*cdeltay+crvaly
linex=liney*0.0+(1750+0.5-crpixx)*cdeltax+crvalx
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
linex=(dindgen(1750)+0.5-crpixx)*cdeltax+crvalx
liney=linex*0.0+(400+0.5-crpixy)*cdeltay+crvaly ;;note
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
linex=(dindgen(1250)+0.5-crpixx)*cdeltax+crvalx
liney=linex*0.0+(1150+0.5-crpixy)*cdeltay+crvaly ;;note
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
linex=(dindgen(2069-1750)+1750+0.5-crpixx)*cdeltax+crvalx
liney=linex*0.0+(500+0.5-crpixy)*cdeltay+crvaly
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
linex=(dindgen(2069-1750)+1750+0.5-crpixx)*cdeltax+crvalx
liney=linex*0.0+(1000+0.5-crpixy)*cdeltay+crvaly
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
;linex=(dindgen(1750-1250)+1250+0.5-crpixx)*cdeltax+crvalx
;liney=linex*0.0+(400+0.5-crpixy)*cdeltay+crvaly;;note
;oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
linex=(dindgen(1750-1250)+1250+0.5-crpixx)*cdeltax+crvalx
liney=linex*0.0+(750+0.5-crpixy)*cdeltay+crvaly
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
linex=(dindgen(1750-1250)+1250+0.5-crpixx)*cdeltax+crvalx
liney=linex*0.0+(1100+0.5-crpixy)*cdeltay+crvaly
oplot,linex,liney,linestyle=0,thick=3,color=0.30*!d.n_colors
xyouts, 73.2, 23, '!17 1',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 73.2, 26, '!17 2',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 73.2, 29, '!17 3',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
;xyouts, 71.5, 23, '!17 4',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 71.5, 26, '!17 4',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
;xyouts, 71.5, 29, '!17 6',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
;xyouts, 70, 23, '!17 7',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 70, 26, '!17 6',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 70, 29, '!17 5',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
;xyouts, 68.2, 22.5, '!17 9',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 68.2, 25.5, '!17 7',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
;xyouts, 68.2, 29, '!17 12',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 66.7, 23, '!17 8',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 66.7, 25, '!17 9',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 66.7, 26, '!17 10',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 66.7, 29, '!17 11',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
;xyouts, 65.5, 24, 'XVI',charsize=1.0,charthick=1.5,color=0.30*!d.n_colors
;xyouts, 65.5, 27, 'XVII',charsize=1.0,charthick=1.5,color=0.30*!d.n_colors
;xyouts, 65.5, 29, 'XVIII',charsize=1.0,charthick=1.5,color=0.30*!d.n_colors
xyouts, 64.0, 23, '!17 12',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 64.0, 26, '!17 13',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors
xyouts, 64.0, 29, '!17 14',charsize=1.0,charthick=1.5,color=0.90*!d.n_colors

major=1.0d
minor=1.0d
angle=300.0d
ra=65.0d
dec=27.0d

convert=3.1415926/180.0
ra0=crvalx
dec0=crvaly

totmass=0.0d
coremass=0.0d
countlong=0
while not eof(lun)do begin
    readf,lun,ID,peak1,peak2,peak3,cen1,cen2,cen3,size1,size2,size3,sum,peak,volume,fwhm1,fwhm2,fwhm3,ra,dec,major,minor,angle
yd=dec
xr=(ra-ra0)*cos(dec*convert)+ra0
nxr=floor((xr-crvalx)/cdeltax+crpixx)
nyd=floor((yd-crvaly)/cdeltay+crpixy)
;print,xr,yd,nxr,nyd
;sum=sum/0.47;correct for the beam efficiency
;peak=peak/0.47;correct for the beam efficiency
;cmass = sum*area*vchannel*1.93e3*nu13^2/Aul13*correct(nxr,nyd)*1.67d-24*2.0/2.0d33/1.7e-6
;virmass = (3*(size3*1.0e2)^2)*(Distance*(minor+major)/2.0*convert)/G/2.0d33
;temperature = tk[nxr,nyd]

;printf,lun2,ra,' & ',dec,' & ',major*60.0,' & ',minor*60.0,' & ',angle,' & ',temperature,' & ',cmass,' & ',virmass,' & ',size3/1000.0*2*sqrt(2.0*alog(2.0)),' \\',format='(2(f8.2,A3),2(f8.3,A3),5(f8.2,A3))'
;printf,lun3,ra,dec,major*60.0,minor*60.0,angle,temperature,cmass,size3/1000.0*2*sqrt(2.0*alog(2.0)),format='(2(f8.2,A3),2(f8.3,A3),4(f8.2,A3))'

;tangra=-sin(dec*convert)*sin((dec-dec0)*convert)
;tangdec=cos(dec*convert)*cos((dec-dec0)*convert)
angle1=(90-angle)

if (major/minor lt 10.0 and peak gt 0.7) then begin
;if (major/minor lt 10.0 and peak gt 2.0) then begin
;if (peak gt 0.5) then begin
;ellipse,major,minor,angle1,0.0,360.0,xr,yd,color=ncolor*4/5
;if (major/minor gt 10.0) then ellipse,major,minor,angle1,0.0,360.0,xr,yd,color=ncolor*4/5
if (major/minor lt 10.0) then ellipse,major,minor,angle1,0.0,360.0,xr,yd,color=ncolor*3/5
coremass=coremass+sum
endif else begin
countlong=countlong+1
endelse
totmass=totmass+sum

endwhile
free_lun,lun

;openr,lun,'txt/catalog_sup.txt',/get_lun
openr,lun,'txt/catalog.txt',/get_lun
while not eof(lun)do begin
    readf,lun,ID,peak1,peak2,peak3,cen1,cen2,cen3,size1,size2,size3,sum,peak,volume,fwhm1,fwhm2,fwhm3,ra,dec,major,minor,angle
yd=dec
xr=(ra-ra0)*cos(dec*convert)+ra0
nxr=floor((xr-crvalx)/cdeltax+crpixx)
nyd=floor((yd-crvaly)/cdeltay+crpixy)
;print,xr,yd,nxr,nyd
;sum=sum/0.47;correct for the beam efficiency
;peak=peak/0.47;correct for the beam efficiency
;cmass = sum*area*vchannel*1.93e3*nu13^2/Aul13*correct(nxr,nyd)*1.67d-24*2.0/2.0d33/1.7e-6
;virmass = (3*(size3*1.0e2)^2)*(Distance*(minor+major)/2.0*convert)/G/2.0d33
;temperature = tk[nxr,nyd]

;printf,lun2,ra,' & ',dec,' & ',major*60.0,' & ',minor*60.0,' & ',angle,' & ',temperature,' & ',cmass,' & ',virmass,' & ',size3/1000.0*2*sqrt(2.0*alog(2.0)),' \\',format='(2(f8.2,A3),2(f8.3,A3),5(f8.2,A3))'
;printf,lun3,ra,dec,major*60.0,minor*60.0,angle,temperature,cmass,size3/1000.0*2*sqrt(2.0*alog(2.0)),format='(2(f8.2,A3),2(f8.3,A3),4(f8.2,A3))'

;tangra=-sin(dec*convert)*sin((dec-dec0)*convert)
;tangdec=cos(dec*convert)*cos((dec-dec0)*convert)
angle1=(90-angle)

if (major/minor lt 10.0 and peak gt 0.5) then begin
;ellipse,major,minor,angle1,0.0,360.0,xr,yd,color=ncolor*4/5
;ellipse,major,minor,angle1,0.0,360.0,xr,yd,color=ncolor*4/5
coremass=coremass+sum
endif else begin
countlong=countlong+1
endelse
totmass=totmass+sum

endwhile

free_lun,lun

device,/CLOSE

;================================================
print,coremass/totmass
print,countlong



END
