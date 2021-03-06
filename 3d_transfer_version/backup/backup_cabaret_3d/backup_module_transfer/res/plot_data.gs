"set display color white"
"clear"
"cpt_gray"
"open tracer.ctl"
*"open tracer1.ctl"
"open upb.ctl"
"open vpb.ctl"
"open wpb.ctl"
"set yflip on"
tcf=1
tplot=20
lat=0
lon=-20
"set lat -90 90"
*"set lat "lat
"set lon 0 360"
"set lon "lon
"set z 3 78"
*"set lev 200000"
"set t 1"
"define u1=v.3"
"define u2=w.4"
*"d u"
"set t "tplot+1
"set gxout shaded"
"set grads off"
*"set clevs 0 0.1 0.2 0.3 0.4 0.5 1 2 3 4 5 10"
*"set ccols 0 36 34 32 30 28 26 24 22 20 19 18 17"
*"d maskout(u,-u)"
"d u/1e6"
"cbarn"
*"d ave(u,lat=-90,lat=90)"
"set gxout contour"
*"d u"
*"set clevs 0 1"
*"d u(t=1)/1e6"
*"set gxout stream"
"set gxout vector"
*"d hcurl(u.2(t=1),v.3(t=1))"
"d skip(u1,2,2);u2;mag(u1,u2)"
"set gxout contour"
"set cint 1"
"d mag(u1,u2)"
*"draw title T="tplot*tcf"h, lat="lat"`3."
"draw title T="tplot*tcf"h, lon="lon"`3."
*"printim adv_yz_lat"lat"_t"tplot".png"
