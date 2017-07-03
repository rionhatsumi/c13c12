# plot penetration factors for carbon reactions

set format y "%.2e"

hbarc=197.32

# Coulomb barrier calculation based on 
# Woods-Saxon Equivalent to a Double Folding Potential
# by Braz J Phys (2016) 46:120–128

ffxx(z1,a1,z2,a2)=27.1*(a1**0.333+a2**0.3333)**2/z1/z2
frb(z1,a1,z2,a2)=1.3*(a1**0.3333+a2**0.3333)+0.65*log(ffxx(z1,a1,z2,a2))
fvb(z1,a1,z2,a2)=z1*z2*1.44/frb(z1,a1,z2,a2)-15/(ffxx(z1,a1,z2,a2)+1)
print "c12c12 Coulomb barrier height : ", fvb(6,12,6,12)

mu_c12c12=6.0
mu_c12c13=12*13/25.0
mu_c13c13=6.5
mu_c12c14=12*14/26.0
mu_c12o16=12*16/28.
mu_c13o16=13*16/29.0
mu_o16o16=8.0
mu_b11c12=11*12/23.
mu_b10c12=10*12/22.

rbc12c12=8.04
rbc12c13=8.15
rbc13c13=8.22
rbc13o16=8.27271

vbb10c12=4.98
vbb11c12=4.91
vbc12c12=5.93
vbc12c13=5.85
vbc13c13=5.79
vbc12c14=5.79
vbc12o16=7.76
vbc13o16=7.65649
vbo16o16=10.17

#########################################
# Coulomb barrier based on formula
# ######################################
vbb10c12=fvb(5,10,6,12)
vbb11c12=fvb(5,11,6,12)

vbc12c12=fvb(6,12,6,12)
vbc12c13=fvb(6,12,6,13)
vbc13c13=fvb(6,13,6,13)
vbc12c14=fvb(6,12,6,14)

vbc12o16=fvb(6,12,8,16)
vbc13o16=fvb(6,13,8,16)
vbo16o16=fvb(8,16,8,16)

vbo16o17=fvb(8,16,8,17)
vbo16o18=fvb(8,16,8,18)

kbc12c12=0.486
kbc12c13=0.573
kbc13c13=0.658
kbc13o16=0.423092

fk(e,mu)=sqrt(2*mu*931.5*e)/hbarc

fref(e,mu)=(31.29*36/sqrt(1000)/sqrt(e)+0.46/mu*e)-16.5

# xsec should be in barn
pene(e,mu,xsec)=xsec*1e2*fk(e,mu)**2/3.14

s(e,mu,xsec,zz)=log(1.0/pene(e,mu,xsec)-0)/sqrt(mu*zz)

# penetration through a Coulomb barrier with Radius R0 (R0>Rb)
peneCoulomb(e,mu,R0)=exp(87.21*sqrt(mu/6)/sqrt(e)+0.46*e)

sigma(e,smod)=smod/(e*exp(87.21/sqrt(e)+0.46*e))

# isotope offset for the normalized penetration function
#offset_cc(vb,zz)=9.4551*(vb/zz) - 1.5648
#offset_co(vb,zz)=4.7075*(vb/zz) - 0.7842

offset_cc(vb,zz)=9.6334*(vb/zz-fvb(6,12,6,12))-0.0099
offset_co(vb,zz)=offset_cc(vb,zz)

#offset_cc(vb,zz)=0
#offset_co(vb,zz)=0

# normalized penetration function
ffeff(e,vb)=e-(vb-fvb(6,12,6,12))
fe(e,zz)=(e/zz)**(-0.5)

pause -1 "Plot potential (carbon isotopes)"
plot [5:15] [0.1:0.17] "c12c12/verify.out" u 1:($4/36) w l lw 2, \
"c12c13/verify.out" u 1:($4/36) w l lw 2, \
"c13c13/verify.out" u 1:($4/36) w l lw 2, \
"c12c14/verify.out" u 1:($4/36) w l lw 2

pause -1 "Plot potential (carbon,oxygen isotopes)"
plot [5:15] [0.1:0.17] \
"c12c12/verify.out" u 1:($4/36) w l lw 2, \
"c12o16/verify.out" u 1:($4/48) w l lw 2, \
"o16o16/verify.out" u 1:($4/64) w l lw 2

pause -1 "Plot potential (B+C isotopes)"
plot [5:15] [0.1:0.17] \
"b10c12/verify.out" u 1:($4/30) w l lw 2, \
"b11c12/verify.out" u 1:($4/30) w l lw 2, \
"c12c12/verify.out" u 1:($4/36) w l lw 2

pause -1 "TEST normalized 12c13)unction(Sao-Paulo)"
ffeff(e,vb)=e-(vb-fvb(6,12,6,12))*0.2
fx(e,mu,xsec,zz,vb)=s(e,mu,xsec,zz)-(0.959935*fe(e,zz) -2.5843) - offset_cc(vb,zz)*0
plot [] [] \
"c12c13/cross.dat" u (fe(ffeff($1,vbc12c13),36)):(fx(ffeff($1,vbc12c13),mu_c12c13,$2,36,vbc12c13)) w l lw 2, \
"c12c12/cross.dat" u (fe(ffeff($1,vbc12c12),36)):(fx(ffeff($1,vbc12c12),mu_c12c12,$2,36,vbc12c12)) w l lw 2, \
"c13c13/cross.dat" u (fe(ffeff($1,vbc13c13),36)):(fx(ffeff($1,vbc13c13),mu_c13c13,$2,36,vbc13c13)) w l lw 2, \
"c12c14/cross.dat"  u (fe(ffeff($1,vbc12c14),36)):(fx(ffeff($1,vbc12c14),mu_c12c14,$2,36,vbc12c14)) w l lw 2,\
"c12c13/cross.dat" u (fe(ffeff($1,vbc12c13),36)):(fx(ffeff($1,vbc12c13),mu_c12c13,$2*1.2,36,vbc12c13)) w l lw 2,\
"c12c13/cross.dat" u (fe(ffeff($1,vbc12c13),36)):(fx(ffeff($1,vbc12c13),mu_c12c13,$2*0.8,36,vbc12c13)) w l lw 2

pause -1 "TEST normalized penetration function(M3Y)"
plot [] [] \
"henning/henning_c12c13_cc.dat" u (fe(ffeff($1,vbc12c13),36)):(fx(ffeff($1,vbc12c13),mu_c12c13,$2,36,vbc12c13)) w l lw 2, \
"henning/henning_c12c12_cc.dat" u (fe(ffeff($1,vbc12c12),36)):(fx(ffeff($1,vbc12c12),mu_c12c12,$2,36,vbc12c12)) w l lw 2, \
"henning/henning_c13c13_cc.dat" u (fe(ffeff($1,vbc13c13),36)):(fx(ffeff($1,vbc13c13),mu_c13c13,$2,36,vbc13c13)) w l lw 2, \
"henning/henning_c12c13_cc.dat" u (fe(ffeff($1,vbc12c13),36)):(fx(ffeff($1,vbc12c13),mu_c12c13,$2*1.2,36,vbc12c13)) w p lw 2,\
"henning/henning_c12c13_cc.dat" u (fe(ffeff($1,vbc12c13),36)):(fx(ffeff($1,vbc12c13),mu_c12c13,$2*0.8,36,vbc12c13)) w p lw 2

pause -1 "normalized penetration function(Sao-Paulo)"
fx(e,mu,xsec,zz,vb)=s(e,mu,xsec,zz)-(0.959935*fe(e,zz) -2.5843) - offset_cc(vb,zz)
plot [] [] \
"c12c13/cross.dat" u (fe($1,36)):(fx($1,mu_c12c13,$2,36,vbc12c13)) w l lw 2, \
"c12c12/cross.dat" u (fe($1,36)):(fx($1,mu_c12c12,$2,36,vbc12c12)) w lp lw 2, \
"c13c13/cross.dat" u (fe($1,36)):(fx($1,mu_c13c13,$2,36,vbc13c13)) w lp lw 2, \
"c12c14/cross.dat"  u (fe($1,36)):(fx($1,mu_c12c14,$2,36,vbc12c14)) w lp lw 2 

pause -1 "normalized penetration function(Sao-Paulo B+C isotopes)"
plot [] [] \
"c12c12/cross.dat" u (fe($1,36)):(fx($1,mu_c12c12,$2,36,vbc12c12)) w lp lw 2, \
"b10c12/cross.dat" u (fe($1,30)):(fx($1,mu_b10c12,$2,30,vbb10c12)) w l lw 5, \
"b11c12/cross.dat" u (fe($1,30)):(fx($1,mu_b11c12,$2,30,vbb11c12)) w l lw 5

pause -1 "normalized penetration function(M3Y)"
plot [] [] \
"henning/henning_c12c13_cc.dat" u (fe($1,36)):(fx($1,mu_c12c13,$2,36,vbc12c13)) w l lw 2, \
"henning/henning_c12c12_cc.dat" u (fe($1,36)):(fx($1,mu_c12c12,$2,36,vbc12c12)) w lp lw 2, \
"henning/henning_c13c13_cc.dat" u (fe($1,36)):(fx($1,mu_c13c13,$2,36,vbc13c13)) w lp lw 2,\
"henning/henning_c12c13_cc.dat" u (fe($1,36)):(fx($1,mu_c12c13,$2*1.15,36,vbc12c13)) w p lw 4,\
"henning/henning_c12c13_cc.dat" u (fe($1,36)):(fx($1,mu_c12c13,$2*0.85,36,vbc12c13)) w p lw 4

pause -1 "normalized penetration function(TDHF)"
plot [] [] \
"tdhf/x_12C_13C.out" u (fe($1,36)):(fx($1,mu_c12c13,$2,36,vbc12c13)) w l lw 2, \
"tdhf/x_12C_12C.out" u (fe($1,36)):(fx($1,mu_c12c12,$2,36,vbc12c12)) w lp lw 2, \
"tdhf/x_12C_14C.out" u (fe($1,36)):(fx($1,mu_c12c14,$2,36,vbc12c14)) w lp lw 2, \
"tdhf/x_12C_13C.out" u (fe($1,36)):(fx($1,mu_c12c13,$2*1.25,36,vbc12c13)) w p lw 4,\
"tdhf/x_12C_13C.out" u (fe($1,36)):(fx($1,mu_c12c13,$2*0.75,36,vbc12c13)) w p lw 4

pause -1 "plot C+C experimental data"
ffit(x)=a*x+b
fit ffit(x) "exp_data/c12c13/imp_experr.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) via a, b

fc13c13_exp_factor=1.2

fx(e,mu,xsec,zz,vb)=-(s(e,mu,xsec,zz)-(0.959935*fe(e,zz) -2.5843) - offset_cc(vb,zz))
plot [2.5:] [-54.8:] "exp_data/c12c12/becker.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
"exp_data/c12c12/spillane.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
"exp_data/c12c13/stockstad.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w lp lw 2,\
"exp_data/c12c13/imp_experr.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w lp lw 2,\
"exp_data/c13c13/trentalange.smod" u (fe($1,36)):(fx($1,mu_c13c13,sigma($1,$2)*fc13c13_exp_factor,36,vbc13c13)) w lp lw 2,\
"esw/c12c13_esw_dayras.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w l lw 4,\
"hindrance/c12o16.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w l lt 1 lw 2,\
"hindrance/o16o16.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2,64,vbo16o16)) w l lt 2 lw 2,\
"hindrance/c12c13.xsec" u (fe($1,36)):(fx($1,mu_c12c13,$2,36,vbc12c13)) w l lt 3 lw 2
#"esw/c13c13_esw.smod"  u (fe($1,36)):(fx($1,mu_c13c13,sigma($1,$2)*fc13c13_exp_factor*1.3,36,vbc13c13)) w l lw 4

###################
# C,O isotopes
##################
fx(e,mu,xsec,zz,vb)=-(s(e,mu,xsec,zz)-(0.959935*fe(e,zz) -2.5843) - offset_co(vb,zz))

pause -1 "normalized penetration function(Carbon,Oxygen isotopes)"
plot [] [] \
"c12c12/cross.dat" u (fe($1,36)):(fx($1,mu_c12c12,$2,36,vbc12c12)) w lp lw 1, \
"c12o16/cross.dat" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w l lw 1, \
"o16o16/cross.dat" u (fe($1,64)):(fx($1,mu_o16o16,$2,64,vbo16o16)) w l lw 1

pause -1 "plot experimental data for C,O fusion systems"
#plot [] [-54.8:] "exp_data/b10c12/stockstad_1977.xsec" u (fe($1,30)):(fx($1,mu_b10c12,$2,30,vbb10c12)) w lp lw 2,\
#"exp_data/b11c12/cujec_1977.xsec" u (fe($1,30)):(fx($1,mu_b11c12,$2/1000,30,vbb11c12)) w lp lw 2,\
#"esw/c12c13_esw_dayras.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w l lw 4,\
#"exp_data/c12c12/becker.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
#"exp_data/c12c12/spillane.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
#"exp_data/c12o16/ch77a.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w lp lw 2,\
#"exp_data/c13o16/dasmaha.xsec" u (fe($1,48)):(fx($1,mu_c13o16,$2/1000,48,vbc13o16)) w lp lw 2,\
#"exp_data/o16o16/wu84.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2/1000,64,vbo16o16)) w lp lw 2,\
#"henning/henning_o16o16_cc.dat" u (fe($1,64)):(fx($1,mu_o16o16,$2,64,vbo16o16)) w lp lw 2,\
#"hindrance/c12o16.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w l lw 2,\
#"hindrance/o16o16.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2,64,vbo16o16)) w l lw 2,\
#"hindrance/c12c13.xsec" u (fe($1,36)):(fx($1,mu_c12c13,$2,36,vbc12c13)) w l lw 2

plot [2.5:] [-54.8:] "esw/c12c13_esw_dayras.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w l lw 4,\
"exp_data/c12c12/becker.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
"exp_data/c12c12/spillane.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
"exp_data/c12o16/ch77a.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w lp lw 2,\
"exp_data/o16o16/wu84.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2/1000,64,vbo16o16)) w lp lw 2,\
"henning/henning_o16o16_20170307_sol2.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2/1000,64,vbo16o16)) w l lw 2,\
"hindrance/c12o16.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w l lw 2,\
"hindrance/o16o16.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2,64,vbo16o16)) w l lw 2,\
"hindrance/c12c13.xsec" u (fe($1,36)):(fx($1,mu_c12c13,$2,36,vbc12c13)) w l lw 2

pause -1 "plot experimental data for C,O fusion systems"
plot [2.5:] [-54.8:] "esw/c12c13_esw_dayras.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w l lw 4,\
"exp_data/c12c12/becker.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
"exp_data/c12c12/spillane.smod" u (fe($1,36)):(fx($1,mu_c12c12,sigma($1,$2),36,vbc12c12)) w lp lw 2,\
"exp_data/c12o16/ch77a.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w lp lw 2,\
"exp_data/o16o16/wu84.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2/1000,64,vbo16o16)) w lp lw 2,\
"tdhf/x_c_o_zzc.out" u (fe($1,48)):(fx($1,mu_c12o16,$2/1000,48,vbc12o16)) w lp lw 4,\
"tdhf/x_o_o_zzc.out" u (fe($1,64)):(fx($1,mu_o16o16,$2/1000,64,vbo16o16)) w lp lw 4,\
"henning/henning_o16o16_20170307_sol2.xsec" u (fe($1,64)):(fx($1,mu_o16o16,$2/1000,64,vbo16o16)) w l lw 4

pause -1 "c13+o16, c12+o16"
plot [2.5:] [-54.8:] "esw/c12c13_esw_dayras.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w l lw 4,\
"exp_data/c12o16/ch77a.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2,48,vbc12o16)) w lp lw 2,\
"exp_data/c12o16/dasmaha.xsec" u (fe($1,48)):(fx($1,mu_c12o16,$2/1000.0,48,vbc12o16)) w lp lw 2,\
"exp_data/c13o16/dasmaha.xsec" u (fe($1,48)):(fx($1,13.*16/29.,$2/1000.0,48,vbc13o16)) w lp lw 2

pause -1 "o+o"
plot [2.5:] [-54.8:] "esw/c12c13_esw_dayras.smod" u (fe($1,36)):(fx($1,mu_c12c13,sigma($1,$2),36,vbc12c13)) w l lw 4,\
"exp_data/o16o16/thomas.xsec" u (fe($1,64)):(fx($1,8.0,$2/1000,64,vbo16o16)) w lp lw 2,\
"exp_data/o16o17/thomas.xsec" u (fe($1,64)):(fx($1,16*17/33.,$2/1000,64,vbo16o17)) w lp lw 2,\
"exp_data/o16o18/thomas.xsec" u (fe($1,64)):(fx($1,16*18/34.,$2/1000,64,vbo16o18)) w lp lw 2

print "c12c12 Coulomb barrier height : ", fvb(6,12,6,12), fvb(6,12,6,13), fvb(6,13,6,13)
print fvb(6,12,6,12), fvb(6,12,8,16), fvb(8,16,8,16)

print fe(1.0,64), fx(1.0,mu_o16o16,4.18e-52,64,vbo16o16)
