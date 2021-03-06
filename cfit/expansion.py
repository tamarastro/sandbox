from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
plt.close("all")

def ezinv(z,om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0):
    ok = 1.-om-ox-orr
    a  = 1./(1.+z)
    f = a**(3*(1.+w0+wa)) * np.exp(3*wa*(1-a))
    ez = np.sqrt( orr/a**4 + om/a**3 + ok/a**2 + ox/f )
    return 1./ez

def xx(z,om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0):
    xx = integrate.quad(ezinv,0.,z,args=(om,ox,w0,wa,orr))
    # Multiply by R_0 (i.e. c/H0 if flat or c/H0*sqrt(abs(1/ok)) to get in units of Mpc.
    # Remember to curvature correct if want proper distance.
    return xx

def vv(xx,c=299792.): # c in km/s
    return xx*c

def vvapprox(z,q=-0.55,j=1.0,c=299792.):
    return c*z/(1+z)*( 1+0.5*(1-q)*z - (1./6)*(1-q-3*q**2+j)*z**2 )

def dist_curve_correct(xx, ok):
    # Corrects the comoving distance, x, for curvature, ok.
    # Result is proper distance / (c/H0)
    if ok < 0.0:
        dk = np.sin(np.sqrt(-ok)*xx)/np.sqrt(-ok)
    elif ok > 0.0:
        dk = np.sinh(np.sqrt(ok)*xx)/np.sqrt(ok)
    else:
        dk = xx
    return dk

def dist_lum(xx,ok,z):
    return dist_curve_correct(xx,ok)*(1+z)

def H0kmsmpc2Gyr(H0kmsmpc):
    # Convert to inverse seconds
    H0s = H0kmsmpc * 3.24e-20 #s-1
    # Convert to inverse Giga-years
    H0Gyr = H0s* 3.156e16 #Gyr-1
    return H0Gyr

def adotinv(a,om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0,H0kmsmpc=70.):
    # Calculates 1/(a*dot(a)).
    ok = 1.-om-ox-orr
    H0Gyr=H0kmsmpc2Gyr(H0kmsmpc)
    adot = H0Gyr * a * np.sqrt ( orr/a**4 + om/a**3 + ok/a**2 + ox/( a**(3*(1+w0+wa))*np.exp(3*wa*(1-a)) ) )
    return 1./adot

def time(a0, a1, om=0.3,ox=0.7,w0=-1.0,wa=0.0,orr=0.0,H0kmsmpc=70.):
    # Calculates the time between two scalefactors.
    time = integrate.quad(adotinv,a0,a1,args=(om,ox,w0,wa,orr,H0kmsmpc))
    return time



#################
n=51
zmax=0.1
zs=np.arange(n)/(n-1)*zmax
xs=np.zeros(n)
xslo=np.zeros(n)
xshi=np.zeros(n)
xslow=np.zeros(n)
xshiw=np.zeros(n)
err=np.zeros(n)


#xxs = np.vectorize(xx)
for i in np.arange(n):
    xs[i],err[i]=xx(zs[i])
    xslo[i],err[i]=xx(zs[i],om=0.29,ox=0.71)
    xshi[i],err[i]=xx(zs[i],om=0.31,ox=0.69)
    xslow[i],err[i]=xx(zs[i],w0=-0.9)
    xshiw[i],err[i]=xx(zs[i],w0=-1.1)


vs = vv(xs)
vslo = vv(xslo)
vshi = vv(xshi)
vslow = vv(xslow)
vshiw = vv(xshiw)

vapproxs = vvapprox(zs)
DL=dist_lum(xs,0.0,zs)

for i in np.arange(n): print(zs[i], vs[i],vapproxs[i],vapproxs[i]-vs[i],vslo[i],vs[i]-vslo[i])

# Plot
f = plt.figure()
plt.subplots_adjust(hspace=0.001)

ax1 = plt.subplot(211)
plt.ylabel('Velocity (km/s)')
plt.xlim(0.0,0.1)

ax1.plot(zs, vs, label='om=0.3,w=-1')
ax1.plot(zs, vslo, label='om=0.29')
ax1.plot(zs, vshi, label='om=0.31')
ax1.plot(zs, vslow, label='w=-0.9')
ax1.plot(zs, vshiw, label='w=-1.1')
ax1.plot(zs, vapproxs, '--',label='approx')
ax1.legend(frameon=False,loc=2)
#plt.yticks(np.arange(-0.9, 1.0, 0.4))
#plt.ylim(-1, 1)

ax2 = plt.subplot(212, sharex=ax1)
ax2.plot(zs, vs-vapproxs, label='vs')
ax2.plot(zs, vslo-vapproxs, label='om=0.29')
ax2.plot(zs, vshi-vapproxs, label='om=0.31')
ax2.plot(zs, vslow-vapproxs, label='w=-0.9')
ax2.plot(zs, vshiw-vapproxs, label='w=1.1')
ax2.plot(zs, vapproxs-vapproxs, '--', label='approx')
#plt.yticks(np.arange(0.1, 1.0, 0.2))
#plt.ylim(0, 1)
plt.xlim(0.0,0.1)

plt.xlabel('Redshift')
plt.ylabel('Delta v')
xticklabels = ax1.get_xticklabels() #+ ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)

plt.savefig('Velocities.png',bbox_inches='tight')
#blah=input('check1:')
plt.show()
plt.close()
#plt.close('all')


#####Second Plot######
omarr = np.arange(0.29,0.31,0.01)
w0arr = np.arange(-1.1,-0.9,0.2)
print(w0arr)
f2=plt.figure()
plt.subplots_adjust(hspace=0.001)

ax1 = plt.subplot(211)
plt.ylabel('Velocity (km/s)')
xticklabels = ax1.get_xticklabels() #+ ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)

ax2 = plt.subplot(212, sharex=ax1)
#plt.xlabel('Redshift')
plt.ylabel('Delta v')
plt.xlim(0,0.1)

#ax3 = plt.subplot(213, sharex=ax1)
#plt.xlabel('Redshift')
#plt.ylabel('Luminosity Distance')
#plt.xlim(0,0.1)

ax1.plot(zs, vapproxs, '--', label='approx')
ax2.plot(zs, vapproxs-vapproxs,'--', label='approx')

for om in omarr:
    for i in np.arange(n):
        xs[i],err[i]=xx(zs[i],om=om,ox=1.0-om)
    vs=vv(xs)
    ax1.plot(zs, vs, label='om='+str(om))
    ax2.plot(zs, vs-vapproxs, label='om='+str(om))
for w0 in w0arr:
    for i in np.arange(n):
        xs[i],err[i]=xx(zs[i],w0=w0)
        vs=vv(xs)
    ax1.plot(zs, vs, label='w='+str(w0))
    ax2.plot(zs, vs-vapproxs, label='w0='+str(w0))
ax1.legend(frameon=False,loc=2)

#plt.show()

#plt.plot(zs,vs)
#plt.plot(zs,vapproxs)
#plt.show()

