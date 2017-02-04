from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

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


#################
n=50
zmax=0.1
zs=np.arange(n+1)/(n)*zmax # set array of redshifts
zs=zs[1:]                  # get rid of z=0
xs=np.zeros(n)             # set array of comoving distances
err=np.zeros(n)            # set array of errors in integral

# Set constants
c=299792.   # km/s
H0fid = 70. # km/s/Mpc

# Calculate approximate velocity
vapproxs = vvapprox(zs)

# Calculate fiducial model velocity
for i in np.arange(n):
    xs[i],err[i]=xx(zs[i])
vs = vv(xs)

# Calculate Luminosity distance in fiducial model
DL=c/H0fid*dist_lum(xs,0.0,zs) #= dist_curve_correct(xs,1.0)*(1+zs)

H0 = vs *(1+zs) / DL      # Array of H0s that would be "measured" if the fiducial model was the correct one
print('H0 fiducial model=',np.mean(H0))


###################################
# Set up a range of possible models
omarr = np.arange(0.29,0.31,0.005)
w0arr = np.arange(-1.1,-0.9,0.05)

###################################
# Calculate Luminosity distance in the range of possible models, and calculating  H_0 values inferred
f=plt.figure()


##### Loop over omega_m models #####
ax1 = plt.subplot(211)
plt.xlabel('om')
plt.ylabel('H0')
ax1.get_yaxis().get_major_formatter().set_useOffset(False)
ax1.get_yaxis().get_major_formatter().set_scientific(False)

H0arrapprox=np.zeros(len(omarr))
H0arr      =np.zeros(len(omarr))
for j in range(len(omarr)):
    for i in np.arange(n):
        xs[i],err[i]=xx(zs[i],om=omarr[j],ox=1.0-omarr[j])
    DL=c/H0fid*dist_lum(xs,0.0,zs)
    H0s = vs * (1+zs) / DL
    H0arr[j] = np.mean(H0s)
    H0s = vapproxs * (1+zs) / DL
    H0arrapprox[j] = np.mean(H0s)
ax1.plot(omarr, H0arr, label='exact')
ax1.plot(omarr, H0arrapprox, '--', label='approx')
ax1.plot([min(omarr),max(omarr)], [H0fid,H0fid],'--')
plt.xlim([min(omarr),max(omarr)])
ax1.annotate('H0 inferred when true model is om, but assumed model is om=0.3',xy=(min(omarr)+(max(omarr)-min(omarr))*0.02,max(H0arr)-(max(H0arr)-min(H0arr))*0.05))


##### Loop over w0 models #####
ax2 = plt.subplot(212)
plt.xlabel('w0')
plt.ylabel('H0')

H0arrapprox=np.zeros(len(w0arr))
H0arr      =np.zeros(len(w0arr))
for j in range(len(w0arr)):
    for i in np.arange(n):
        xs[i],err[i]=xx(zs[i],w0=w0arr[j])
    DL=c/H0fid*dist_lum(xs,0.0,zs)
    H0s = vs * (1+zs) / DL
    H0arr[j] = np.mean(H0s)
    H0s = vapproxs * (1+zs) / DL
    H0arrapprox[j] = np.mean(H0s)
ax2.plot(w0arr, H0arr, label='exact')
ax2.plot(w0arr, H0arrapprox, '--', label='approx')
ax2.plot([min(w0arr),max(w0arr)], [H0fid,H0fid],'--')
ax2.annotate('H0 inferred when true model is w0, but assumed model is w0=-1.0',xy=(min(w0arr)+(max(w0arr)-min(w0arr))*0.02,max(H0arr)-(max(H0arr)-min(H0arr))*0.05))
plt.xlim([min(w0arr),max(w0arr)])

# Save figure
plt.savefig('H0test.png',bbox_inches='tight')
plt.show()
plt.close(f)
