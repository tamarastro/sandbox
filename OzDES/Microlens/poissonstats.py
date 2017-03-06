import numpy as np
import matplotlib.pyplot as plt

def poisson(r,k):
    return r**k * np.exp(-r) / np.math.factorial(k)

k=0
r=4
print(poisson(r,k))
print(1-poisson(r,k))
exit()

N=range(0,16)
print(N)
count = [57,203,383,525,532,408,273,139,45,27,10,4,0,1,1,0]
#plt.plot(N,count)
#plt.show()
tot = np.sum(count)
rate = 4.0

rate = 25
n=16
karr = np.arange(0,n)
#print(karr)

poissonarr = np.empty(shape=n)
cumularr = np.empty(shape=n)
#print(len(poissonarr))
#print(len(karr))
i=0
for k in karr:
    poissonarr[i]=poisson(rate,k)
    cumularr[i]=np.sum(poissonarr[0:i+1])
    i+=1


probge4 = np.sum(poissonarr[4:11])
proble3 = np.sum(poissonarr[0:4])

f=plt.figure()
ax1 = plt.subplot(211)

ax1.plot(karr,poissonarr)#,label=str(rate))
ax1.plot(N,count/tot,'o')

#ax1.plot(karr,cumularr)#,label=str(rate))
plt.xlabel('k')
plt.ylabel('Probability')
#plt.ylim(0,0.2)


ax2 = plt.subplot(212)

ax2.plot(karr,cumularr)#,label=str(rate))
plt.xlabel('k')
plt.ylabel('Cumulative')

plt.show()
