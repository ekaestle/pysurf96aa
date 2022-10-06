import numpy as np
from pysurf96aa import surf96aa,surf96,layermod2depthmod
import matplotlib.pyplot as plt

def makemod(): # constant thickness layermodel
    thickness = np.ones(80)*2
    thickness[-1] = 0.
    vs = np.zeros(len(thickness))
    vs[:4] = 2.5
    vs[4:8] = 3.2
    vs[8:12] = 3.5
    vs[12:16] = 3.9
    vs[16:28] = 4.4
    vs[28:70] = 4.3
    vs[70:] = 4.8
    psi2amp = np.ones_like(vs)*0.05
    psi2amp[40:] = 0.01
    psi2dir = np.ones_like(vs)*-30/180.*np.pi
    psi2dir[:4] = 30/180.*np.pi
    psi2dir[12:] = 60/180.*np.pi
    vp = vs * 1.73
    rho = vp * .32 + .77
    return thickness,vp,vs,rho,psi2amp,psi2dir

# Define the velocity model in km and km/s
thickness,vp,vs,rho,psi2amp,psi2dir = makemod()
depth,vpm,vsm,rhom,psi2ampm,psi2dirm = layermod2depthmod(
    thickness,(vp,vs,rho,psi2amp,psi2dir))

# Plot the model
plt.figure(figsize=(12,6))
ax1 = plt.subplot(151)
plt.plot(vpm,depth)
ax1.set_ylim(np.max(depth),0)
plt.xlabel("Vp [km/s]")
plt.ylabel("Depth [km]")
ax2 = plt.subplot(152,sharey=ax1)
ax2.plot(vsm,depth)
plt.xlabel("Vs [km/s]")
ax3 = plt.subplot(153,sharey=ax1)
plt.plot(rhom,depth)
plt.xlabel("Density [g/cm³]")
ax4 = plt.subplot(154,sharey=ax1)
ax4.plot(psi2ampm,depth)
plt.xlabel("Anisotropy amplitude")
ax5 = plt.subplot(155,sharey=ax1)
ax5.plot(psi2dirm/np.pi*180.,depth,'o')
plt.xlabel("Anisotropy fastdir [deg]")
plt.show()

# Periods we are interested in
periods = np.around(np.logspace(np.log10(2),np.log10(80),8))

# split layers into sub-layers (not necessary, the sensitivity kernels look
# better, but  the resulting anisotropy parameters should not change)
nrefine = 1 

v1 = surf96(thickness, vp, vs, rho, periods,
            wave='rayleigh', mode=1, velocity='phase', flat_earth=False)
velocities,c_aa_amp,c_aa_ang,Lsen,dcrda,dcrdl = surf96aa(
    thickness, vp, vs, rho, psi2amp, psi2dir, periods, 
    nrefine=nrefine, wave='rayleigh', mode=1, velocity='phase', 
    flat_earth=False, return_sensitivities=True)

# if nrefine = 1, thickness_refine is identical to thickness
thickness_refine = np.append(np.repeat(thickness[:-1],nrefine)/nrefine,0)

# plot of the layerwise sensitivities
plt.figure(figsize=(12,6))
ax1 = plt.subplot(141)
plt.plot(vsm,depth)
plt.xlabel('Vs')
plt.ylabel('Depth [km]')
plt.ylim(np.max(depth),0)
plt.subplot(142,sharey=ax1)
plt.title("Layerwise sensitivity (has to be divided by layerthickness to "+
          "get the continuous sensitivity kernels)",loc='left',fontsize=8)
for i,line in enumerate(dcrdl):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.ylim(np.max(depth_refine),0)
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$dC/dL$')
plt.subplot(143,sharey=ax1)
for i,line in enumerate(dcrda):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$dC/dA$')
plt.subplot(144,sharey=ax1)
for i,line in enumerate(Lsen):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$Lsen$')
plt.show()

# plot of resulting phase velocity and anisotropic parameters
fig = plt.figure(figsize=(8,7))
ax1 = fig.add_subplot(311)
ax1.plot(periods,velocities)
ax1.set_ylabel("Phase velocity [km/s]")
ax2 = fig.add_subplot(312,sharex=ax1)
ax2.plot(periods,c_aa_amp)
ax2.set_ylim(0,np.max(c_aa_amp)*1.1)
ax2.set_ylabel("Anisotropic amplitude")
ax3 = fig.add_subplot(313,sharex=ax1)
ax3.plot(periods,c_aa_ang/np.pi*180,'o')
ax3.set_ylabel("Anisotropic angles [deg]")
ax3.set_ylim(-90,90)
ax3.set_xlabel("Period [s]")
plt.show()

#%% Test 2
# use the same model, but with thicker layers. The layers are refined using the
# nrefine parameter. The result is the same as above
thickness,vp,vs,rho,psi2amp,psi2dir = makemod()
thickness = thickness[::4]*4
thickness[-1] = 0. # halfspace thickness
depth = np.cumsum(thickness)-thickness/2.
vs = vs[::4]
vp = vp[::4]
rho = rho[::4]
psi2amp = psi2amp[::4]
psi2dir = psi2dir[::4]

depth,vpm,vsm,rhom,psi2ampm,psi2dirm = layermod2depthmod(
    thickness,(vp,vs,rho,psi2amp,psi2dir))

nrefine=4 # refine layers by splitting each layer in 4 sublayers
velocities,c_aa_amp,c_aa_ang,Lsen,dcrda,dcrdl = surf96aa(
    thickness, vp, vs, rho, psi2amp, psi2dir, periods,
    nrefine=nrefine, wave='rayleigh', mode=1, velocity='phase', 
    flat_earth=False, return_sensitivities=True)

thickness_refine = np.append(np.repeat(thickness[:-1],nrefine)/nrefine,0)

# plot of the layerwise sensitivities
plt.figure(figsize=(12,6))
ax1 = plt.subplot(141)
plt.plot(vsm,depth)
plt.xlabel('Vs')
plt.ylabel('Depth [km]')
plt.ylim(np.max(depth),0)
plt.subplot(142,sharey=ax1)
plt.title("Layerwise sensitivity (has to be divided by layerthickness to "+
          "get the continuous sensitivity kernels)",loc='left',fontsize=8)
for i,line in enumerate(dcrdl):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.ylim(np.max(depth_refine),0)
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$dC/dL$')
plt.subplot(143,sharey=ax1)
for i,line in enumerate(dcrda):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$dC/dA$')
plt.subplot(144,sharey=ax1)
for i,line in enumerate(Lsen):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$Lsen$')
plt.show()

# plot of resulting phase velocity and anisotropic parameters
fig = plt.figure(figsize=(8,7))
ax1 = fig.add_subplot(311)
ax1.plot(periods,velocities)
ax1.set_ylabel("Phase velocity [km/s]")
ax2 = fig.add_subplot(312,sharex=ax1)
ax2.plot(periods,c_aa_amp)
ax2.set_ylim(0,np.max(c_aa_amp)*1.1)
ax2.set_ylabel("Anisotropic amplitude")
ax3 = fig.add_subplot(313,sharex=ax1)
ax3.plot(periods,c_aa_ang/np.pi*180,'o')
ax3.set_ylabel("Anisotropic angles [deg]")
ax3.set_ylim(-90,90)
ax3.set_xlabel("Period [s]")
plt.show()

#%% Test 3
# this is again the same model as above, but the layers have varying thicknesses
layers = []
thickness,vp,vs,rho,psi2amp,psi2dir = makemod()
reftest = 0.
for i in range(len(thickness)):
    test = thickness[i]+vp[i]+vs[i]+rho[i]+psi2amp[i]+psi2dir[i]
    if test!=reftest:
        layers.append([thickness[i],vp[i],vs[i],rho[i],psi2amp[i],psi2dir[i]])
    else:
        layers[-1][0] += thickness[i]
    reftest = test
layers = np.array(layers)
thickness = layers[:,0]
#thickness[-1] = 0.
vp = layers[:,1]
vs = layers[:,2]
rho= layers[:,3]
psi2amp = layers[:,4]
psi2dir = layers[:,5]

depth,vpm,vsm,rhom,psi2ampm,psi2dirm = layermod2depthmod(
    thickness,(vp,vs,rho,psi2amp,psi2dir))


nrefine=1 # nrefine=1 : do not refine layers
v2 = surf96(thickness, vp, vs, rho, periods,
            wave='rayleigh', mode=1, velocity='phase', flat_earth=False)
velocities,c_aa_amp,c_aa_ang,Lsen,dcrda,dcrdl = surf96aa(
    thickness, vp, vs, rho, psi2amp, psi2dir, periods,
    nrefine=nrefine, wave='rayleigh', mode=1, velocity='phase', 
    flat_earth=False, return_sensitivities=True)

thickness_refine = np.append(np.repeat(thickness[:-1],nrefine)/nrefine,0)

# plot of the layerwise sensitivities
plt.figure(figsize=(12,6))
ax1 = plt.subplot(141)
plt.plot(vsm,depth)
plt.xlabel('Vs')
plt.ylabel('Depth [km]')
plt.ylim(np.max(depth),0)
plt.subplot(142,sharey=ax1)
plt.title("Layerwise sensitivity (has to be divided by layerthickness to "+
          "get the continuous sensitivity kernels)",loc='left',fontsize=8)
for i,line in enumerate(dcrdl):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.ylim(np.max(depth_refine),0)
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$dC/dL$')
plt.subplot(143,sharey=ax1)
for i,line in enumerate(dcrda):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$dC/dA$')
plt.subplot(144,sharey=ax1)
for i,line in enumerate(Lsen):
    depth_refine,layersens = layermod2depthmod(thickness_refine,(np.append(line,0),))
    plt.plot(layersens,depth_refine,label="%.1f" %periods[i])
plt.legend(loc='lower right',title='Periods')
plt.xlabel('$Lsen$')
plt.show()

# plot of resulting phase velocity and anisotropic parameters
fig = plt.figure(figsize=(8,7))
ax1 = fig.add_subplot(311)
ax1.plot(periods,velocities)
ax1.set_ylabel("Phase velocity [km/s]")
ax2 = fig.add_subplot(312,sharex=ax1)
ax2.plot(periods,c_aa_amp)
ax2.set_ylim(0,np.max(c_aa_amp)*1.1)
ax2.set_ylabel("Anisotropic amplitude")
ax3 = fig.add_subplot(313,sharex=ax1)
ax3.plot(periods,c_aa_ang/np.pi*180,'o')
ax3.set_ylabel("Anisotropic angles [deg]")
ax3.set_ylim(-90,90)
ax3.set_xlabel("Period [s]")
plt.show()
#%% compare velocities from fine- and rough-layered models
plt.figure()
plt.plot(periods,v1,label='model with fine layers')
plt.plot(periods,v2,'--',label='model with few layers')
plt.xlabel("period [s]")
plt.ylabel("phase velocity [km/s]")
plt.legend()
plt.show()

#%% unused
# from matplotlib.collections import LineCollection

# fig = plt.figure(figsize=(12,10))
# ax0 = fig.add_subplot(411)
# ax0.plot(periods,velocities)
# ax = fig.add_subplot(412)
# lines = []
# for i,d in enumerate(propdirs):
#     lines.append(np.column_stack((periods,dc[:,i])))
# linecoll = LineCollection(lines, array=propdirs/np.pi*180,
#                           cmap=plt.cm.rainbow, linewidths=1)
# ax.add_collection(linecoll)
# ax.plot(periods,np.zeros_like(periods),linewidth=0.1,color='k')
# #plt.legend()
# axcb = plt.colorbar(linecoll,shrink=0.8,pad=0.02,fraction=0.08)
# ax2 = fig.add_subplot(413)
# ax2.plot(periods,c_aa_amp)
# ax2.set_ylim(0,np.max(c_aa_amp)*1.1)
# ax3 = fig.add_subplot(414)
# ax3.plot(periods,c_aa_ang/np.pi*180)
# ax3.set_ylim(-90,90)
# plt.show()



