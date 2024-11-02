import numpy as np
import matplotlib.pyplot as plt
import struct
#-----------------------
from matplotlib import animation     # animation function
from IPython.display import HTML     # this displays the animation within the notebook window
from matplotlib import cm

#-----------------------
#-----------------------

N = 256         # Resolution
Nframes = 100     # Number of frames/snapshots
comp = 'zet'    # Component; zet, eta, xii

def readFile(file_name):
    Zin = np.zeros((N,N,N))
    with open(file_name, 'rb') as inh:
        indata = inh.read()

    # There may be more Pythonic way to do this?
    for iz in range(N):
        for iy in range(N):
            for ix in range(N):
                ik = ((iz*N + iy)*N+ix)*8
                Z = struct.unpack("d", indata[ik:ik+8])
                Zin[ix,iy,iz] = float(Z[0])
    return Zin


x = np.linspace(0,2*np.pi,N)
y = np.linspace(0,2*np.pi,N)
[xx,yy] = np.meshgrid(x,y)

fig3,ax3 = plt.subplots()            # this needs to be empty subplots for animating
fig3.set_size_inches(4, 5, True)
#--------------------------------------
# We define an animation function. 
# This is called sequentially to draw
# each frame of the animation
#--------------------------------------
def animate(i):
    
    file_name = comp+f"{(i+1):03d}"+'.raw'
    Z = readFile(file_name)
    ax3.clear()
    ax3.set_aspect(1)
    ax3.set_xlabel('$x$')
    ax3.set_ylabel('$y$')
    ax3.pcolor(xx,yy,Z[:,:,N//2],cmap=cm.twilight)
    ax3.set_title('$i=$%2d'%(i)) 

anim = animation.FuncAnimation(fig3, animate,
                               frames=Nframes-1, interval=50, blit=False)

FFwriter = animation.FFMpegWriter(fps=10)

anim.save('animation.mp4', writer = FFwriter)
