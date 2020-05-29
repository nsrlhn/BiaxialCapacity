# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 11:12:53 2020

@author: ensar
"""
import math,numpy
import matplotlib.pyplot as plt

def sind(teta):
    return math.sin(math.radians(teta))
def cosd(teta):
    return math.cos(math.radians(teta))
def tand(teta):
    return math.tan(math.radians(teta))
def cotd(teta):
    return 1/tand(teta)

def Steel_force(h,b,ky,teta,xs,ys,As,e_cu,c):
    fyk=365         #in MPa
    Es= 200000      #in MPa
    e_sy= fyk/Es    #strain where steel yields
    e_su= 0.1       #strain where steel fails
    
    P=0
    Fs=0
    Msy=0
    Msx=0
    for i1 in range(len(ys)):
        P=-xs[i1]*sind(teta)+ys[i1]*cosd(teta)-b/2*sind(teta)+(ky*h-h/2)*cosd(teta)   #new coordinate system
        e_s=e_cu/c*P
        if abs(e_s) <= e_sy:
            sigma_s = Es*e_s
        elif e_sy < abs(e_s) and abs(e_s) <= e_su:
            sigma_s = fyk*e_s/abs(e_s)
        else:
            sigma_s = 0
        F=sigma_s*As[i1]
        Msy=Msy+F*xs[i1]
        Msx=Msx+F*ys[i1]
        Fs=Fs+F
    return Fs,Msx,Msy

def Shaded_area(h,b,ky,teta,k1):
    x=k1*ky*h
    if x < h:
        if b*tand(teta) < x:
            Ac=b*x-b**2*tand(teta)/2
            k=x-b*tand(teta)
            yc=h/2-(b*k*k/2+b/2*(x-k)*((x-k)/3+k))/Ac
            xc=(b*k*b/2+(x-k)*b/2*b/3)/Ac-b/2
            #Case=1
        else:
            Ac=(x)**2*cotd(teta)/2
            yc=h/2-x/3
            xc=x*cotd(teta)/3-b/2
            #Case=2
    else:
        if b*tand(teta) < x:
            k=x-b*tand(teta)
            if k > h:
                Ac=b*h
                yc=0
                xc=0
            else:
                Ac=b*h-(h-k)**2*cotd(teta)/2
                yc=(k*(h-k)*cotd(teta)*(h/2-k/2)+(2/3*(h-k)-h/2)*(h-k)**2*cotd(teta)/2)/Ac
                x2=b-(h-k)*cotd(teta)
                xc=((h-k)*x2*(x2/2-b/2)+(h-k)**2*cotd(teta)/2*(b/2-2/3*(h-k)*cotd(teta)))/Ac
            #Case=3
        else:
            Ac=h**2*(2*x/h-1)*cotd(teta)/2
            yc=((h/2-h/3)*h**2/2*cotd(teta))/Ac
            xc=(h*(x-h)*cotd(teta)*((x-h)*cotd(teta)-b)/2+h**2/2*cotd(teta)*(h*cotd(teta)/3+(x-h)*cotd(teta)-b/2))/Ac;
            #Case=4
    return Ac,xc,yc


#Section and Material Properties     
fck=16          #in Mpa
k1=0.85
e_cu= 0.003     #strain where concrete fails
h=500           #width of the column in mm
b=500           #depth of the column in mm

xs=[0, 0, 0, 0]       #location of steel members with respect to center of gravity in mm
ys=[220, 73.3, -73.3, -220]      #location of steel members with respect to center of gravity in mm
As=[2123, 1061, 1061, 2123]      #area of the steel in mm2

N_max=0.85*fck*b*h+sum(As)*365
N_min=0
loop1=50
N = numpy.arange(N_min, N_max, (N_max-N_min)/(loop1-1))
Mx = []
My = []
N_plot = []
for i1 in range(loop1-1):
    error=0.000001
    if N[i1]*error > 1:
        errorlimit = N[i1]*error
    else:
        errorlimit = 1
    teta_min=0
    teta_mcleax=85
    loop2=86
    teta=teta_min
    Mx.append([None]*loop2)
    My.append([None]*loop2)
    N_plot.append([None]*loop2)
    for i2 in range(loop2):
        if teta < 10:
            ky_max=1
        else:
            ky_max=4
        ky_min=0
        ky=(ky_max+ky_min)/2
        c = ky*h*cosd(teta)
        Fs,Msx,Msy = Steel_force(h,b,ky,teta,xs,ys,As,e_cu,c)
        Ac,xc,yc = Shaded_area(h,b,ky,teta,k1)
        Fc = 0.85*fck*Ac
        error=N[i1]-Fc-Fs
        j=0
        while abs(error)>errorlimit:
            if error<0:
                ky_max=ky
                ky=0.5*(ky_min+ky)
            else:
                ky_min=ky
                ky=0.5*(ky+ky_max)
            c = ky*h*cosd(teta)
            Fs,Msx,Msy = Steel_force(h,b,ky,teta,xs,ys,As,e_cu,c)
            Ac,xc,yc = Shaded_area(h,b,ky,teta,k1)
            Fc = 0.85*fck*Ac
            error=N[i1]-Fc-Fs #sum(Fs) !!!!!!!!!!!!
            j=j+1;
            if j > 1000:
                break
        if ky > 3.9:
            Mx[i1][i2]= 0
            My[i1][i2]= 0
            N_plot[i1][i2]= 0
        else:
            if j > 999:
                Mx[i1][i2]= 0
                My[i1][i2]= 0
                N_plot[i1][i2]= 0
            else:
                Mx[i1][i2]=(Fc*yc+Msx)/1000000 #kNm
                My[i1][i2]=(Fc*xc+Msy)/1000000 #kNm
                N_plot[i1][i2]=N[i1]/1000      #kN
        teta=teta+1

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
z = numpy.array(N_plot)
ax.plot_wireframe(Mx, My, z, rstride=10, cstride=10)

plt.show()
"""
surf(Mx,My,N_plot);
xlabel 'Mx (kNm)'
ylabel 'My (kNm)'
zlabel 'N (kN)'
xL = xlim;
yL = ylim;
zL = zlim;
line([0 0], yL);
line(xL, [0 0]);
line([0,0],[0,0],zL);
"""