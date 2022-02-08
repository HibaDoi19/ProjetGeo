from math import *
import numpy as np

def geocentrique2geodesique(x,y,z,a,b):
    x=float(x)
    y=float(y)
    z=float(z)
    a=float(a)
    b=float(b)



    e=sqrt(1-b**2/a**2)




    if x ==0 and y==0 and z>0 :
        return 90,0,z-a/sqrt(1-e**2*sin(pi/2)**2)*(1-e**2)
    elif x ==0 and y==0 and z<0:
        return -90,0,z-a/sqrt(1-e**2*sin(-pi/2)**2)*(1-e**2)
    elif x ==0 and z==0 and y>0:
        return 0,90,sqrt(y**2)-a
    elif x ==0 and z==0 and y<0:
        return 0,-90,sqrt(y**2)-a

    else:

        if x==0 :
            g=float("inf")
        else:
            g=y/x

        long=np.arctan(g)
        phi0=np.arctan(z/(x**2+y**2)*(1/1-e**2))
        N0=a/sqrt(1-e**2*sin(phi0)**2)
        if z!=0 :

            phi1=np.arctan(z/sqrt(x**2+y**2)*(1+(N0*e**2*sin(phi0))/z))
            while abs(phi1-phi0)>=0.001 :
                phi0=phi1
                N1=a/sqrt(1-e**2*sin(phi1)**2)
                phi1=np.arctan(z/sqrt(x**2+y**2)*(1+(N1*e**2*sin(phi0))/z))
        else:
            phi1=0
            N1=a


        if phi1 !=pi/2:
            h=sqrt(x**2+y**2)/cos(phi1)-N1
        if phi1 !=0:
            h=z/sin(phi1)-N1*(1-e**2)
        if x>=0 :
            long=np.arctan(g)*180/pi
        elif y>0 :
            long=np.arctan(g)*180/pi+180
        else :
            long=np.arctan(y/x)*180/pi-180


    return phi1*180/pi,long,h