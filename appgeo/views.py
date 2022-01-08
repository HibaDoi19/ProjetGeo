from typing import Text
from django.shortcuts import render
from django.http import HttpResponse
from math import *
from plotly.offline import iplot, init_notebook_mode
import plotly.graph_objects as go
from numpy import sin, cos, pi
import  numpy as np
import plotly.express as px 
import datetime
import glob
import logging
import os


def direct(latitude1, longitude1, alpha1, s, a, b):
                if a<0 or b<0:
                        texte="a et b doivent être positifs"
                        return texte
                else:
                        f = (a-b)/a
                        ep = sqrt((a**2-b**2)/b**2)
                        #Calcul de la latitude réduite beta1
                        TanBeta1 = (1-f) * tan(latitude1)
                        Beta1 = atan( TanBeta1 )
                        #Calcul de la latitude réduite beta0
                        CosBeta0 = cos(Beta1)*sin(alpha1)
                        Beta0 = acos(CosBeta0)
                        #Calcul de la constante w2
                        w2 = (sin(Beta0)*ep)*2
                        #Calcul de la distance angulaire sigma1 sur la sphère auxiliaire
                        sigma1 = atan2( TanBeta1, cos(alpha1) )
                        #Calcul de l'azimut de la géodésique à l'équateur alphae
                        Sinalphae = CosBeta0
                        alphae = asin(Sinalphae)
                        cosalphae_2= 1.0 - Sinalphae**2
                        #Calcul des constantes de Vincenty A' et B'
                        A = 1.0 + (w2 / 16384) * (4096 + w2 * (-768 + w2 * (320 - 175 * w2) ) )
                        B = (w2 / 1024) * (256 + w2 * (-128 + w2 * (74 - 47 * w2) ) )
                        #Calcul de la distance angulaire sigma sur le grand cercle de la sphère auxiliaire
                        sigma0 = (s / (b * A))
                        sigma_m2 = (2 * sigma1 + sigma0)
                        delta_sigma = B * sin(sigma0) * ( cos(sigma_m2)+ (B/4) * (cos(sigma0)*(-1 + 2 * pow( cos(sigma_m2), 2 ) - (B/6) * cos(sigma_m2) *(-3 + 4 * pow(sin(sigma0), 2 )) *(-3 + 4 * pow( cos (sigma_m2), 2 )))))
                        sigma = sigma0  + delta_sigma
                        while ( abs(sigma-sigma0)>0.00001) :
                                sigma0=sigma
                                sigma_m2 = (2 * sigma1 + sigma0)
                                delta_sigma = B * sin(sigma0) * ( cos(sigma_m2)+ (B/4) * (cos(sigma0) *(-1 + 2 * pow( cos(sigma_m2), 2 ) -(B/6) * cos(sigma_m2) *(-3 + 4 * pow(sin(sigma0), 2 )) *(-3 + 4 * pow( cos (sigma_m2), 2 )))))
                                sigma = (s / (b * A)) + delta_sigma
                        #Calcul de la latitude 2
                        TanBeta2=(sin(Beta1) * cos(sigma) + cos(Beta1) * sin(sigma) * cos(alpha1))/(sqrt( pow(Sinalphae, 2)+pow(sin(Beta1) * sin(sigma) - cos(Beta1) * cos(sigma) * cos(alpha1), 2)))
                        latitude2 =atan(TanBeta2/(1-f))
                        #Calcul de la différence de longitude sur la sphère auxiliaire deltau
                        deltau=(sin(sigma)*sin(alpha1))/(cos(Beta1)*cos(sigma)-sin(Beta1)*sin(sigma)*cos(alpha1))
                        #Calcul de la constante C de Vincenty
                        C = (f/16) * cosalphae_2 * (4 + f * (4 - 3 * cosalphae_2 ))
                        #Calcul de la différence de longitude dellambda
                        deltalambda = deltau - (1-C) * f * Sinalphae *(sigma + C * sin(sigma) * (cos(sigma_m2) +C * cos(sigma) * (-1 + 2 * pow(cos(sigma_m2),2) )))
                        longitude2 = longitude1 + deltalambda
                        #Calcul de l'azimut alpha2 et de l'azimut inverse
                        alpha2 = atan2 ( Sinalphae, (cos(Beta1) * cos(sigma) * cos(alpha1)-sin(Beta1) * sin(sigma)))
                        if ( alpha2 < pi ) :
                                alpha21 = alpha2 + pi
                        if ( alpha2 > pi ) :
                                alpha21 = alpha2 - pi
                        #transformation des angles en degre
                        latitude2= latitude2*180/pi
                        longitude2= longitude2*180/pi
                        alpha21 = alpha21*180/pi
                        #if longitude2>180:
                        #longitude2 = longitude2-360
                        #if longitude2<0:
                        # longitude2 = longitude2+360'''
                        return round(longitude2), round(latitude2), round(alpha21)

def say_hello(request):
        return render(request ,'direct.html')

def say_hello2(request):
        return render(request ,'inverse.html')

def acceuil(request):
        return render(request ,'acceuil.html')



def add(request):
        latitude1d=int(request.GET['d'])
        latitude1mm=int(request.GET['min'])
        latitude1ss=int(request.GET['ss'])
        latitude1=(int(latitude1d) +int(latitude1mm)/60+int(latitude1ss)/3600)*pi/180
        longitude1d=int(request.GET['dlo'])
        longitude1mm=int(request.GET['minlo'])
        longitude1ss=int(request.GET['sslo'])
        longitude1=(int(longitude1d) +int(longitude1mm)/60+int(longitude1ss)/3600)*pi/180
        s=float(request.GET['s'])
        a=float(request.GET['a'])
        b=float(request.GET['b'])
        alpha1=float(request.GET['az'])
        res1,res2,res3=direct(latitude1, longitude1, alpha1, s, a, b)
        return render( request   ,'result.html', {'result1': res1 ,'result2': res2 ,'result3': res3 } )




def inverse(request):
        def inversee(phi1,phi2,lam1,lam2,a,b):
                if a<0 or b<0:
                        texte="a et b doivent être positifs"
                        return texte
                        #hhhh
                else:
                        f = (a-b)/a		
                        ep = sqrt((a*2-b2)/b*2)
                        #Calcul de la latitude réduite beta1 
                        TanBeta1 = (1-f) * tan(latitude1)
                        Beta1 = atan( TanBeta1 )
                        #Calcul de la latitude réduite beta0
                        CosBeta0 = cos(Beta1)*sin(alpha1)
                        Beta0 = acos(CosBeta0)
                        #Calcul de la constante w2 
                        w2 = (sin(Beta0)*ep)*2
                        #Calcul de la distance angulaire sigma1 sur la sphère auxiliaire 
                        sigma1 = atan2( TanBeta1, cos(alpha1) )
                        #Calcul de l'azimut de la géodésique à l'équateur alphae
                        Sinalphae = CosBeta0
                        alphae = asin(Sinalphae)
                        cosalphae_sq = 1.0 - Sinalphae**2
                        #Calcul des constantes de Vincenty A' et B'
                        A = 1.0 + (w2 / 16384) * (4096 + w2 * (-768 + w2 * (320 - 175 * w2) ) )
                        B = (w2 / 1024) * (256 + w2 * (-128 + w2 * (74 - 47 * w2) ) )
                        #Calcul de la distance angulaire sigma sur le grand cercle de la sphère auxiliaire 
                        sigma = (s / (b * A))
                        sigma_m2 = (2 * sigma1 + sigma)
                        while ( abs(sigma-sigma0)>0.00001) :
                                sigma_m2 = (2 * sigma1 + sigma)
                                delta_sigma = B * sin(sigma) * ( cos(sigma_m2) 
                                        + (B/4) * (cos(sigma) * 
                                        (-1 + 2 * pow( cos(sigma_m2), 2 ) - 
                                        (B/6) * cos(sigma_m2) * \
                                        (-3 + 4 * pow(sin(sigma), 2 )) *  \
                                        (-3 + 4 * pow( cos (sigma_m2), 2 ))))) 
                                sigma = (s / (b * A)) + delta_sigma
                        #Calcul de la latitude 2
                        TanBeta2=(sin(Beta1) * cos(sigma) + cos(Beta1) * sin(sigma) * cos(alpha1)) \
                        /(sqrt( pow(Sinalphae, 2)+pow(sin(Beta1) * sin(sigma) - cos(Beta1) * cos(sigma) * cos(alpha1), 2)))
                        latitude2 =atan(TanBeta2/(1-f))
                        #Calcul de la différence de longitude sur la sphère auxiliaire deltau
                        deltau=(sin(sigma)*sin(alpha1))/(cos(Beta1)*cos(sigma)-sin(Beta1)*sin(sigma)*cos(alpha1))
                        #Calcul de la constante C de Vincenty 
                        C = (f/16) * cosalphae_sq * (4 + f * (4 - 3 * cosalphae_sq ))
                        #Calcul de la différence de longitude dellambda
                        deltalambda = deltau - (1-C) * f * Sinalphae *\
                                (sigma + C * sin(sigma) * (cos(sigma_m) + \
                                C * cos(sigma) * (-1 + 2 * pow(cos(sigma_m),2) )))
                
                        longitude2 = longitude1 + deltalambda
                        #Calcul de l'azimut alpha2 et de l'azimut inverse 
                        alpha2 = atan2 ( Sinalphae, (cos(Beta1) * cos(sigma) * cos(alpha1)\
                                -sin(Beta1) * sin(sigma)))
                
                        alpha21 = alpha2 + 2*pi / 2.0
                        if ( alpha21 < 0.0 ) :
                                alpha21 = alpha21 + 2*pi
                        if ( alpha21 > 2*pi ) :
                                alpha21 = alpha21 - 2*pi
                #transformation des angles en degre 
                        latitude2= latitude2*180/pi
                        longitude2= longitude2*180/pi
                        alpha21 = alpha21*180/pi
                        if longitude2>180:
                                longitude2 = longitude2-360
                        if longitude2<180:
                                longitude2 = longitude2+360
                        return round(longitude2), round(latitude2), round(alpha21)


                """f = (a-b)/a
                if phi1==90 or phi1==-90:
                        beta1=phi1
                else:
                        beta1=atan((1-f)*tan(phi1))
                if phi2==90 or phi2==-90:
                        beta2=phi2
                else:
                        beta2=atan((1-f)*tan(phi2))
                '''Calcul de la différence de longitude sur l'ellipsoïde'''
                delta_lam=lam2-lam1
                delta_u0= delta_lam
                #Calcul par itération
                #calcul des grandeurs suivantes
                sin2_sigma=pow(cos(beta2)*sin(delta_lam), 2) + pow((cos(beta1)*sin(beta2) - sin(beta1)*cos(beta2)*cos(delta_lam)),2)
                Sin_sigma = sqrt(sin2_sigma )
                Cos_sigma = sin(beta1) * sin(beta2) + cos(beta1)*cos(beta2)*cos(delta_lam)
                sigma = atan2(Sin_sigma, Cos_sigma)
                Sin_alphaE = cos(beta1)* cos(beta2)*sin(delta_lam)/ sin(sigma)
                alphaE = asin(Sin_alphaE)
                Cos2sigma_m = cos(sigma) - (2 * sin(beta1) * sin(beta2) / pow(cos(alphaE), 2) )
                #C Constante de Vincenty
                C = (f/16) * pow(cos(alphaE), 2) * (4 + f * (4 - 3 * pow(cos(alphaE), 2)))
                delta_u = delta_lam + (1-C)*f*sin(alphaE)*(sigma + C * sin(sigma)*(Cos2sigma_m + C*cos(sigma)*(-1 + 2 * pow(Cos2sigma_m, 2))))
                while delta_u - delta_u0>1.0e-5:
                # calcul des grandeurs suivantes
                        sin2_sigma = pow(cos(beta2) * sin(delta_u), 2) + pow((cos(beta1) * sin(beta2) - sin(beta1) * cos(beta2) * cos(delta_u)), 2)
                        Sin_sigma = sqrt(sin2_sigma)
                        Cos_sigma = sin(beta1) * sin(beta2) + cos(beta1) * cos(beta2) * cos(delta_u)
                        sigma = atan2(Sin_sigma, Cos_sigma)
                        Sin_alphaE = cos(beta1) * cos(beta2) * sin(delta_u) / sin(sigma)
                        alphaE = asin(Sin_alphaE)
                        Cos2sigma_m = cos(sigma) - (2 * sin(beta1) * sin(beta2) / pow(cos(alphaE), 2))
                        # C Constante de Vincenty
                        C = (f / 16) * pow(cos(alphaE), 2) * (4 + f * (4 - 3 * pow(cos(alphaE), 2)))
                        delta_u = delta_u0
                        delta_u0= delta_lam + (1 - C) * f * sin(alphaE) * (sigma + C * sin(sigma) * (Cos2sigma_m + C * cos(sigma) * (-1 + 2 * pow(Cos2sigma_m, 2))))
                '''Calcul de la latitude réduite de vertex de la géodésie'''
                beta0=acos(Sin_alphaE)
                #Calcul de la constante w_carre
                w_carre = pow(cos(alphaE),2) * (a*a-b*b) / (b*b)
                #Calcul des constantes de Vincenty A et B
                A = 1 + (w_carre/16384) * (4096 + w_carre * (-768 + w_carre * (320 - 175 * w_carre)))
                B = (w_carre/1024) * (256 + w_carre * (-128+ w_carre * (74 - 47 * w_carre)))
                #Calcul de delta_sigma
                delta_sigma = B * Sin_sigma * (Cos2sigma_m + (B/4)*(Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2) )-(B/6) * Cos2sigma_m * (-3 + 4 *sin2_sigma)*(-3 + 4 * pow(Cos2sigma_m,2 ) )))
                #Calcul de la distance géodésique s
                s = b * A * (sigma - delta_sigma)
                #Calcul de l'azimut directe alpha12
                alpha12 = atan2( (cos(beta2) * sin(delta_lam)),(cos(beta1) * sin(beta2) - sin(beta1) * cos(beta2) * cos(delta_lam)))
                #Calcul de l'azimut inverse alpha21
                alpha21 = atan2( (cos(beta1) * sin(delta_lam)),(-sin(beta1) * cos(beta2) + cos(beta1) * sin(beta2) * cos(delta_lam)))
                if ( alpha12 < 0.0 ) :
                        alpha12 =  alpha12 + 2*pi
                if ( alpha12 > 2*pi ) :
                        alpha12 = alpha12 - 2*pi
                        alpha21 = alpha21 + 2*pi / 2.0
                if ( alpha21 < 0.0 ) :
                        alpha21 = alpha21 + 2*pi
                if ( alpha21 > 2*pi ) :
                        alpha21 = alpha21 - 2*pi


                alpha12 = alpha12*180/pi
                alpha21 = alpha21*180/pi
                return(alpha12,alpha21,s)"""
        aa=float(request.GET['aa'])
        bb=float(request.GET['bb'])
        phi1d=float(request.GET['d'])
        phi1mm=float(request.GET['min'])
        phi1ss=float(request.GET['ss'])
        phi2d=float(request.GET['d2'])
        phi2mm=float(request.GET['min2'])
        phi2ss=float(request.GET['ss2'])
        lam1d=float(request.GET['dlo'])
        lam1mm=float(request.GET['minlo'])
        lam1ss=float(request.GET['sslo'])
        lam2d=float(request.GET['dlo2'])
        lam2mm=float(request.GET['minlo2'])
        lam2ss=float(request.GET['sslo2'])
        phi1= (phi1d + phi1mm/60+ phi1ss/3600)*pi/180
        phi2= (phi2d + phi2mm/60+ phi2ss/3600)*pi/180
        lam1= (lam1d + lam1mm/60+ lam1ss/3600)*pi/180
        lam2= (lam2d + lam2mm/60+ lam2ss/3600)*pi/180
        
        systeme=request.GET['sys']
        if systeme =="1":
                a=6378137
                b=a*(1-1/298.257223563)
        elif systeme =="2":
                a=6378135
                b=a*(1-1/298.26)
        elif systeme =="3":
                a=6378145
                b=a*(1-1/298.25)
        elif systeme =="4":
                a=6378165
                b=a*(1-1/298.3)
        elif systeme =="5":
                a=6378160
                b=a*(1-1/298.25)
        elif systeme =="6":
                a=6378245
                b=a*(1-1/298.3)
        elif systeme =="7":
                a=6378288
                b=a*(1-1/297)
        elif systeme =="8":
                a=6378270
                b=a*(1-1/297)
        elif systeme =="9":
                a=6378137
                b=a*(1-1/298.257222101)
        elif systeme =="10":
                a=6378140
                b=a*(1-1/298.257) 
        elif systeme =="11":
                a=6378160
                b=a*(1-1/298.247167427)  
        elif systeme =="12":
                a=6378150
                b=a*(1-1/298.3)
        elif systeme =="13":
                a=6378166
                b=a*(1-1/298.3)
        elif systeme =="14":
                a=6377276.345
                b=a*(1-1/300.8017)
        elif systeme =="15":
                a=6378249.145
                b=a*(1-1/293.465)
        elif systeme =="16":
                a=6378206.4
                b=a*(1-1/294.9786982)
        elif systeme =="17":
                a=6377397.155
                b=a*(1-1/299.1528128)
        elif systeme =="18":
                a=6377563.396
                b=a*(1-1/299.3249646)


        if  aa==0 or bb==0 or not aa or not bb:
                x=a
                y=b
        else:
                x=aa
                y=bb
        res=inversee(phi1,phi2,lam1,lam2,x,y)
        return render( request ,'resultinv.html', {'result': res } )


def visualisation(longitude1 , latitude1,alpha1,s,a,b):
        #visualisation de l'ellipsoide
        n=500
        t=s/(n-1)
        az_dir=[]
        lat=[]
        long=[]
        az_dir.append(alpha1)
        lat.append(latitude1)
        long.append(longitude1) 
        for i in range(0,n-1) :
                alpha1=az_dir[i]
                longitude1=long[i] 
                latitude1=lat[i]
                aa,bb,cc=direct(latitude1,longitude1,alpha1,s,a,b)
                lat.append(aa)
                long.append(bb)
                az_dir.append(cc)
        x=[]
        y=[]
        z=[]
        for i in range(0,n): 
                x.append((a*cos(long[i])*cos(lat[i]))/sqrt(1-(a**2-b**2)/(b**2)*(sin(lat[i]))**2))
                y.append((a*sin(long[i])*cos(lat[i]))/sqrt(1-(a**2-b**2)/(b**2)*(sin(lat[i]))**2))
                z.append((a*(1-(a**2-b**2)/(b**2))*sin(lat[i]))/sqrt(1-(a**2-b**2)/(b**2)*(sin(lat[i]))**2))
        return x,y,z 
def ellipsoide(longitude1, latitude1,alpha1,s,a,b):
        # some math: generate points on the surface of ellipsoid
        
        phi = np.linspace(0, 2*pi,800)
        theta = np.linspace(-pi/2, pi/2,800)
        phi, theta=np.meshgrid(phi, theta)

        x = cos(theta) * sin(phi) * a
        y = cos(theta) * cos(phi) * a

        z = sin(theta)*b

        # to use with Jupyter notebook

        init_notebook_mode()

        fig = go.Figure(data=[go.Mesh3d(
        x=x.flatten(), y=y.flatten(), z=z.flatten(), color='green', opacity=1.0, alphahull=0)])
        fig.write_html('first_figure.html', auto_open=True)

        fig.add_scatter3d(
        x1,y1,z1=visualisation(longitude1 , latitude1,alpha1,s,a,b),
        mode="lines"
        )
        plot_div = plot(fig, output_type='div', include_plotlyjs=False)
        return plot_div
        

















def essaiA(request):
        return render(request ,'essaiA.html')
        
def essaiB(request):
        
       
        
        f1=int(request.GET["f1"])
        
         
        res=f1
        return render( request ,'essaiB.html', {'result': res } )