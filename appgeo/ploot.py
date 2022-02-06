from re import template
from typing import Text
from _plotly_utils.utils import template_doc
from django.shortcuts import render
from django.http import HttpResponse
from math import *
from plotly.offline import iplot, init_notebook_mode, plot 
import plotly.graph_objects as go
from numpy import FLOATING_POINT_SUPPORT, sin, cos, pi
import  numpy as np
import plotly.express as px 
import datetime
import glob
import logging
import os


def visualisation(longitude1 , latitude1,alpha1,s,a,b,h):

        #visualisation de l'ellipsoide
        h=int(h)
        n=500
        t=s/(n-1)
        az_dir=[]
        lat=[]
        long=[]
        az_dir.append(alpha1)
        lat.append(latitude1)
        long.append(longitude1) 
        for i in range(1,n) :
                aa,bb,cc=directe(latitude1,longitude1,alpha1,t,a,b)
                lat.append(bb)
                long.append(aa)
                az_dir.append(cc)
                alpha1=az_dir[i]
                longitude1=long[i] 
                latitude1=lat[i]
        x=[]
        y=[]
        z=[]
        for i in range(2,n): 
                e=sqrt((a**2-b**2)/(a**2))
                N=a/sqrt(1-e**2*sin(lat[i])**2)
                x.append((N+h)*cos(long[i])*cos(lat[i]))
                y.append((N+h)*sin(long[i])*cos(lat[i]))
                z.append((N*(1-e**2)+h)*sin(lat[i]))
        return x,y,z 

def directe(latitude1, longitude1, alpha1, s, a, b):
        #VERIFICATION DES PARAMèTRE DE L'ELLIPSOIDE

                f = (a-b)/a
                ep = sqrt((a**2-b**2)/b**2)
                #Calcul de la latitude réduite beta1
                if (abs(latitude1) != pi/2) :
                        TanBeta1 = (1-f) * tan(latitude1)
                        Beta1 = atan(TanBeta1 )
                else:
                        Beta1=latitude1
                #Calcul de la latitude réduite beta0
                CosBeta0 = cos(Beta1)*sin(alpha1)
                Beta0 = acos(CosBeta0)
                #Calcul de la constante w2
                w2 = (sin(Beta0)*ep)**2
                #Calcul de la distance angulaire sigma1 sur la sphère auxiliaire
                if latitude1 == pi/2 or latitude1 == -pi/2 :
                        sigma1=Beta1
                elif alpha1==pi/2 or alpha1==3*pi/2 :
                        if latitude1==0 :
                                sigma1=0
                        else :
                                sigma1=pi/2
                else:
                        sigma1 = atan2( TanBeta1, cos(alpha1))
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
                deltau=atan((sin(sigma)*sin(alpha1))/(cos(Beta1)*cos(sigma)-sin(Beta1)*sin(sigma)*cos(alpha1)))
                #Calcul de la constante C de Vincenty
                C = (f/16) * cosalphae_2 * (4 + f * (4 - 3 * cosalphae_2 ))
                #Calcul de la différence de longitude dellambda
                deltalambda = deltau - (1-C) * f * Sinalphae *(sigma*pi/180 + C * sin(sigma) * (cos(sigma_m2) +C * cos(sigma) * (-1 + 2 * pow(cos(sigma_m2),2) )))
                longitude2 = (longitude1 + deltalambda)
                #Calcul de l'azimut alpha2 et de l'azimut inverse
                alpha2 = atan2 ( Sinalphae, (cos(Beta1) * cos(sigma) * cos(alpha1)-sin(Beta1) * sin(sigma)))
                alpha21=0
                if ( alpha2 < pi ) :
                        alpha21 = alpha2 + pi
                if ( alpha2 > pi ) :
                        alpha21 = alpha2 - pi
                if alpha1<pi and longitude2<(longitude1*180/pi) and longitude1<0:
                        longitude2 = longitude2+180
                if alpha1<pi and longitude1>0 and longitude2<(longitude1*180/pi) : 
                        longitude2=longitude2-180  
                if alpha1>pi and longitude2>(longitude1*180/pi) and longitude1<0: 
                        longitude2 = longitude2+180
                if alpha1>pi and longitude1>0 and longitude2>(longitude1*180/pi): 
                        longitude2=longitude2+180
                if longitude2>180:
                        longitude2=longitude2-360
                if longitude2<-180:
                        longitude2=longitude2+360
                                        

                return longitude2,latitude2,alpha2

def ellipsoide(longitude1 , latitude1,alpha1,s,a,b,h):                      
                phi = np.linspace(0, 2*pi)
                theta = np.linspace(-pi/2, pi/2)
                phi, theta=np.meshgrid(phi, theta)

                X = cos(theta) * sin(phi) * a
                Y = cos(theta) * cos(phi) * a
                Z = sin(theta)*b
                
                layout = go.Layout(width = 800, height =800,title_text='Géodesique')

                fig = go.Figure(data=[go.Surface(x = X, y = Y, z=Z, colorscale = 'Blues')], layout=layout)

                fig.update_traces(contours_z=dict(show=True, usecolormap=True,highlightcolor="limegreen", project_z=True))

                x,y,z=visualisation(longitude1,latitude1,alpha1,s,a,b,h)

                fig.add_scatter3d(x=x,y=y,z=z,mode='lines',marker ={'color':'black'})

                fig.add_scatter3d(x=[x[0]],y=[y[0]],z=[z[0]],mode='markers',marker ={'color':'yellow'})
                fig.add_scatter3d(x=[x[-1]],y=[y[-1]],z=[z[-1]],mode='markers',marker ={'color':'green'})
                ############################################
                for i in [-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80]:
                        phi = np.linspace(0, 360,100)
                        theta = np.linspace(i, i,100)
    
                        phi = phi*pi/180
                        theta = theta*pi/180
                        #phi, theta=np.meshgrid(phi, theta)
                        x = cos(theta) * sin(phi) * a
                        y = cos(theta) * cos(phi) * a
                        z = sin(theta)*b
                        t="latitude"+str(i)
                        fig.add_scatter3d(x =x, y = y, z=z,mode='lines', marker={'color':'red'},text=t)
                        
                        phi = np.linspace(0, 360,100)
                theta = np.linspace(0, 0,100)
                phi = phi*pi/180
                theta = theta*pi/180
                #phi, theta=np.meshgrid(phi, theta)
                x = cos(theta) * sin(phi) * a
                y = cos(theta) * cos(phi) * a
                z = sin(theta)*b
                t="équateur"
                fig.add_scatter3d(x =x, y = y, z=z,mode='lines', marker={'color':'#77004D'},text=t)
    
                for i in [-170,-160,-150,-140,-130,-120,-110,-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180]:
                        phi = np.linspace(i, i,100)
                        theta = np.linspace(-90,90,100)
                        
                        phi = phi*pi/180
                        theta = theta*pi/180
                        #phi, theta=np.meshgrid(phi, theta)
                        x = cos(theta) * sin(phi) * a
                        y = cos(theta) * cos(phi) * a
                        z = sin(theta)*b
                        t="longitude"+str(i)
                        fig.add_scatter3d(x =x, y = y, z=z,mode='lines', marker={'color':'red'},text=t)
                phi = np.linspace(0, 0,100)
                theta = np.linspace(-90,90,100)
                
                phi = phi*pi/180
                theta = theta*pi/180
                #phi, theta=np.meshgrid(phi, theta)
                x = cos(theta) * sin(phi) * a
                y = cos(theta) * cos(phi) * a
                z = sin(theta)*b
                t="Greenwich"
                fig.add_scatter3d(x =x, y = y, z=z,mode='lines', marker={'color':'#28B463'},text=t)
                        
                ###########################################

                plot_div = fig.to_html(full_html=False)

               
                
                return plot_div

