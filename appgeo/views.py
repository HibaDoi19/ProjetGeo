
from re import template
from tkinter import E, W
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
from . import ploot
from . import direct 
from . import innverse
from . import geocentrique2geodesique

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
                aa=float(request.GET['aa'])
                bb=float(request.GET['bb'])
                alpha1=float(request.GET['az'])
                alpha1=alpha1*pi/180

                drc_lat=int(request.GET["drc_lat"])
                drc_long=int(request.GET["drc_long"])

                X_centrique= request.GET["cx"]
                Y_centrique= request.GET["cy"]
                Z_centrique= request.GET["cz"]


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
                latitude1=latitude1*drc_lat
                longitude1=longitude1*drc_long

                long_des,lat_des,h_des=geocentrique2geodesique.geocentrique2geodesique(X_centrique,Y_centrique,Z_centrique,x,y)
                long_des=long_des*pi/180
                lat_des=lat_des*pi/180

                systeme=request.GET["multi_note"]
                if systeme =="1":
                        latitude1=latitude1
                        longitude1=longitude1

                                
                elif systeme =="2":
                        latitude1=lat_des
                        longitude1=long_des

                                


                #if request.methode =="post" and "Calculer" in request.POST :
                res1,res2,res3=direct.direct(latitude1, longitude1, alpha1, s, x, y)
                long_d=(int(res1))
                long_m=(int((res1-long_d)*60))
                long_s=(int(((res1-long_d)*60-long_m)*60)) 
                lat_d=(int(res2))
                lat_m=(int((res2-lat_d)*60))
                lat_s=(int(((res2-lat_d)*60-lat_m)*60))
                if long_d >=0 and long_d >=0 and long_d >=0 :
                        jj="E"
                else :
                        jj="W"
                        long_d=abs(long_d)
                        long_m=abs(long_m)
                        long_s=abs(long_s)
                if lat_d >=0 and lat_d >=0 and lat_d >=0 :
                        hh="N"
                else :
                        hh="S"
                        lat_d=abs(lat_d)
                        lat_m=abs(lat_m)
                        lat_s=abs(lat_s)


                long_d=str(long_d)+"°"
                long_m=str(long_m)+"'"
                long_s=str(long_s)+"''"+jj
                lat_d=str(lat_d)+"°"
                lat_m=str(lat_m)+"'"
                lat_s=str(lat_s)+"''"+hh

                res3=str(res3)+"°"






                return render( request   ,'result.html', {'longituded': long_d ,'longitudem': long_m ,'longitudes': long_s ,'latituded': lat_d ,'latitudem': lat_m,'latitudes': lat_s,'angle_azhimutale': res3 ,'plot':ploot.ellipsoide(longitude1,latitude1, alpha1, s, x, y)} )
                #return render( request   ,'result.html', {'longitude': res1 ,'latitude': res2 ,'angle azhimutale': res3 ,'plot':ploot.ellipsoide(longitude1,latitude1, alpha1, s, x, y)} )
        
def inverse(request):         
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
        
        drc_lat_1=int(request.GET["drc_lat_1"])
        drc_long_1=int(request.GET["drc_long_1"])
        drc_lat_2=int(request.GET["drc_lat_2"])
        drc_long_2=int(request.GET["drc_long_2"])

        X1= request.GET["cx1"]
        Y1= request.GET["cy1"]
        Z1= request.GET["cz1"]

        X2= request.GET["cx2"]
        Y2= request.GET["cy2"]
        Z2= request.GET["cz2"]

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
        phi1=phi1*drc_lat_1
        phi2=phi2*drc_lat_2
        lam1=lam1*drc_long_1
        lam2=lam2*drc_long_2

        lat1,long1,h1=geocentrique2geodesique.geocentrique2geodesique(X1,Y1,Z1,x,y)
        lat2,long2,h2=geocentrique2geodesique.geocentrique2geodesique(X2,Y2,Z2,x,y)

        systeme=request.GET["multi_note"]
        if systeme =="1":
                latitude1=lam1
                longitude1=phi1
                latitude2=lam2
                longitude2=phi2

                                
        elif systeme =="2":
                latitude1=lat1
                longitude1=long1
                latitude2=lat2
                longitude2=long2

        xx,yy,zz=innverse.inversee(latitude1,latitude2,longitude1,longitude2,x,y)
        xx=float(xx)
        yy=float(yy)
        zz=float(zz)

        #return render( request ,'resultinv.html', {'ALPHA12': round(xx),'ALPHA21': round(yy),'s': round(zz) })
        return render( request ,'resultinv.html', {'ALPHA12': round(xx),'ALPHA21': round(yy),'s': round(zz) ,'plot':ploot.ellipsoide(lam1,phi1, xx*pi/180, zz, x, y)})

        