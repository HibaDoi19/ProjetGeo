from math import *
def direct(latitude1, longitude1, alpha1, s, a, b):

    #VERIFICATION DES PARAMèTRE DE L'ELLIPSOIDE
    if a<0 or b<0 or b>a:
        texte="a et b doivent être positifs et b doit être inférieur à a"
        return texte
    #vérification de s
    elif s<0 :
        texte="erreur, distance doit être positive"
        return texte
    #vérification de l'azimut
    elif alpha1<0 or alpha1>360 :
        texte="Erreur : azimut doit être entre 0° et 360°"
        return texte
    elif abs(latitude1)>pi/2 :
        texte="Erreur: valeur de la latitude hors rang"
        return texte
    elif abs(longitude1)>pi :
        texte= "erreur : longitude hors rang ! "
        return texte
    else:
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
        if ( alpha2 <= pi ) :
                alpha21 = alpha2 + pi
        if ( alpha2 > pi ) :
                alpha21 = alpha2 - pi
    #transformation des angles en degre
        latitude2= latitude2*180/pi
        #if latitude2>90:
             #latitude2 =90-(latitude2-90)
        #if latitude2<180:
             #latitude2 =-90-(latitude2+90)  
        longitude2= longitude2*180/pi
        alpha21 = alpha21*180/pi
        #if alpha21<0:
            #alpha21=alpha21+2*pi
        #if alpha21>360:
           # alpha21=alpha21-2*pi
       # if longitude2>180:
          #   longitude2 = longitude2-360
        #if longitude2<0:
            # longitude2 = longitude2+360
        if alpha1<pi and longitude2<(longitude1*180/pi) and longitude1<0:
                longitude2 = longitude2+180
        if alpha1<pi and longitude1>0 and longitude2<(longitude1*180/pi) : 
                longitude2=longitude2-180  
        return longitude2,round(latitude2),round(alpha21)