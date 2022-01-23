from math import *       
def inversee(phi1,phi2,lam1,lam2,a,b):
                if a<0 or b<0:
                        texte="a et b doivent être positifs"
                        return texte
                else:
                        f = (a-b)/a
                #condition sur phi1
                if phi1 >(pi/2) or phi1<-(pi/2) :
                        print('la valeur de la latitude doit être inférieur ou égal à 90°')
                elif phi1==(pi/2) or phi1==-(pi/2):
                        beta1=phi1
                else:
                        beta1=atan((1-f)*tan(phi1))
                #condition sur phi2
                if phi2 > (pi/2) or phi2 < -(pi/2):
                        print('la valeur de la latitude doit être inférieur ou égal à 90°')
                elif phi2 == (pi/2) or phi2 == -(pi/2):
                        beta2=phi2
                else:
                        beta2=atan((1-f)*tan(phi2))
                        #condition sur lam1
                if lam1>pi or lam1<-pi :
                        print('la valeur de la longétude doit être inférieur ou égal à 180°')
                #condition sur lam2
                if lam2>pi or lam2<-pi :
                        print('la valeur de la longétude doit être inférieur ou égal à 180°')
                #Calcul de la différence de longitude sur l ellipsoïde
                delta_lam=lam2-lam1
                delta_u0= delta_lam
                #Calcul par itération
                #calcul des grandeurs suivantes
                sin2_sigma=pow(cos(beta2)*sin(delta_lam), 2) + pow((cos(beta1)*sin(beta2) - sin(beta1)*cos(beta2)*cos(delta_lam)),2)
                Sin_sigma = sqrt(sin2_sigma )
                Cos_sigma = sin(beta1) * sin(beta2) + cos(beta1)*cos(beta2)*cos(delta_lam)
                sigma = atan2(Sin_sigma, Cos_sigma)
                Sin_alphaE = cos(beta1)* cos(beta2)*sin(delta_lam)/ Sin_sigma
                alphaE = asin(Sin_alphaE)
                Cos2sigma_m = Cos_sigma - (2 * sin(beta1) * sin(beta2) / pow(cos(alphaE), 2) )
                #C Constante de Vincenty
                C = (f/16) * pow(cos(alphaE), 2) * (4 + f * (4 - 3 * pow(cos(alphaE), 2)))
                delta_u = delta_lam + (1-C)*f* sin(alphaE)*(sigma + C * Sin_sigma *(Cos2sigma_m + C* Cos_sigma *(-1 + 2 * pow(Cos2sigma_m, 2))))
                while abs(delta_u - delta_u0)>1.0e-5:
                        # calcul des grandeurs suivantes
                        sin2_sigma = pow(cos(beta2) * sin(delta_u), 2) + pow((cos(beta1) * sin(beta2) - sin(beta1) * cos(beta2) * cos(delta_u)), 2)
                        Sin_sigma = sqrt(sin2_sigma)
                        Cos_sigma = sin(beta1) * sin(beta2) + cos(beta1) * cos(beta2) * cos(delta_u)
                        sigma = atan2(Sin_sigma, Cos_sigma)
                        Sin_alphaE = cos(beta1) * cos(beta2) * sin(delta_u) / Sin_sigma
                        alphaE = asin(Sin_alphaE)
                        Cos2sigma_m = Cos_sigma - (2 * sin(beta1) * sin(beta2) / pow(cos(alphaE), 2))
                        # C Constante de Vincenty
                        C = (f / 16) * pow(cos(alphaE), 2) * (4 + f * (4 - 3 * pow(cos(alphaE), 2)))
                        delta_u0 = delta_u
                        delta_u= delta_lam + (1 - C) * f * sin(alphaE) * (sigma + C * Sin_sigma * (Cos2sigma_m + C * Cos_sigma * (-1 + 2 * pow(Cos2sigma_m, 2))))
                #Calcul de la latitude réduite de vertex de la géodésie
                beta0=acos(Sin_alphaE)
                #Calcul de la constante w_carre
                w_carre =((a*a-b*b) / (b*b))*pow(sin(beta0),2)
                #Calcul des constantes de Vincenty A' et B'
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
                return alpha12, alpha21, s
                #return round(alpha12), round(alpha21), round(s)