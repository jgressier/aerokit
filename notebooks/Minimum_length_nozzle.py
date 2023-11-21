import math as m
import matplotlib.pyplot as plt

# Toute constante est exprimée en unités S.I

gamma = 1.27

section_col = 1/10000*m.pi

largeur = 0.01

rayon_col = m.sqrt(section_col / m.pi)


def pmeyer(M): #calcule la fonction de Prandtl-Meyer

    c = (gamma-1)/(gamma+1)

    v = m.sqrt(1/c) * m.atan( m.sqrt( c*(M**2.0 -1) ) )

    v = v - m.atan( m.sqrt(M**2.0 - 1.0) )

    return v


def pmeyer_prime(M) : #dérivée de la fonction de Prandtl-Meyer

    if (M==1):

        return 0.0

    else:
        c = (gamma-1)/(gamma+1)

        v_prime = ( ( M / (1 + c*(M**2 - 1)) ) - 1/M ) / m.sqrt(M**2 - 1)

        return v_prime

def mach_via_pmeyer(v): #réciproque de Prandtl-Meyer via méthode de Newton

    if(v==0):

        return 1.0

    # On définit la fonction dont on va chercher le zero
    # ainsi que sa dérivée

    def f(m):

        return ( v - pmeyer(m) )

    def f_prime(m):

        return ( -1* pmeyer_prime(m) )

    M = 2.0

    while(abs(f(M)) > 0.00001):

        M = M - (f(M)/f_prime(M))

    return M


def mach_angle(M): # calcule l'angle de mach mu

    return m.asin(1/M)



# On défini une classe représentative d'un point du maillage

class point:

    # coordonnées spatiales

    x = 0.0
    y = 0.0

    # Angles informant sur l'écoulement (isentropique)

    theta = 0.0
    nu_pmeyer = 0.0

    # Constantes liées à la C+ et la C- passant par le point

    def k_plus(self):

        return self.theta - self.nu_pmeyer

    def k_moins(self):

        return self.theta + self.nu_pmeyer

    # Permet de copier le point, usage pratique

    def copie(self):

        copy = point()
        copy.x = self.x
        copy.y = self.y
        copy.theta = self.theta
        copy.nu_pmeyer = self.nu_pmeyer

        return copy


def new_point(caratab, a, b): # Crée le point numéro b de la C- numéro a et le renvoie

    # Parents du nouveau point

    father_Cmoins = caratab[a][b-1]
    father_Cplus = caratab[a-1][b]

    # Création du point à l'aide de l'eq de Prandtl-Meyer

    new = point()

    new.theta =  ( father_Cmoins.k_moins() + father_Cplus.k_plus() ) / 2

    new.nu_pmeyer = ( father_Cmoins.k_moins() - father_Cplus.k_plus() ) / 2

    # Angles de busemann associés aux zones d'émission des caractéristiques

    Beta_plus = ( father_Cplus.theta + mach_angle(mach_via_pmeyer(father_Cplus.nu_pmeyer)) + new.theta + mach_angle(mach_via_pmeyer(new.nu_pmeyer)) ) / 2
    Beta_moins = ( father_Cmoins.theta - mach_angle(mach_via_pmeyer(father_Cmoins.nu_pmeyer)) + new.theta - mach_angle(mach_via_pmeyer(new.nu_pmeyer)) ) / 2
    print("Beta_plus : ", Beta_plus * 180 / m.pi, "°")
    print("Beta_moins : ", Beta_moins * 180 / m.pi, "°")

    # On utilise ces angles pour déterminer la position du nouveau point

    new.y = (father_Cplus.y - father_Cmoins.y * m.tan(Beta_plus) / m.tan(Beta_moins) + (father_Cmoins.x - father_Cplus.x) * m.tan(Beta_plus)) / (1 - m.tan(Beta_plus) / m.tan(Beta_moins))

    new.x = father_Cmoins.x + (new.y - father_Cmoins.y) / m.tan(Beta_moins)

    return new


# La fonction suivante va générer un profil de divergent
# via la méthode de Busemann. Le résultat sera un tableau
# avec les coordonnées des points sur la tuyère

# Le delta_theta est l'incrément angulaire du divergent,
# qui sera la variable de précision de la modélisation.
# Le delta_l sera la longueur de chaque segment modélisant
# le divergent. Plus cette variable sera grande, plus la détente
# se fera progressivement, mais plus le divergent sera long et donc
# la tuyère lourde. Il s'agira du paramètre d'optimisation.

def Minimum_length_nozzle(mach_final, area_throat, n):

    theta0 = pmeyer(mach_final) / 2

    print("Theta0 = ", theta0*180/m.pi, "°")

    dtheta = theta0 / n / 100 # Angle de la première caractéristique, proche de la verticale

    delta_theta = (theta0 - dtheta) / n # Incrément angulaire

    print("Dtheta = ", dtheta*180/m.pi, "° & Delta theta = ", delta_theta*180/m.pi, "°")

    # Chaque point à une intersection de charactéristiques
    # sera indicé par la caractéristique descendante C-
    # sur laquelle il se trouve (a, indicé à 1) ainsi que son numéro sur cette C-
    # en partant de la paroi (b, indicé à 0). On trace a+1 points par C-

    points_ecoulement = []

    profil_divergent = [(0.0, m.sqrt(area_throat/m.pi))] # Liste contenant le profil à générer

    point_actuel = point()

    point_actuel.x = 0.0
    point_actuel.y = m.sqrt(area_throat/m.pi)


    # On crée la première caractéristique, proche de la verticale (car proche de la ligne sonique)

    point_actuel.x = 0.0
    point_actuel.y = m.sqrt(area_throat/m.pi)
    point_actuel.theta = dtheta
    point_actuel.nu_pmeyer = dtheta

    points_ecoulement.append([point_actuel.copie()])

    point_actuel.theta = dtheta
    point_actuel.nu_pmeyer = dtheta
    Beta = points_ecoulement[0][0].theta - mach_angle(mach_via_pmeyer(points_ecoulement[0][0].nu_pmeyer))
    point_actuel.x = points_ecoulement[0][0].x - points_ecoulement[0][0].y / m.tan(Beta)
    point_actuel.y = 0

    points_ecoulement[0].append(point_actuel.copie())

    plt.plot(point_actuel.x, point_actuel.y, 'ro')

    print("")
    print("Beta : ", Beta * 180 / m.pi, "°")
    print("abcisse axe : ", point_actuel.x)
    print("")

    for a in range(1,n+1): # On parcours le faisceau de détente, et donc les charactéristiques


        point_actuel.x = 0.0
        point_actuel.y = m.sqrt(area_throat/m.pi)
        point_actuel.theta = points_ecoulement[a-1][0].theta + delta_theta
        point_actuel.nu_pmeyer = point_actuel.theta

        points_ecoulement.append( [point_actuel.copie()] )


        for b in range(1,a+1): # On parcours les points, et donc les zones, sauf ceux sur le profil et sur l'axe

            point_actuel = new_point(points_ecoulement, a, b)
            points_ecoulement[a].append( point_actuel.copie() )

            print("")
            print("a = ", a, "b = ", b)
            print("x = ", points_ecoulement[a][b].x, "y = ", points_ecoulement[a][b].y)
            print("sigma = ", (point_actuel.theta + mach_angle(mach_via_pmeyer(point_actuel.nu_pmeyer))) *180/m.pi , "°")
            print("")

        # On ajoute le point sur l'axe

        point_actuel.theta = 0
        point_actuel.nu_pmeyer = points_ecoulement[a][a].k_moins()
        Beta = points_ecoulement[a][a].theta - mach_angle(mach_via_pmeyer(points_ecoulement[a][a].nu_pmeyer))
        point_actuel.x = points_ecoulement[a][a].x - points_ecoulement[a][a].y / m.tan(Beta)
        point_actuel.y = 0

        points_ecoulement[a].append(point_actuel.copie())

        print("")
        print("Beta : ", Beta * 180 / m.pi, "°")
        print("abcisse axe : ", point_actuel.x)
        print("")



    # Il y a n+1 ondes de détente réfléchies à annuler

    print("")
    print("")


    # On crée le point annulant la caractéristique initiale

    pere_mur = points_ecoulement[-1][0]
    pere_ecoulement = points_ecoulement[n][1]
    point_actuel = pere_ecoulement.copie() # Le point sur la paroi garde les propriétés du point duquel il est issu

    sigma = pere_ecoulement.theta + mach_angle(mach_via_pmeyer(pere_ecoulement.nu_pmeyer))

    theta_moy = (point_actuel.theta + pere_mur.theta) / 2

    point_actuel.y = (pere_ecoulement.y - pere_mur.y * m.tan(sigma) / m.tan(theta_moy) + (pere_mur.x - pere_ecoulement.x) * m.tan(sigma)) / (1 - m.tan(sigma) / m.tan(theta_moy))
    point_actuel.x = pere_ecoulement.x + (point_actuel.y - pere_ecoulement.y) / m.tan(sigma)

    points_ecoulement.append([point_actuel.copie()])

    profil_divergent.append((point_actuel.x, point_actuel.y))

    for k in range(2,n+2):

        # On va créer le point à la paroi annulant la k-ième détente

        pere_mur = points_ecoulement[-1][0]
        pere_ecoulement = points_ecoulement[n][k]
        point_actuel = pere_ecoulement.copie() # Le point sur la paroi garde les propriétés du point duquel il est issu

        sigma = pere_ecoulement.theta + mach_angle(mach_via_pmeyer(pere_ecoulement.nu_pmeyer))

        theta_moy = (point_actuel.theta + pere_mur.theta) / 2

        point_actuel.y = ( pere_ecoulement.y - pere_mur.y * m.tan(sigma) / m.tan(theta_moy) + (pere_mur.x - pere_ecoulement.x) * m.tan(sigma) ) / (1 - m.tan(sigma) / m.tan(theta_moy))
        point_actuel.x = pere_mur.x + (point_actuel.y - pere_mur.y) / m.tan(theta_moy)

        points_ecoulement.append( [point_actuel.copie()] )

        profil_divergent.append( (point_actuel.x, point_actuel.y) )

        print("")
        print("pere ecoulement x = ", pere_ecoulement.x, "y = ", pere_ecoulement.y, "theta = ", pere_ecoulement.theta * 180/m.pi, "°")
        print("pere mur x = ", pere_mur.x, "y = ", pere_mur.y, "theta = ", pere_mur.theta * 180 /m.pi, "°")

        print("")
        print("theta mur = ", theta_moy*180/m.pi, "°", " sigma = ", sigma*180/m.pi, "°")
        print("")
        print("Nouveau point paroi x = ", profil_divergent[-1][0], "y = ", profil_divergent[-1][1])

    print("nombre d'increments :", n*2)

    print("")
    print("Angle final :", points_ecoulement[-1][0].theta*180/m.pi, "°")

    # On va tracer le profil

    list_x = []
    list_y = []
    list_ymoins = []
    for k in range(len(profil_divergent)):
        list_x.append(profil_divergent[k][0])
        list_y.append(profil_divergent[k][1])
        list_ymoins.append(-1*profil_divergent[k][1])

    plt.plot(list_x, list_y)
    plt.plot(list_x, list_ymoins)


    for i in range(n+1):

        list_cara_x = []
        list_cara_y = []
        list_cara_ymoins = []

        for j in range(i+2):

            list_cara_x.append( points_ecoulement[i][j].x )
            list_cara_y.append( points_ecoulement[i][j].y )
            list_cara_ymoins.append(-1*points_ecoulement[i][j].y)


        plt.plot(list_cara_x, list_cara_y)
        plt.plot(list_cara_x, list_cara_ymoins)

    c = 0

    for i in range(1,n+2):

        list_cara_x = []
        list_cara_y = []
        list_cara_ymoins = []


        for j in range(i-1, n+1):

            list_cara_x.append( points_ecoulement[j][i].x )
            list_cara_y.append( points_ecoulement[j][i].y )
            list_cara_ymoins.append(-1 *points_ecoulement[j][i].y)

            if points_ecoulement[j][i].x != 0 :

                c += 1
                print("Point ", c, "theta = ", points_ecoulement[j][i].theta*180/m.pi, "° pmeyer = ", points_ecoulement[j][i].nu_pmeyer*180/m.pi)

        c += 1
        print("Point ", c, "theta = ", points_ecoulement[j][i].theta * 180 / m.pi, "° pmeyer = ", points_ecoulement[j][i].nu_pmeyer * 180 / m.pi)

        list_cara_x.append( points_ecoulement[n + i][0].x )
        list_cara_y.append( points_ecoulement[n + i][0].y )
        list_cara_ymoins.append(-1 *points_ecoulement[n + i][0].y)
        plt.plot(list_cara_x, list_cara_y)
        plt.plot(list_cara_x, list_cara_ymoins)


    plt.axis([ 0.0, profil_divergent[-1][0], -1 * profil_divergent[-1][0], profil_divergent[-1][0] ])


    plt.show()


    return (profil_divergent, mach_via_pmeyer(point_actuel.nu_pmeyer))


def rapport_expansion(mach): # Calcule le rapport d'aire par rapport au col pour aobtenir le Mach souhaité

    return m.pow( (2 + (gamma-1) * (mach**2)) / (gamma + 1), (gamma + 1) / 2 / (gamma - 1)) / mach


mach_sortie = 6

result = Minimum_length_nozzle(mach_sortie, section_col, 20)
print("")

print("Longueur du divergent :", result[0][-1][0]*100, "cm")
print("")
print("Hauteur de sortie : ", result[0][-1][1]*100*2, "cm")
print("Hauteur de sortie théorique : ", rapport_expansion(mach_sortie) * rayon_col * 2 * 100, "cm")
print("")
print("Mach à la sortie :",result[1])
print("")
print("Rapport d'expansion théorique : ", rapport_expansion(mach_sortie))
print("Rapport d'epansion obtenu : ", result[0][-1][1] / rayon_col)

print("")
print("Rayon au col : ", rayon_col*100, "cm")
