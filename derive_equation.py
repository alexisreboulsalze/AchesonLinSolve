# -*-A coding: utf-8 -*-
import sympy as sym
import numpy as np

sym.init_printing(use_unicode=True, wrap_line=True)

#Variables definition
Bp = sym.Symbol('Bp',positive=True,real=True) #magnetic field
gr=sym.Symbol('gr',real=True)
gz=sym.Symbol('gz',real=True)
rho=sym.Symbol('rho',positive=True,real=True)
cs2=sym.Symbol('cs2',positive=True,real=True) #careful compared to acheson a^2 = \gamma cs^2
Omega=sym.Symbol('Omega',positive=True,real=True)
Pressure=sym.Symbol('P',positive=True,real=True)
gamma=sym.Symbol('gamma',real=True)
tilde_kappa=sym.Symbol('tilde_kappa',positive=True,real=True)
tilde_nu=sym.Symbol('tilde_nu',positive=True,real=True)
tilde_eta=sym.Symbol('tilde_eta',positive=True,real=True)
m=sym.Symbol('m',integer=True)
n=sym.Symbol('n',real=True)
l=sym.Symbol('l',real=True)
r=sym.Symbol('r',real=True)
I=sym.I
i=sym.I #in case still mistake
#Derivative
dOmedh=sym.Symbol('dOmegadh',real=True)
dOmer2dh=sym.Symbol('dOmegar2dh',real=True)
#Composed variables
s2=sym.Symbol('s2',positive=True,real=True)#l**2+n**2
V=sym.Symbol('V',positive=True,real=True)#Bp/sym.sqrt(4*sym.pi*rho)
F=sym.Symbol('dF',real=True)#sym.log(Bp/(rho*r))
Q=sym.Symbol('dQ',real=True)#sym.log(Bp*r)
E=sym.Symbol('dE',real=True)#sym.log(Pressure*rho**(-gamma))
p=sym.Symbol('p',real=True)#l**2+n**2
#ouput
w=sym.Symbol("w",complex=True)

G= sym.Symbol('G',real=True)#gr - (l/n)*gz

#For polynomial in w: multiply by expr by (w*gamma+I*tilde_kappa)*(w+I*eta*s2)*(w*gamma+I*tilde_kappa)
expr1_poly=V**2*(2*Omega*(w*gamma+I*tilde_kappa)*m/r + (w+I*tilde_nu)*(2*(w*gamma+I*tilde_kappa)/r-G*(w+I*tilde_kappa)/(cs2/gamma)))*(m*(w*gamma+I*tilde_kappa)*dOmedh+F*(w+I*tilde_eta)*(w*gamma+I*tilde_kappa) -w*(w+I*tilde_eta)*E)
expr2_poly=(s2/n**2*((w+I*tilde_nu)*(w+I*tilde_eta)*(w*gamma+I*tilde_kappa)  - (m**2*V**2/r**2)*(w*gamma+I*tilde_kappa))- G*(w+I*tilde_eta)*E)*((w+I*tilde_nu)*(w+I*tilde_eta)*(w*gamma+I*tilde_kappa)-(m*V/r)**2*(w*gamma+I*tilde_kappa) +V**2/(cs2/gamma)*(w+I*tilde_kappa)*w*(w+I*tilde_nu))
expr3_poly=-((dOmer2dh*(w+I*tilde_eta)+m*V**2*Q)*(2*Omega/r*((w+I*tilde_eta)*(w*gamma+I*tilde_kappa)+V**2/(cs2/gamma)*(w+I*tilde_kappa)*w)+m*V**2/r**2*(2*(w*gamma+I*tilde_kappa)/r - G/(cs2/gamma)*(w+I*tilde_kappa)))*(w*gamma+I*tilde_kappa))

expr_poly=expr1_poly+expr2_poly+expr3_poly
test=sym.poly_from_expr(expr_poly,w)

degree=sym.polys.polytools.degree(expr_poly, gen=w)
coeffs=[]
terms=[]
expr_tmp=expr_poly
f=open("Acheson_equation.py","w")
f.write("import numpy as        np\n")
f.write(" \n")
f.write("def Acheson_coeffs(Omega,q,N,rho,g,r,theta,Bp,eta,nu,kappa,cs2,gamma,p,s,l,n,s2,m):\n")
f.write("    dOmegadh=q*Omega/r\n")
f.write("    dOmegar2dh=(q+2)*Omega*r\n")
f.write("    c=3e10\n    V_non_rel=Bp/np.sqrt(4*np.pi*rho)\n    V=V_non_rel*c/(c**2+V_non_rel**2)**(1/2)\n")
f.write("    dE=gamma*N**2/g*(np.sin(theta)-l/n*np.cos(theta))\n")
#f.write("    H= cs2/g\n")
f.write("    m=1\n")
f.write("    I=1j\n")
f.write("    z=r*np.cos(theta)\n")
f.write("    dQ=(p+1)/r -l/n*s/z #no z dependence\n")
f.write("    dF=(p-1)/r -l/n*s/z #no z dependence\n")
f.write("    gr=g*np.sin(theta)\n")
f.write("    gz=g*np.cos(theta)\n")
f.write("    tilde_eta=eta*s2\n")
f.write("    tilde_kappa=kappa*s2\n")
f.write("    tilde_nu=nu*s2\n")
f.write("    G= gr - gz*(l/n)\n")
#f.write("    G= g*np.sin(theta) - g*np.cos(theta)*(l/n)\n")
f.write("    coeffs=[]\n")
for iii in range(degree+1):
    print(sym.polys.polytools.degree(expr_tmp, gen=w))
    LT_tmp=sym.polys.polytools.LT(expr_tmp,w)
    print(LT_tmp)
    terms.append(LT_tmp)
    coeffs.append(sym.simplify(sym.polys.polytools.LC(expr_tmp,w)))
    expr_tmp=expr_tmp-LT_tmp
    f.write("    coeffs.append(")
    f.write(str(coeffs[iii]))
    f.write(")\n")

f.write("    return(coeffs)\n")
f.close()

print(sym.latex(expr_poly))


#SPRUIT A2 equation *(w*gamma+i*kappa*s2)*(w+i*eta*s2)*(w*gamma+i*kappa*s2)
#Hypothesis Pole (gr=0, gz=g), nu=0, dOmegadh=0, dOmer2dh=2*Omega*r
#H= sym.Symbol('H',positive=True,real=True)#Pressure scaleheight
H=cs2/gz
expr1_S99 = (V**2/r*(2*Omega*m*(w*gamma+I*tilde_kappa) + w*(2*(w*gamma+I*tilde_kappa)+l*r/n/H*(w+I*tilde_kappa)))*(F*(w*gamma+I*tilde_kappa)-w*E))*(w+I*tilde_eta)
expr2_S99 = (s2/n**2*(w*(w+I*tilde_eta)*(w*gamma +I*tilde_kappa) - (w*gamma +I*tilde_kappa)*(m*V/r)**2) + l/n*(gz*(w+I*tilde_eta)*E))*(w*(w+I*tilde_eta)*(w*gamma +I*tilde_kappa) - (w*gamma +I*tilde_kappa)*(m*V/r)**2 + V**2/cs2*(w+I*tilde_kappa)*w**2)
expr3_S99 = - (2*Omega*(w+I*tilde_eta) + m*V**2/r**2*(p+1))*(2*Omega*((w+I*tilde_eta)*(w*gamma +I*tilde_kappa)+V**2/cs2*(w+I*tilde_kappa)*w)+ m*V**2/r**2*(2*(w*gamma +I*tilde_kappa)/r + l/n*r/H*(w+I*tilde_kappa)))*(w*gamma +I*tilde_kappa)

expr_S99_poly = expr1_S99 +expr2_S99+expr3_S99
degree_S99=sym.polys.polytools.degree(expr_S99_poly, gen=w)
coeffs_S99=[]
terms_S99=[]
expr_tmp=expr_S99_poly
f=open("Spruit_equation.py","w")
f.write("import numpy as	np\n")
f.write(" \n")
f.write("def Spruit_coeffs(Omega,q,N,rho,g,r,theta,Bp,eta,nu,kappa,cs2,gamma,p,s,l,n,s2,m):\n")
f.write("    dOmegadh=q*Omega/r\n")
f.write("    dOmegar2dh=(q+2)*Omega*r\n")
f.write("    c=3e10\n    V_non_rel=Bp/np.sqrt(4*np.pi*rho)\n    V=V_non_rel*c/(c**2+V_non_rel**2)**(1/2)\n")
f.write("    dE=gamma*N**2/g*(np.sin(theta)-l/n*np.cos(theta))\n") #polytropic gaz hypothesis for now
f.write("    H= cs2/g\n")
f.write("    m=1\n")
f.write("    I=1j\n")
f.write("    z=r*np.cos(theta)\n")
f.write("    dQ=(p+1)/r -l/n*s/z #no z dependence\n")
f.write("    dF=(p-1)/r -l/n*s/z #no z dependence\n")
f.write("    tilde_eta=eta*s2\n")
f.write("    tilde_kappa=kappa*s2\n")
f.write("    tilde_nu=nu*s2\n")
f.write("    gz=g\n")
f.write("    coeffs=[]\n")
for iii in range(degree_S99+1):
    print(sym.polys.polytools.degree(expr_tmp, gen=w))
    LT_tmp_S99=sym.polys.polytools.LT(expr_tmp,w)
    #print(LT_tmp_S99)
    terms_S99.append(LT_tmp_S99)
    coeffs_S99.append(sym.simplify(sym.polys.polytools.LC(expr_tmp,w)))
    expr_tmp=expr_tmp-LT_tmp_S99
    f.write("    coeffs.append(")
    f.write(str(coeffs_S99[iii]))
    f.write(")\n")

f.write("    return(coeffs)\n")
f.close()

print(sym.latex(expr_S99_poly))


#Simplified Acheson
#gr=0.0
#Omega=0.0
tilde_nu=0.0
tilde_eta=0.0
tilde_kappa=0.0
E=0.0
#dOmedh=0.0
#G= gr - (l/n)*gz
#Q=(p+1)/r
#dOmer2dh=2*Omega
#tilde_eta=0.0

#For polynomial in w: multiply by expr by (w*gamma+I*tilde_kappa)*(w+I*eta*s2)*(w*gamma+I*tilde_kappa)
expr1_poly=V**2*(2*Omega*(w*gamma+I*tilde_kappa)*m/r + (w+I*tilde_nu)*(2*(w*gamma+I*tilde_kappa)/r-G*(w+I*tilde_kappa)/(cs2/gamma)))*(m*(w*gamma+I*tilde_kappa)*dOmedh+F*(w+I*tilde_eta)*(w*gamma+I*tilde_kappa) -w*(w+I*tilde_eta)*E)
expr2_poly=(s2/n**2*((w+I*tilde_nu)*(w+I*tilde_eta)*(w*gamma+I*tilde_kappa)  - (m**2*V**2/r**2)*(w*gamma+I*tilde_kappa))- G*(w+I*tilde_eta)*E)*((w+I*tilde_nu)*(w+I*tilde_eta)*(w*gamma+I*tilde_kappa)-(m*V/r)**2*(w*gamma+I*tilde_kappa) +V**2/(cs2/gamma)*(w+I*tilde_kappa)*w*(w+I*tilde_nu))
expr3_poly=-((dOmer2dh*(w+I*tilde_eta)+m*V**2*Q)*(2*Omega/r*((w+I*tilde_eta)*(w*gamma+I*tilde_kappa)+V**2/(cs2/gamma)*(w+I*tilde_kappa)*w)+m*V**2/r**2*(2*(w*gamma+I*tilde_kappa)/r - G/(cs2/gamma)*(w+I*tilde_kappa)))*(w*gamma+I*tilde_kappa))

expr_poly=expr1_poly+expr2_poly+expr3_poly
test=sym.poly_from_expr(expr_poly,w)

degree=sym.polys.polytools.degree(expr_poly, gen=w)
coeffs=[]
terms=[]
expr_tmp=expr_poly
f=open("Acheson_equation_diffusionless.py","w")
f.write("import numpy as        np\n")
f.write(" \n")
f.write("def Acheson2_coeffs(Omega,q,N,rho,g,r,theta,Bp,eta,nu,kappa,cs2,gamma,p,s,l,n,s2,m):\n")
f.write("    dOmegadh=q*Omega/r\n")
f.write("    dOmegar2dh=(q+2)*Omega*r\n")
f.write("    c=3e10\n    V_non_rel=Bp/np.sqrt(4*np.pi*rho)\n    V=V_non_rel*c/(c**2+V_non_rel**2)**(1/2)\n")
f.write("    dE=gamma*N**2/g*(np.sin(theta)-l/n*np.cos(theta))\n") #polytropic gaz hypothesis for now
#f.write("    H= cs2/g\n")
f.write("    m=1\n")
f.write("    I=1j\n")
f.write("    z=r*np.cos(theta)\n")
f.write("    dQ=(p+1)/r -l/n*s/z #no z dependence\n")
f.write("    dF=(p-1)/r -l/n*s/z #no z dependence\n")
f.write("    gr=g*np.sin(theta)\n")
f.write("    gz=g*np.cos(theta)\n")
f.write("    tilde_eta=eta*s2\n")
f.write("    tilde_kappa=kappa*s2\n")
f.write("    tilde_nu=nu*s2\n")
f.write("    G= gr - gz*(l/n)\n")
#f.write("    G= g*np.sin(theta) - g*np.cos(theta)*(l/n)\n")
f.write("    coeffs=[]\n")
for iii in range(degree+1):
    tmp_degree=sym.polys.polytools.degree(expr_tmp, gen=w)
    print("degree=",tmp_degree)
    if tmp_degree != -np.infty:
        if tmp_degree== degree-iii:
            LT_tmp=sym.polys.polytools.LT(expr_tmp,w)
            print(LT_tmp)
            terms.append(LT_tmp)
            f.write("    coeffs.append(")
            f.write(str(sym.simplify(sym.polys.polytools.LC(expr_tmp,w))))
            f.write(")\n")
            expr_tmp=expr_tmp-LT_tmp
        else:
            f.write("    coeffs.append(0.0)\n")
f.write("    return(coeffs)\n")
f.close()

print(sym.latex(expr_poly))
