from pygost import gost34112012256
import random
from Cryptodome.Util.number import inverse, GCD
from Cryptodome.Random.random import randint
import os
from elliptic_curve import *
from sys import argv
import asnGenerator
from sympy import gcd
from sympy import *
from Cryptodome.Util.number import isPrime, getPrime
import time
#from sympy import *


# ----------------------------------------------------------------
# Parameters
# Характеристика поля
#p = 307
# Порядок группы
#q = 167

#  Коэффициенты уравнения кривой
#a = 1
#b = 11
#p_x = 306
#p_y = 304

#q_x = 146
#q_y = 65



#d = 64921815105644748118349040703417588965511539353917797572335798592964210429984 

def check_cmth(a,b):
    tmp = 4*pow(a,3)+27*pow(b,2)
    if tmp == 0:
        return False
    return True

def sum(x_P, y_P, x_Q, y_Q, a, b, p):
    if (x_P == 0 or x_P is None):
        return x_Q, y_Q
    if (x_Q == 0 or x_Q is None):
        return x_P, y_P

    obr_y_P = (p - y_P)

    if (x_P == x_Q and y_P == y_Q):
      
        #if (GCD(2*y_P,p)!=1):
        #    return None, None
        #print("2*y_P = ", 2*y_P)
        m = ((3*x_P**2 + a) * invert(2*y_P,p))
        m = int(m)
        x_R = (m ** 2 - x_P - x_Q) % p
        #x_R = int(x_R)
        y_R = (m*(x_P - x_R)-y_P)%p

        return x_R, y_R

    if (x_P == x_Q and obr_y_P == y_Q):
        #print("P_BESK!!")
        return None, None
    
    
    #print("x_P-x_Q = ",x_P-x_Q, ", p = ",p)
    #print("GCD (x_P-x_Q,p) = ", GCD(x_P-x_Q,p))
    #print("invert(x_P-x_Q,p)) = ", invert(x_P-x_Q,p))
    

    #print("x_P-x_Q = ", x_P-x_Q % p)
    m = ((y_Q-y_P) * invert(x_Q-x_P,p)) #%p#((y_P-y_Q) * invert(x_P-x_Q,p)) %p
    m = int(m)
    x_R = (m ** 2 - x_P-x_Q) %p
    #x_R = int(x_R)
    y_R = (m*(x_P - x_R) - y_P)%p

    return x_R, y_R

    
    return None

def _mul_points (x_P, y_P, number, a, b, p):

    if number == 0:
        return None, None
    k = number
    k_binary =  [int(i) for i in bin(k)[2:]]
    t = len(k_binary)
    x_Q = None
    y_Q = None
    i = 0
    while i<t:
        #if (x_Q is not None or y_Q is not None):
            #print("Q <- 2Q  = (",x_Q,", ",y_Q,") + (",x_Q,", ",y_Q,")")
        x_Q, y_Q = sum(x_Q, y_Q, x_Q, y_Q, a, b, p)
        
        if k_binary[i]!=0:
            #print("Q <- Q + P = (",x_Q,", ",y_Q,") + (",x_P,", ",y_P,")")
            x_Q, y_Q = sum(x_P,y_P,x_Q,y_Q,a,b,p )

        i = i+1

    return x_Q, y_Q



def H(P_x, L):
    return P_x % L

def _pollard(L,  p, q, a, b, p_x, p_y, q_x, q_y):

    print("q = ",q)
    print("a = ",a)
    print("p = ",p)
    print("p_x = ",p_x)
    print("p_y = ",p_y)
    a_j = [None] * L
    b_j = [None] * L
    R_j_x = [None] * L
    R_j_y = [None] * L
    alg = 1
    while True:
        #print(" ---> STEP #3: Genetating a_j, b_j, R_j...")
        #j=1
        #a_j = [0 for i in range(L+1)]
        #b_j = [0 for i in range(L+1)]
        #R_j_x = [0 for i in range(L+1)]
        #R_j_y = [0 for i in range(L+1)]
        for j in range(L):#while j<=L:
            a_j[j] = randint(2,q-1)
            b_j[j] = randint(2,q-1)
            #print("               a[{j}] = ", a_j[j])
            #print("               b[{j}] = ", b_j[j])
            ###print("a_j*P = (",p_x,", ",p_y,") * ",a_j[j])
            aP_x, aP_y = _mul_points(p_x,p_y,a_j[j], a, b, p)
            #print("GOt aP = (",aP_x, " ", aP_y, ")")
            ###print("b_j*Q = (", q_x,", ",q_y,") * ",b_j[j])
            bQ_x, bQ_y = _mul_points(q_x,q_y,b_j[j], a, b, p)
            #print("GOt bQ = (",bQ_x, " ", bQ_y, ")")
            ###print("R_j = a_j*P + b_j*Q = (",aP_x,", ",aP_y,") + (",bQ_x,", ",bQ_y,")")
            R_j_x[j], R_j_y[j] = sum(aP_x,aP_y,bQ_x,bQ_y, a, b, p)
            #print("               R[{j}] = <", R_j_x[j],", ",R_j_y[j],">")
            j = j + 1

        #print(" ---> STEP #4: Genetating alpha, beta, T...")
        alpha = randint(2,q-1)
        betta = randint(2,q-1)
        #print("               alpha' = ", alpha)
        #print("               betta' = ", betta)

        alphaP_x, alphaP_y = _mul_points(p_x,p_y, alpha, a, b, p)
        bettaQ_x, bettaQ_y = _mul_points(q_x,q_y, betta, a, b, p)
        
        T_x, T_y = sum(alphaP_x, alphaP_y,bettaQ_x, bettaQ_y, a, b, p)
        #print("               T' = <", T_x,", ",T_y,">") 

        _T_x = T_x
        _T_y = T_y 
        _alpha = alpha
        _betta = betta

        #print(" ---> STEP #5: While T' != T''...")
        cycle = 1
        while True:

            j = H(T_x, L)
            T_x, T_y = sum(T_x, T_y, R_j_x[j], R_j_y[j], a, b, p)
            ###print("#",cycle,":   f^i(P) = (",T_x,", ",T_y,")")
            alpha = (alpha + a_j[j]) % q
            betta = (betta + b_j[j]) % q
            #print("               T' = <", T_x,", ",T_y,">") 
            #print("               alpha' = ", alpha)
            #print("               betta' = ", betta)

            j = H(_T_x, L)
            _T_x, _T_y = sum(_T_x, _T_y, R_j_x[j], R_j_y[j], a, b, p)
            _alpha = (_alpha + a_j[j]) % q
            _betta = (_betta + b_j[j]) % q
            j = H(_T_x, L)
            _T_x, _T_y = sum(_T_x, _T_y, R_j_x[j], R_j_y[j], a, b, p)
            _alpha = (_alpha + a_j[j]) % q
            _betta = (_betta + b_j[j]) % q
            ###print("#",cycle,":   f^2i(P) = (",_T_x,", ",_T_y,")")
            #print("               T'' = <", _T_x,", ",_T_y,">") 
            #print("               alpha'' = ", _alpha)
            #print("               betta'' = ", _betta)

            if T_x == _T_x and T_y == _T_y:
                #print("               At step ",cycle," T' == T'' !!!")
                break
            cycle = cycle + 1
        
        #print(" ---> STEP #6: Checking if alpha'==alpha'' and betta'==betta''...")
        if alpha == _alpha or betta == _betta:
            #print("               alpha'==alpha'' or betta'==betta'' !!!")
            #print("               Chose another L !!!")
            d = -1
            #return False
            alg = alg + 1
        else:
            if GCD(_betta - betta, q)==1:
                d = (alpha - _alpha)*invert(_betta - betta, q) % q
                print("               alpha'!=alpha'' and betta'!=betta'' !!!")
                print("               FOUND D !!!")
                print("               d = ", d)

                print("\n ---> Checking...")
                print("Calculating Q=dP...")
                tmp_x, tmp_y = _mul_points(p_x, p_y, d, a, b, p)
                print("Initial Q = <", q_x, ", ",q_y,">")
                print("Calculd Q = <", tmp_x, ", ",tmp_y,">")

                if (q_x==tmp_x and q_y==tmp_y):
                    print("d found on #", alg)
                    return True
            alg = alg + 1

            
        

    return False



def generate_prime(q):

    while True:
        k = randint(1, q - 1)

        if gcd(k, q) == 1:
            return k


def main():

    # Parameters
    # Характеристика поля
    p = 307
    # Порядок группы
    q = 167

    #  Коэффициенты уравнения кривой
    a = 1
    b = 11
    p_x = 306
    p_y = 304

    q_x = 146
    q_y = 65

    
    p = 2774052499
    a = 2552774921
    b = 1967537144
    d = 629204683
    q = 924688404
    p_x = 2083077157
    p_y = 1745053455
    q_x = 2508372582
    q_y = 1108100667

    
    #p = 79
    # Порядок группы
    #q = 83

    #  Коэффициенты уравнения кривой
    #a = 3
    #b = 9
    #p_x = 31
    #p_y = 6

    #q_x = 24
    #q_y = 53

    
    print(" ---> STEP #1: Choosing L...")
    L = 32

    print("               L = ", L)

    print(" ---> STEP #2: Choosing H(<P>)...")
    print("               H(<P>) = P.x % L")
    
    if _pollard(L, int(p), int(q), int(a), int(b), int(p_x), int(p_y), int(q_x), int(q_y)):
        print("all good :)")
    else:
        print("oh oh :(")

if __name__ == "__main__":
   main()




