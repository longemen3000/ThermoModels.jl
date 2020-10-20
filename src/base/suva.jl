








#==

sussesive Van der Waals aproximations

idea: helmolhtz equations are costly to evaluate
what if, instead of using a root finder, we aproximate via interpolation via van der waals?


VdW: P(1+a/v2)(v-b) = RT 

P = RT/(v-b) - a/v2 





procedure: obtaining v₁

relax factor = 0.9 = λ

if aproximation is from liquid:

    using: a/v2 ≈ 0
    b = RT/P₀ - v₀ #estimated covolume, do not test volumes lesser than this
    v₁ = λv₀ + (1-λ)b
    P₁ = P(v₁)
    if P₁ < Pspec:
        P₀,v₀ = P₁,v₁, REPEAT
    else:
        return  P₀,v₀,P₁,v₁

if aproximation is from gas:
    at high volumes:
        P ≈ RT/v - a/v2
        Pv2 = RTv - a 
        a = RTv - Pv2 = v(RT-Pv)
        v2 - RT/P v + a/P = 0
        v = 0.5*(RT/P + sqrt((RT/P)^2 - 4a/P))
    v₁ = v₀/λ
    P₁ = P(v₁)
    if P₁ > Pspec:
        P₀,v₀ = P₁,v₁, REPEAT
    else:
        return  P₀,v₀,P₁,v₁
    
with  P₀,v₀,P₁,v₁, we can start aproximating P

the basis is the following:

1.  P₀ = RT/(v₀-b) - a/v₀^2  # * v₀^2
    P₁ = RT/(v₁-b) - a/v₁^2  #  *v₁^2 

2.  P₀v₀^2 = RTv₀^2 /(v₀-b) - a  
    P₁v₁^2 = RTv₁^2 /(v₁-b) - a   

joining:
3.  P₁v₁^2 - P₀v₀^2   = RTv₁^2 /(v₁-b) - RTv₀^2 /(v₀-b) # * (v₁-b)*(v₀-b)

4. (P₁v₁^2 - P₀v₀^2)*(v₁-b)*(v₀-b) = RTv₁^2(v₀-b) - RTv₀^2(v₁-b)

5. (P₁v₁^2 - P₀v₀^2)*(v₁-b)*(v₀-b) = RTv₁^2(v₀-b) - RTv₀^2(v₁-b)

# M = (P₁v₁^2 - P₀v₀^2)

6. M*(v₁-b)*(v₀-b) = -(RTv₁^2 - RTv₀^2)*b + (RTv₀v₁^2 - RTv₁v₀^2)

# N = RT(v₁^2 - v₀^2), L = RT(v₀v₁^2 - v₁v₀^2)

7. = (v₁-b)*(v₀-b) + (N/M)b - L/M == 0

#expanding:
b^2 (- N/M-v₁-v₀) b + v₁v₀ - L/M

using cuadratic solving: b is obtained:

#

then a is obtaining as follows:

(RT/(v₁-b)-P₁)=  a #using ₁ as is the most recent aproximation

with b and , now we have a local VdW aproximation, we can obtain the new P as the following:

A = a*Pspec/(RT)^2
B = b*Pspec/RT

#generalized cubic compressibility factor polynomial
Z3 + ((c1+c2-1)B -1)Z2 + (c1c2B^2 - (c1+c2)*(B2+B)+A)Z - B(c1c2(B2+B)+A)

#in VdW: c1 = c2 = 0
Z3 + ((-1)B -1)Z2 + (+A)Z - B(+A)
Z3 + -(B +1)Z2 + AZ - AB
s
using cubic solving, a trio of z (zsols) is obtained

if liquid: 
    z = minimum(reals(zsols))
if gas
    z = maximum(reals(zsols))

vnew = zRT/Pspec
Pnew = P(vnew)

check convergence, if not:
P₀,v₀,P₁,v₁ = P₁,v₁,Pnew,vnew,  REPEAT

==#
"""
uses van der waals aproximations to calculate the function p(v) = p0

"""
function suva_fzero(p,Pspec,v0,t;
    liquid=true,
    relax=0.9*one(float(Pspec)),
    maxevals=100,
    atol = zero(float(Pspec)),
    rtol = 8*eps(one(float(Pspec))))
    pspec = Pspec

    _Pspec,_v0,_t,_atol,_rtol,_relax = promote(pspec,v0,t,atol,rtol,relax)
    return _suva_fzero(p,_Pspec,_v0,t,liquid,_relax,maxevals,_atol,_rtol)
end

function _suva_fzero(P,pspec::TT,v0::TT,T::TT,liquid=true,relax = TT(0.9),maxevals=100, atol=zero(TT), rtol=8*eps(TT)) where TT
    Pspec = pspec
    R = RGAS #constant defined outside, 8.314...
    v₀ = v0
    P₀= P(v0)

    λ = relax
    _1 = one(TT)
    _0 = zero(TT)
    P₁ = _0
    v₁ = _0
    nan = _0/_0
    RT = R*T
    RT_P = RT/pspec
    invRT = _1/RT
    preA = Pspec*invRT*invRT
    preB = Pspec*invRT
    count = 0

    if liquid
        while count < 110
            # at liquid volumes, a/v2 ≈ 0
            # p ≈ RT/(v-b)
            b = RT/P₀ - v₀ #estimated covolume, do not test volumes lesser than this
            #@show b
            if b < 0 #patological v0, changing to z aprox
                P₀= P(v0)
                b = RT/P₀ - v₀
                #@show v₀,P₀
                #@show b
            end
      
            v₁ = 1.5*b
            P₁ = P(v₁)
            #@show v₁,P₁

            if (P₁ < Pspec < P₀) | (P₀ < Pspec < P₁) #bounded P, excelent
                break
            elseif (P₀ < P₁< Pspec) #step in the right direction, P1 increases, but is lower than Pspec
               # @show P₁
                P₀,v₀ = P₁,v₁
                count +=1 
            elseif (P₁ < P₀ < Pspec) #step in the wrong direction P1 decreased, danger zone
                v₁ = v₀*0.9
                P₁ = P(v₁)
               # @show P₁
            elseif min(Pspec,P₀,P₁) == Pspec #we have increasing pressures, good
                P₀,v₀ = P₁,v₁
                b = RT/P₀ - v₀
                v₁ = RT_P - b
                P₁ = P(v₁)
               # @show P₁
                break
            elseif count == 100 
                return error("cannot find initial value for gas")
            else
                P₀,v₀ = P₁,v₁
               # @show P₁
                count +=1 
            end
        end
    else
        while count < 110
            #at high volumes, P ≈ RT/v - a/v2
            a = v₀*(RT-P₀*v₀)
           # @show a
            v₁ = 0.5*(RT_P + sqrt(RT_P^2 - 4*a/Pspec))
            P₁ = P(v₁)
            if (P₁ < Pspec < P₀) | (P₀ < Pspec < P₁) #bounded P, excelent
                break
            elseif (Pspec < P₁ < P₀) #step in the right direction, P₁ < P₀, but Pspec < P₁
                P₀,v₀ = P₁,v₁
                count +=1 
            elseif (Pspec < P₀ < P₁) #step in the wrong direction, reduce volume directly
                v₁ = v₀*1.1
                P₁ = P(v₁)
            elseif max(Pspec,P₀,P₁) == Pspec #very good, and extra iteration to refine, as (P₀,v₀) could be very bad
                P₀,v₀ = P₁,v₁
                a = v₀*(RT-P₀*v₀)
               # @show a
                v₁ = 0.5*(RT_P + sqrt(RT_P^2 - 4*a/Pspec))
                P₁ = P(v₁)
                break
            elseif count == 100 
                return error("cannot find initial value for gas")
            else
                P₀,v₀ = P₁,v₁
                count +=1 
            end
        end
    end

    evalscount = 0 
    uatol = atol / oneunit(atol) * oneunit(real(v₁))
    adjustunit = oneunit(real(P₀))/oneunit(real(v₀))
    vnew = v₁
    Pnew = P₁
    
    while evalscount < maxevals

        v₁2 = v₁*v₁
        v₀2 = v₀*v₀
        v₀3 = v₀2*v₀
        v₁3 = v₁2*v₁
   
        M = (P₁*v₁2 - P₀*v₀2)
        _a = M
        _b = M*(v₀+v₁)
        _c =  RT*(v₀ - v₁) - M*v₁*v₀
        sqrtΔ = sqrt(_b*_b - 4*_a*_c)
        adiv = 0.5/_a
        b1 = (-_b + sqrtΔ)*adiv#covolume
        b2 = (-_b - sqrtΔ)*adiv
        b = b1
        a = (RT/(v₁-b) - P₁)*v₁2 #temperature coefficient
        
        A = a*preA
        B = b*preB
        #@show -A*B
        #Z3 + ((-1)B -1)Z2 + (+A)Z - B(+A)
        #-AB + AZ +(-B-1)Z2 + Z3
        poly = (-A*B, A, -B-_1, _1)
        zvals = cardan(poly)
        #realvals = cardan2(poly)
        realvals = cardan_reals(zvals)
        #@show zvals
        #@show realvals
        #realvals = filter(>(_0),_realvals)
        z = maximum(realvals)
        
        vnew = z*RT_P #stops negative volumes
        Pred = RT/(vnew-b) - a/(vnew*vnew)
        #@show Pred
        #@show (a,b)
        if vnew < 0 
            vnew= 0.5*(v₀+v₁)
        end
        Pnew = P(vnew)
        ΔP = Pnew-Pspec
        ΔP₀ = P₀-Pspec
        ΔP₁ = P₁-Pspec
        #@show ΔP,ΔP₁,ΔP₀
        

        if (abs(ΔP) > abs(ΔP₀)> abs(ΔP₁))
            correction = "corrected spec,P₀ is better"
         #   @show correction
            vnew = v₀ - (v₀-v₁)*ΔP₀/(P₀ - P₁)
            #@show vnew
            if vnew < 0 #correction again for negative volumes
                vnew= 0.5*(v₀+v₁)
            end
            Pnew = P(vnew)
            ΔP = Pnew-Pspec
            #@show ΔP,ΔP₁,ΔP₀
        elseif (abs(ΔP) > abs(ΔP₁))
            correction ="corrected spec P₁ is better"
          #  @show correction
            vnew = v₁ - (v₀-v₁)*ΔP₁/(P₀ - P₁)
            if vnew < 0 #correction again for negative volumes
                vnew= 0.5*(v₀+v₁)
            end
            Pnew = P(vnew)
            ΔP = Pnew-Pspec
            #@show ΔP,ΔP₁,ΔP₀
        else  
            correction ="no correction"
           # @show correction
        end
 
        vcor = (vnew-b)
        pcor = (Pnew+a/(vnew*vnew))
        #@show  vcor,pcor
        #@show  vcor*pcor
        #@show  RT
        if ΔP == eps(Pnew)
            return vnew
        end
        #@show ΔP,ΔP₁,ΔP₀
      #  @show Pnew,P₁,P₀
        #progress = abs(ΔP)/abs(ΔP₁)
        #@show progress
        abs(ΔP) <= adjustunit * max(uatol, abs(vnew) * rtol) && return vnew
        iszero(ΔP) && return vnew
        isnan(ΔP) && return nan
        if (ΔP == ΔP₁)
            sdp = sign(ΔP)
            v_sprev = sign(P(prevfloat(vnew))-Pspec)
            v_snext = sign(P(prevfloat(vnew))-Pspec)
            p_sprev = sign(prevfloat(ΔP))
            p_snext = sign(nextfloat(ΔP))
            sdp * v_sprev <= 0 && return vnew
            sdp * v_snext <= 0 && return vnew
            sdp * p_sprev <= 0 && return vnew
            sdp * p_snext <= 0 && return vnew
            return nan
        end

            P₀,v₀,P₁,v₁ = P₁,v₁,Pnew,vnew
        evalscount += 1
        #@show evalscount
    end
    return vnew
end


# P = RT/v-b - a/v2
#  a= v2(RT/v-b - P)
#(p+a/v2)(v-b) > 0


# P (1- aV/RT) = b

#==
v₀2 = v₀^2
v₁2 = v₁^2
v₀3 = v₀2*v₀
v₁3 = v₁2*v₁
M = (P₁*v₁2 - P₀*v₀2)
p₀*v₀2()
_b = (v₁+v₀)
b = -b +1- (1/2M)*sqrt(_b^2 - 4M(k(v₀ - v₁) + P₀*v₀3*v₁ - P₁*v₀*v₁3)) + P₀ v₀3 + P₀ v₀^2 v₁ - p₁ v₀ v₁^2 - p₁ v₁3)/(2 (P₀ v₀^2 - p₁ v₁^2))
==#

# Z = PV/RT

