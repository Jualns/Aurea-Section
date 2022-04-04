function intervalo(p, f)
    a, s, b = BigFloat(0.0), p, 2*p
    while f(b) < f(s)
        a, s, b = s, b, 2b
    end
    return a,b
end

function  aurea(a, b, ϵ, f; θ1 = (BigFloat(3.0) - BigFloat(sqrt(5)))/BigFloat(2.0), θ2 = (BigFloat(-1.0) + BigFloat(sqrt(5)))/BigFloat(2.0))
    a, s, b = BigFloat(0.0), ρ, 2ρ
    while f(b) < f(s)
        a, s, b = s, b, 2b
    end

    u, v = a + θ1*(b - a), a + θ2*(b - a)
    while (b - a) > ϵ
        if f(u) < f(v)
            b, v, u = v , u, a + θ1*(b - a)
        else
            a, u, v = u, v, a + θ2*(b - a)
        end
    end
    return (u + v)/2
end
ϵ = 1e-8 #algum ϵ
ρ = 0.1 #algum ρ
a,b = intervalo(ρ, x -> 4 + x^2 )#f(x)
print(aurea(a, b, ϵ, x -> 4 + x^2))
