using LinearAlgebra, Printf, SparseArrays

function lu_tri(d, e, f)
    for k = 1:length(e)
        e[k] = e[k]/d[k]
        d[k+1] = (d[k+1]-e[k]*f[k])
    end
end

function resolve_tri(d, e, f, b)
    for k = 2:length(b)
        b[k] = b[k] - e[k-1]*b[k-1]
    end
    b[length(b)] = b[length(b)]/d[length(b)]
    for i = length(b)-1 :-1:1
        b[i] = (b[i] - f[i]*b[i+1])/d[i]
    end
end

####################
# Não mexer abaixo #
####################

function main()
  E_tri = 0.0
  E_sis = 0.0
  all = 0
  contador = 0
  for n = [5; 10; 50; 100; 500; 1000; 5000; 10000; 50000; 100000]
    @printf("Avaliando com n = %7d", n)
    Δt = time()
    failed = false
    for t = 1:100
      t % 5 == 0 && print(".")
      contador += 1
      d = 10 .+ rand(n)
      e = randn(n-1)
      f = randn(n-1)

      A = spdiagm(0 => d, -1 => e, 1 => f)
      b = A * ones(n)

      all += @allocated lu_tri(d, e, f)
      L = spdiagm(0 => ones(n), -1 => e)
      U = spdiagm(0 => d, 1 => f)
      all += @allocated resolve_tri(d, e, f, b)

      E_tri += norm(L * U - A)
      E_sis += norm(b .- 1)
      if time() - Δt > 10
        failed = true
        break
      end
    end
    if !failed
      @printf(" Pronto. Erros até agora: %8.2e  %8.2e\n", E_tri / contador, E_sis / contador)
    end
    Δt = time() - Δt
    if Δt > 10
      println("\033[31mExcesso de tempo ocorrido, verifique que seu código está otimizado\033[0m")
      break
    end
  end
  E_tri /= contador
  E_sis /= contador
  println("Erro no LU: $E_tri - Deveria ser perto de 1e-16")
  println("Erro na resolução do sistema: $E_sis - Deveria ser perto de 1e-16")
  println("Alocações: $all - Deveria ser 0")
end

main()
