#Neural Networks implementation for optimization

using Optimization, OptimizationFlux, DifferentialEquations, LinearAlgebra,DiffEqFlux

include("dx_dt!.jl")

function optimiseForMass(x0,tspan,solver,p)

  tspan=(tspan[1],round(tspan[2]))
  #Initialising the training data
  ts = collect(0.0:1:tspan[2])
  prob_train = ODEProblem{true}(dx_dt!, x0, tspan, p)
  data_train = Array(solve(prob_train, solver, saveat = ts)) #solver is the ODE solver - Tsit5() e.g.

  #Initialising the neural network
  A = [LegendreBasis(2), LegendreBasis(1)]
  nn = TensorLayer(A, 6)

  f = x -> min(30one(x),x) #Why?


  #We are getting the error that x has length 100 but A has second dimension 12

  function dxdt_pred(dx,x,p,t)
    dx[1] = x[2]
    dx[2] = f(nn(x,p)[1])
    dx[3] = x[4]
    dx[4] = f(nn(x,p)[2])
    dx[5] = x[6]
    dx[6] = f(nn(x,p)[3])
    dx[7] = f(nn(x,p)[4])
    dx[8] = f(nn(x,p)[5])
    dx[9] = f(nn(x,p)[6])
  end

  prob_pred = ODEProblem{true}(dxdt_pred,x0,tspan,p)

  function predict_adjoint(θ)
    x = Array(solve(prob_pred,Tsit5(),p=θ,saveat=ts,
                    sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true))))
  end
  
  function loss_adjoint(θ)
    x = predict_adjoint(θ)
    loss = sum(norm.(x - data_train))
    return loss
  end
  
  iter = 0
  function callback(θ,l)
    global iter
    iter += 1
    if iter%10 == 0
      println(l)
    end
    return false
  end

  α = zeros(length(ts)+1)

  adtype = Optimization.AutoZygote()
  optf = Optimization.OptimizationFunction((x,p) -> loss_adjoint(x), adtype)
  optprob = Optimization.OptimizationProblem(optf, α)
  res1 = Optimization.solve(optprob, ADAM(0.05), callback = callback, maxiters = 150)

  optprob2 = Optimization.OptimizationProblem(optf, res1.u)
  res2 = Optimization.solve(optprob2, ADAM(0.001), callback = callback,maxiters = 150)
  opt = res2.u

  data_pred = predict_adjoint(res1.u)
  return data_pred

end

#=

f = x -> min(30one(x),x)

function dxdt_pred(du,u,p,t)
  du[1] = u[2]
  du[2] = -p[1]*u[1] - p[2]*u[2] + f(nn(u,p[3:end])[1]) #Some kind of linearisation?
end

#Essentially having an actual trained solution, comparing it to a linearised predictor and then defining a loss function. Deriving over the loss function to optimise.

α = zeros(102)

prob_pred = ODEProblem{true}(dxdt_pred,u0,tspan)

function predict_adjoint(θ)
    x = Array(solve(prob_pred,Tsit5(),p=θ,saveat=ts,
                    sensealg=InterpolatingAdjoint(autojacvec=ReverseDiffVJP(true))))
  end
  
  function loss_adjoint(θ)
    x = predict_adjoint(θ)
    loss = sum(norm.(x - data_train))
    return loss
  end
  
  iter = 0
  function callback(θ,l)
    global iter
    iter += 1
    if iter%10 == 0
      println(l)
    end
    return false
  end

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p) -> loss_adjoint(x), adtype)
optprob = Optimization.OptimizationProblem(optf, α)
res1 = Optimization.solve(optprob, ADAM(0.05), callback = callback, maxiters = 150)

optprob2 = Optimization.OptimizationProblem(optf, res1.u)
res2 = Optimization.solve(optprob2, ADAM(0.001), callback = callback,maxiters = 150)
opt = res2.u

using Plots
data_pred = predict_adjoint(res1.u)
plot(ts, data_train[1,:], label = "X (ODE)")
plot!(ts, data_train[2,:], label = "V (ODE)")
plot!(ts, data_pred[1,:], label = "X (NN)")
plot!(ts, data_pred[2,:],label = "V (NN)")

end

=#
