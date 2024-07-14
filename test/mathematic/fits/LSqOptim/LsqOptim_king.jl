using LeastSquaresOptim

"""
  ∑ₖ nₖ = 1

  kwargs:
    autodiff::Symbol ∈ [:forward, :central]
    show_trace::Bool ∈ [true, false]
    factorMethod::Symbol ∈ [:QR, :Cholesky]

"""
show_trace = true
show_trace = false


nModel = nMod[isp3]
# Ml0 = MsnnE[1:3nModel]
Ml0 = MsnnE3[1:3nModel,1]
Ml1 = MsnnE3[1:3nModel,2]

function king!(out, x; Ml0=Ml0, l=0, nModel=nModel)

    nh = x[1:3:end]
    uh = x[2:3:end]
    vhth = x[3:3:end]
    # j = l + 0
    nj = 1
    out[nj] = sum(nh) - Ml0[nj]
    j = l + 2
    nj += 1
    out[nj] = sum((nh .* vhth .^ j) .* (1 .+ 2 / 3 * ((uh ./ vhth) .^ j))) - Ml0[nj]
    j = l + 4
    nj += 1
    out[nj] = sum((nh .* vhth .^ j) .* (1 .+ 4 / 3 * ((uh ./ vhth) .^ 2) .+ 4 / 15 * ((uh ./ vhth) .^ j))) - Ml0[nj]
    for kM in 2:nModel
        for i in 1:3
            j = l + 6(kM-1) + 2(i-1)
            nj += 1
            # N = (j - l) / 2 |> Int
            N = 3(kM-1) + (i-1)
            k = 1:N
            ck = [2^k * binomial(N, k) / prod(3:2:2k+1) for k in 1:N]
            out[nj] = sum((nh .* vhth .^ j) .* (1 .+ [sum(ck .* (uh[s] ./ vhth[s]) .^ (2k)) for s in 1:nModel])) - Ml0[nj]
        end
    end
end

function king_g!(J, x; Ml0=Ml0, l=0, nModel=nModel)

    fill!(J, 0.0)
    nh = x[1:3:end]
    uh = x[2:3:end]
    vhth = x[3:3:end]
    # j = l + 0
    nj = 1
    for s in 1:nModel
        J[nj, 3(s-1)+1] = 1.0
    end
    j = l + 2
    nj += 1
    for s in 1:nModel
        s3 = 3(s-1)
        J[nj, s3+1] = vhth[s]^2 + 2 / 3 * uh[s]^2
        J[nj, s3+2] = 4 / 3 * nh[s] * uh[s]
        J[nj, s3+3] = 2nh[s] * vhth[s]
    end
    j = l + 4 
    nj += 1
    for s in 1:nModel
        s3 = 3(s-1)
        J[nj, s3+1] = vhth[s]^4 + 4 / 3 * uh[s]^2 * vhth[s]^2 + 4 / 15 * uh[s]^4
        J[nj, s3+2] = 8 / 3 * nh[s] * uh[s] * (vhth[s]^2 + 2 / 5 * uh[s]^2)
        J[nj, s3+3] = 4nh[s] * vhth[s] * (vhth[s]^2 + 2 / 3 * uh[s]^2)
    end
    for kM in 2:nModel
        for i in 1:3
            j = l + 6(kM-1) + 2(i-1)
            nj += 1
            N = (j - l) / 2 |> Int
            k = 1:N
            ck = [2^k * binomial(N, k) / prod(3:2:2k+1) for k in 1:N]
            for s in 1:nModel
                s3 = 3(s-1)
                J[nj, s3+1] = vhth[s]^j * (1 + sum(ck .* (uh[s] / vhth[s]) .^ (2k)))
                J[nj, s3+2] = 2nh[s] * vhth[s]^(j - 1) * sum(k .* ck .* (uh[s] / vhth[s]) .^ (2k .- 1))
                J[nj, s3+3] = j * nh[s] * vhth[s]^(j - 1) * (1 + sum((1 .- 2 / j * k) .* ck .* (uh[s] / vhth[s]) .^ (2k)))
            end
        end
    end
end


# The initial solution
p_noise_rel = 1e-6      # `p_noise_rel ≤ 1e-3` default. Relative error of the initial guess parameters `p`
                        # 
p_noise_abs = 1e-15     # Absolute error of the initial guess parameters `p`.
p_noise_ratio = rand(3nModel) * p_noise_rel
x0 = zeros(3nModel)
x0[1:3:end] .= nai0[isp3]
x0[2:3:end] .= uai0[isp3]
x0[3:3:end] .= vthi0[isp3]
x0 .*= (1.0 .+ p_noise_ratio ) + p_noise_abs * rand(3nModel)
# x0[2:3:end] .= epsT1000
@show x0

# The parameter limits for MCF plasma.
lbs = zeros(3nModel)
lbs[1:3:end] .= 0.0
lbs[2:3:end] .= -uhMax
lbs[3:3:end] .= 1/vhthMax
ubs = zeros(3nModel)
ubs[1:3:end] .= nhMax
ubs[2:3:end] .= uhMax
ubs[3:3:end] .= vhthMax

p_tol=1e-37
f_tol=1e-37
g_tol=1e-37
maxIterKing = 1000
(optimizer, optimizer_abbr) = (LeastSquaresOptim.Dogleg, :dl)
# (optimizer, optimizer_abbr) = (LeastSquaresOptim.LevenbergMarquardt, :lm)
factor = LeastSquaresOptim.QR()
# factor = LeastSquaresOptim.Cholesky()
# factor = LeastSquaresOptim.LSMR()
autodiffs = :central
# autodiffs = :forward


# nls = LeastSquaresProblem(x=x0,f! = king!,output_length=length(x0),autodiff=:central)
nls = LeastSquaresProblem(x=x0,f! = king!,g! = king_g!,output_length=length(x0),autodiff=:central)
# res = optimize!(nls, Dogleg(),iterations=maxIterKing,show_trace=show_trace,
#                 x_tol=p_tol,f_tol=f_tol,g_tol=g_tol,lower=lbs,upper=ubs)
# res = optimize!(nls,optimizer(factor),iterations=maxIterKing,show_trace=show_trace,
#                 x_tol=p_tol,f_tol=f_tol,g_tol=g_tol,lower=lbs,upper=ubs)
res = optimize!(nls,optimizer(),iterations=maxIterKing,show_trace=show_trace,
                x_tol=p_tol,f_tol=f_tol,g_tol=g_tol,lower=lbs,upper=ubs)

xfit = res.minimizer         # the vector of best model1 parameters
niter= res.iterations
is_converged = res.converged
xssr = res.ssr                         # sum(abs2, fcur)
is_xconverged = res.x_converged
is_fconverged = res.f_converged
is_gconverged = res.g_converged
optim_method = res.optimizer

naifit = xfit[1:3:end]
uaifit = xfit[2:3:end]
vthifit = xfit[3:3:end]
datafit = DataFrame(ni=nai0[isp3],nfit=naifit,ui=uai0[isp3],ufit=uaifit,vthi=vthi0[isp3],vthifit=vthifit)

dataerrfit = DataFrame(Dni=nai0[isp3]-naifit,Dui=uai0[isp3]-uaifit,Dvthi=vthi0[isp3]-vthifit)

yfit = zeros(3nModel)
king!(yfit, xfit; Ml0=Ml0,l=0)

Jfit = zeros(3nModel,3nModel)
king_g!(Jfit, xfit; Ml0=Ml0,l=0)

fLn1efit = zero.(fLn1e)
fLn1efit = fvLDMz(fLn1efit,vGe,ℓ,naifit,uaifit,vthifit,nModel)
DfLn1e = fLn1e - fLn1efit

@show x0
@show datafit;
@show dataerrfit;
@show norm(DfLn1e)
@show p_noise_rel,p_noise_abs;
@show is_xconverged,is_fconverged,is_gconverged,norm(yfit);
is_converged,niter,xssr,fmtf12.(xfit)