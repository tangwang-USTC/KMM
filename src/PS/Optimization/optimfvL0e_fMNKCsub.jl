

"""
  When `is_nai_const = false/true`
       `∑(nai) = nh = Const`
       `∑(nai .* uai) = uh = Const`
       `∑(nai .* vthi^2) = Th = Const`

  `nModel1 = nModel - 1`

  strategy: nh9 ≤ nhi, ∀i

  Optimization of the amplitude function `fhl0` by using the `King` functions 
    to finding the optimized parameters `(n̂ₛ,ûₛ,v̂ₜₕₛ) = (nai,uai,vthi)`. 

  The general moments is renormalized as: 
    
    `M̂ⱼₗᵐ*`= M̂ⱼₗᵐ / CjLL2(j,L)`.
 
  Notes: `{M̂₁}/3 = Î ≠ û`, generally. Only when `nMod = 1` gives `Î = û`.

  Inputs:
    Mhck1: = M̂ⱼₗᵐ*, which is the renormalized general kinetic moments.
  
  Outputs:
    naik,vthik,nModel = fl0king01_fMNKC!(naik1,vthik1,nModelk1,
                naik,vthik,nModelk,DnuTh,Mhck1,njMs,
                Rdtsabk1,Rdtsabk,NL_solve,rtol_DnuTi,kt,tk,dtk;
                optimizer=optimizer,factor=factor,autodiff=autodiff,
                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi,
                is_optim_CnIK=is_optim_CnIK,is_nhnMod_adapt=is_nhnMod_adapt)
  
"""

# [ite,nMod] 
function fl0king01_fMNKC!(naik1::AbstractVector{T}, vthik1::AbstractVector{T}, nModelk1::Int64,
    LMk::Int64,naik::AbstractVector{T}, vthik::AbstractVector{T}, nModelk::Int64,
    DnuTh::AbstractVector{T}, Mhck1::Matrix{T}, njMs::Int64, edtnIKTs::AbstractVector{T},
    Rdtsabk1::T,Rdtsabk::T, NL_solve::Symbol, rtol_DnuTi::T,kt::Int64,tk::T,dtk::T; 
    optimizer=Dogleg, factor=QR(), autodiff::Symbol=:central,
    is_Jacobian::Bool=true, show_trace::Bool=false, maxIterKing::Int64=200,
    p_tol::Float64=epsT, f_tol::Float64=epsT, g_tol::Float64=epsT, 
    Nspan_optim_nuTi::AbstractVector{T}=[1.1,1.1,1.1],
    is_optim_CnIK::Bool=false,is_nhnMod_adapt::Bool=false) where {T}
    
    # Mode recognition
    if LMk ≤ 1
        if sum(abs.(Mhck1[:,1] .- 1)) / njMs ≤ atol_d1Mhj0
            Mode_Mhck1 = :KM
        else
            Mode_Mhck1 = :KMM
        end
    
        # if Mode_Mhck1 == :KMM
        #     @show Mhck1[:,1] .- 1
        # else
        #     @show Mhck1[:,1] .- 1
        # end
    
        # Tendency according to `Rdtsabk1` and `Rdtsabk`
        if Rdtsabk1 ≤ atol_fM_Rdtsab
            if Rdtsabk ≤ Rdtsabk1               # Deviation from equilibrium,`Mhj0c → ∞`
                Mode_Rdtsab = :KM9
            else                                # Tending to be in equilibrium,`Mhj0c → 1`
                Mode_Rdtsab = :KM1
            end
        else
            if Rdtsabk ≤ Rdtsabk1               # Deviation from equilibrium,`Mhj0c → ∞`
                Mode_Rdtsab = :KMM9
            else                                # Tending to be in equilibrium,`Mhj0c → 1`
                Mode_Rdtsab = :KMM1
            end
        end

        # @show kt,tk, dtk, nModelk
        # @show naik
        # errfff
        # 
        if is_fixed_NK
            nModelk1 = 1NK
            if nModelk1 == 1
                yfitNK = zeros(T,2nModelk1)
                fl0king01_fMNK!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                        optimizer=optimizer,factor=factor,autodiff=autodiff,
                        is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                        p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                
            else
                if is_optim_CnIK
                    yfitNK = zeros(T,2nModelk1-2)
                    fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                            optimizer=optimizer,factor=factor,autodiff=autodiff,
                            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                else
                    yfitNK = zeros(T,2nModelk1)
                    fl0king01_fMNK!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                            optimizer=optimizer,factor=factor,autodiff=autodiff,
                            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                end
            end


            # if nModelk1 == 1
            #     yfitNK = zeros(T,2nModelk1)
            #     fl0king01_fMNKf!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
            #             optimizer=optimizer,factor=factor,autodiff=autodiff,
            #             is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
            #             p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                
            # else
            #     if is_optim_CnIK
            #         yfitNK = zeros(T,2nModelk1-2)
            #         fl0king01_fMNKfC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
            #                 optimizer=optimizer,factor=factor,autodiff=autodiff,
            #                 is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
            #                 p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
            #     else
            #         yfitNK = zeros(T,2nModelk1)
            #         fl0king01_fMNKf!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
            #                 optimizer=optimizer,factor=factor,autodiff=autodiff,
            #                 is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
            #                 p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
            #     end
            # end
            yfitNKM = maximum(abs.(yfitNK))
            if yfitNKM == NaN && yfitNKM == Inf
                @show yfitNKM
                edrf11111
            end
            errNormNK = norm([DnuTh[1], DnuTh[3], yfitNKM])
            if errNormNK ≤ atol_nuTi_optim
                printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3], yfitNKM), color=:green,"\n")
            else
                if yfitNKM ≤ atol_nuTi_optim
                    printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:red,"\n")
                    printstyled("21NK, yfitNKM=", fmtf2.(yfitNKM), color=:green,"\n")
                    if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                        @error("21NK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                    end
                else
                    if norm([DnuTh[1], DnuTh[3]]) ≤ atol_nuTi_optim
                        printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:green,"\n")
                        printstyled("21NK, yfitNK=", fmtf2.(yfitNK), color=:red,"\n")
                    else
                        printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3],yfitNKM), color=:red,"\n")
                        if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                            @error("21NK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                        end
                    end
                end
        
                @error("21NK: Underfitting: Optimization process for `naik,vthik` is not converged when",nModelk1)
            end
        else
            if 12 == 1
                if abs(Rdtsabk) ≤ 1e-4 || tk + dtk ≥ 0.9
                    is_nai_const_inner = true
                else
                    is_nai_const_inner = false
                end
                if Mode_Rdtsab == :KM1                  # Tending to be in equilibrium,`Mhj0c → 1`
                    if nModelk == 1
                        nModelk1 = 2
                    else
                        if nModelk < NK
                            nModelk1 = nModelk + 1
                        else
                            nModelk1 = 1NK
                        end
                    end
                elseif Mode_Rdtsab == :KM9              # Deviation from equilibrium,`Mhj0c → ∞` 
                    if nModelk == 1
                        nModelk1 = 1
                    else
                        nModelk1 = nModelk - 1
                    end
                else
                    if Mode_Rdtsab == :KMM1
                        nModelk1 = 1NK
                    elseif Mode_Rdtsab == :KMM9
                        nModelk1 = 1NK
                    end
                end
            else
            end

            # @show Brag_condition_edtnIKTs, Brag_condition_Rdtsabk, Brag_condition_tk, tk
            if isnan(Brag_condition_edtnIKTs)
                if isnan(Brag_condition_Rdtsabk)
                    if isnan(Brag_condition_tk)
                        sdfgftgd
                    else
                        is_Brag_condition = (tk > Brag_condition_tk)
                        if is_Brag_condition == false
                            is_Brag_condition = norm(edtnIKTs[1:3]) ≥ 1e-8
                        end
                    end
                else
                    is_Brag_condition = abs(Rdtsabk) ≤ Brag_condition_Rdtsabk
                end
            else
                is_Brag_condition = norm(edtnIKTs[1:3]) ≥ Brag_condition_edtnIKTs
            end
            # @show "11111111111112222222222222233333333333333333", is_Brag_condition
            if is_Brag_condition
                is_nai_const_inner = true
            else
                is_nai_const_inner = false
            end

            if is_nai_const_inner
                @show nModelk1
                naik1,vthik1,nModelk1 = fl0king01_fMnhC!(naik1,vthik1,nModelk1,nModelk,
                            DnuTh,Mhck1,NL_solve,rtol_DnuTi,tk,dtk;
                            optimizer=optimizer,factor=factor,autodiff=autodiff,
                            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi,
                            is_optim_CnIK=is_optim_CnIK,is_nhnMod_adapt=is_nhnMod_adapt)
            else
                nModelk1 = 1NK
                if nModelk1 == 1
                    yfitNK = zeros(T,2nModelk1)
                    fl0king01_fMNK!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                            optimizer=optimizer,factor=factor,autodiff=autodiff,
                            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    
                else
                    if is_optim_CnIK
                        yfitNK = zeros(T,2nModelk1-2)
                        fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                                optimizer=optimizer,factor=factor,autodiff=autodiff,
                                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    else
                        yfitNK = zeros(T,2nModelk1)
                        fl0king01_fMNK!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                                optimizer=optimizer,factor=factor,autodiff=autodiff,
                                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    end
                end
                yfitNKM = maximum(abs.(yfitNK))
                if yfitNKM == NaN && yfitNKM == Inf
                    @show yfitNKM
                    edrf11111
                end
                errNormNK = norm([DnuTh[1], DnuTh[3], yfitNKM])
                if errNormNK ≤ atol_nuTi_optim
                    printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3], yfitNKM), color=:green,"\n")
                else
                    if yfitNKM ≤ atol_nuTi_optim
                        printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:red,"\n")
                        printstyled("21NK, yfitNKM=", fmtf2.(yfitNKM), color=:green,"\n")
                        if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                            @error("21NK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                        end
                    else
                        if norm([DnuTh[1], DnuTh[3]]) ≤ atol_nuTi_optim
                            printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:green,"\n")
                            printstyled("21NK, yfitNK=", fmtf2.(yfitNK), color=:red,"\n")
                        else
                            printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3],yfitNKM), color=:red,"\n")
                            if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                                @error("21NK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                            end
                        end
                    end
            
                    @error("21NK: Underfitting: Optimization process for `naik,vthik` is not converged when",nModelk1)
                end
            end
        end
    else
        edfgvbn
        if Rdtsabk ≤ Rdtsabk1               # Deviation from equilibrium,`Mhj0c → ∞`
            Mode_Rdtsab = :KMDM9
        else                                # Tending to be in equilibrium,`Mhj0c → 1`
            Mode_Rdtsab = :KMDM1
        end
        # 
        if Mode_Mhck1 == :KM
            nModelk1 = 1NK
            if nModelk1 == 1
                yfitNK = zeros(T,2nModelk1)
                fl0king01_fMNK!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                        optimizer=optimizer,factor=factor,autodiff=autodiff,
                        is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                        p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
            else
                if is_fixed_NK
                    if is_optim_CnIK
                        yfitNK = zeros(T,2nModelk1-2)
                        fl0king01_fMNKfC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                                optimizer=optimizer,factor=factor,autodiff=autodiff,
                                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    else
                        yfitNK = zeros(T,2nModelk1)
                        fl0king01_fMNKfC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                                optimizer=optimizer,factor=factor,autodiff=autodiff,
                                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    end
                else
                    if is_optim_CnIK
                        yfitNK = zeros(T,2nModelk1-2)
                        fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                                optimizer=optimizer,factor=factor,autodiff=autodiff,
                                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    else
                        yfitNK = zeros(T,2nModelk1)
                        fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                                optimizer=optimizer,factor=factor,autodiff=autodiff,
                                is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                                p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                    end
                end
            end
            yfitNKM = maximum(abs.(yfitNK))
            if yfitNKM == NaN && yfitNKM == Inf
                @show yfitNKM
                edrf11111
            end
            errNormNK = norm([DnuTh[1], DnuTh[3], yfitNKM])
            if errNormNK ≤ atol_nuTi_optim
                printstyled("11KMNK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3], yfitNKM), color=:green,"\n")
            else
                if yfitNKM ≤ atol_nuTi_optim
                    printstyled("12KMNK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:red,"\n")
                    printstyled("12KMNK, yfitNKM=", fmtf2.(yfitNKM), color=:green,"\n")
                    if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                        @error("12KMNK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                    end
                else
                    if norm([DnuTh[1], DnuTh[3]]) ≤ atol_nuTi_optim
                        printstyled("12KMNK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:green,"\n")
                        printstyled("12KMNK, yfitNK=", fmtf2.(yfitNK), color=:red,"\n")
                    else
                        printstyled("12KMNK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3],yfitNKM), color=:red,"\n")
                        if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                            @error("12KMNK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                        end
                    end
                end
        
                @error("11KMNK: Underfitting: Optimization process for `naik,vthik` is not converged when",nModelk1)
            end
        else
            nModelk1 = 1NK
            if nModelk1 == 1
                yfitNK = zeros(T,2nModelk1)
                fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                        optimizer=optimizer,factor=factor,autodiff=autodiff,
                        is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                        p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
            else
                if is_optim_CnIK
                    yfitNK = zeros(T,2nModelk1-2)
                    fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                            optimizer=optimizer,factor=factor,autodiff=autodiff,
                            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                else
                    yfitNK = zeros(T,2nModelk1)
                    fl0king01_fMNKC!(naik1,vthik1,nModelk1,DnuTh,yfitNK,Mhck1,NL_solve;
                            optimizer=optimizer,factor=factor,autodiff=autodiff,
                            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
                            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
                end
            end
            yfitNKM = maximum(abs.(yfitNK))
            if yfitNKM == NaN && yfitNKM == Inf
                @show yfitNKM
                edrf11111
            end
            errNormNK = norm([DnuTh[1], DnuTh[3], yfitNKM])
            if errNormNK ≤ atol_nuTi_optim
                printstyled("21NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3], yfitNKM), color=:green,"\n")
            else
                if yfitNKM ≤ atol_nuTi_optim
                    printstyled("22NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:red,"\n")
                    printstyled("22NK, yfitNKM=", fmtf2.(yfitNKM), color=:green,"\n")
                    if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                        @error("22NK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                    end
                else
                    if norm([DnuTh[1], DnuTh[3]]) ≤ atol_nuTi_optim
                        printstyled("22NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]), color=:green,"\n")
                        printstyled("22NK, yfitNK=", fmtf2.(yfitNK), color=:red,"\n")
                    else
                        printstyled("22NK, (Dnh,DTh)=", (DnuTh[1], DnuTh[3],yfitNKM), color=:red,"\n")
                        if is_optim_CnIK && yfitNKM ≥ atol_nuTi_optim_err
                            @error("22NK, Optimization is falure, (Dnh,DTh)=", (DnuTh[1], DnuTh[3]))
                        end
                    end
                end
        
                @error("21NK: Underfitting: Optimization process for `naik,vthik` is not converged when",nModelk1)
            end
        end
    end
    
    return naik1, vthik1, nModelk1
end

"""
  Inputs:
    Mhck1: = M̂ⱼₗᵐ*, which is the renormalized general kinetic moments.
  
  Outputs:
    fl0king01_fMNKC!(naik,vthik,nModel,DnuTh,yfit,Mhck1,NL_solve;
            optimizer=optimizer,factor=factor,autodiff=autodiff,
            is_Jacobian=is_Jacobian,show_trace=show_trace,maxIterKing=maxIterKing,
            p_tol=p_tol,f_tol=f_tol,g_tol=g_tol,Nspan_optim_nuTi=Nspan_optim_nuTi)
  
"""

# [nMod]
function fl0king01_fMNKC!(naik::AbstractVector{T}, vthik::AbstractVector{T}, nModel::Int64,
    DnuTh::AbstractVector{T}, yfit::AbstractVector{T}, Mhck1::Matrix{T}, NL_solve::Symbol; 
    optimizer=Dogleg, factor=QR(), autodiff::Symbol=:central,
    is_Jacobian::Bool=true, show_trace::Bool=false, maxIterKing::Int64=200,
    p_tol::Float64=epsT, f_tol::Float64=epsT, g_tol::Float64=epsT, 
    Nspan_optim_nuTi::AbstractVector{T}=[1.1,1.1,1.1]) where {T}

    # The parameter limits for MCF plasma.
    nModel2 = 2nModel
    nModel1 = nModel - 1
    nModel12 = 2nModel1
    x0 = zeros(nModel12)      # [nai1, vthi1, nai2, vthi2, ⋯]
    lbs = zeros(nModel12)
    ubs = zeros(nModel12)

    vecMod = 1:nModel
    vecMod1 = 1:nModel1

    # nai
    lbs[1:2:nModel12] .= 0.0
    ubs[1:2:nModel12] .= 1.0
    ############################################### I
    # x0[1:2:nModel12] = deepcopy(naik[vecMod1])
    ############################################### II
    # for i in vecMod1
    #     x0[2i-1] = deepcopy(naik[i])
    # end
    ############################################### III
    # x0[1] = deepcopy(naik[1])
    # # x0[3:2:nModel12] .= n0_c
    ############################################### IV
    # x0[1] = deepcopy(naik[1])
    # x0[3:2:nModel12] .= 0.1
    ############################################### IIV
    x0[1:2:nModel12] .= nhInitial

    # # vthi
    # for i in vecMod1
    #     ubs[2i] = min(vhthMax, Nspan_optim_nuTi[3] * vthik[i])
    #     lbs[2i] = max(vhthMin, vthik[i] / Nspan_optim_nuTi[3])
    # end
    lbs[2:2:nModel12] .= vhthMin
    ubs[2:2:nModel12] .= vhthMax
    if is_NKC_vhthInitial
        for i in vecMod1
            x0[2i] = max(vhthInitial,lbs[2i])
        end
    else
        for i in vecMod1
            x0[2i] = max(vthik[i] * vhthRatio,lbs[2i])
        end
    end

    vhth,vhth2 = zeros(T,nModel),zeros(T,nModel)
    res = fl0king01_fMNKC(deepcopy(x0), naik[vecMod], vthik[vecMod], nModel1, 
        Mhck1[1:nModel2,1], NL_solve; 
        vhth2=vhth2, lbs=lbs, ubs=ubs,
        optimizer=optimizer, factor=factor, autodiff=autodiff,
        is_Jacobian=is_Jacobian, show_trace=show_trace, maxIterKing=maxIterKing,
        p_tol=p_tol, f_tol=f_tol, g_tol=g_tol, NL_solve_method=NL_solve_method)

    if NL_solve == :LeastSquaresOptim
        xfit = res.minimizer         # the vector of best model1 parameters
        naik[vecMod1] = xfit[1:2:nModel12]          # the vector of best model1 parameters
        vthik[vecMod1] = xfit[2:2:nModel12]         # the vector of best model1 parameters
        niter = res.iterations
        is_converged = res.converged
        xssr = res.ssr                         # sum(abs2, fcur)
    elseif NL_solve == :NLsolve
        xfit = res.zero         # the vector of best model1 parameters
        naik[vecMod1] = xfit[1:2:nModel12]          # the vector of best model1 parameters
        vthik[vecMod1] = xfit[2:2:nModel12]
        niter = res.iterations
        is_converged = res.f_converged
        xssr = res.residual_norm                         # sum(abs2, fcur)
    elseif NL_solve == :JuMP
        fgfgg
    end

    if nModel == 2
        naik[nModel], vthik[nModel] = king_fMNKC2!(yfit, xfit, naik[vecMod], Mhck1[1:nModel2,1];vhth=vhth,vhth2=vhth2)
    else
        naik[nModel], vthik[nModel] = king_fMNKC3!(yfit, xfit, naik[vecMod], nModel1, Mhck1[1:nModel2,1];vhth=vhth,vhth2=vhth2)
    end

    DnuTh[1] = sum(naik[vecMod]) - 1
    # # T̂ = ∑ₖ (n̂ₖv̂ₜₕₖ²) - 1.
    DnuTh[3] = sum(naik[vecMod] .* vthik[vecMod] .^ 2) - 1                                   # Thfit .- 1
end

"""
  Inputs:
    Mhck1: = M̂ⱼₗᵐ*, which is the renormalized general kinetic moments.
  
  Outputs:
  res = fl0king01_fMNKC(x0, naik, vthik, Mhck1, NL_solve; 
                    vhth=vhth,vhth2=vhth2,lbs=lbs, ubs=ubs,
                    optimizer=optimizer, factor=factor, autodiff=autodiff,
                    is_Jacobian=is_Jacobian, show_trace=show_trace, maxIterKing=maxIterKing,
                    p_tol=p_tol, f_tol=f_tol, g_tol=g_tol, NL_solve_method=NL_solve_method)
  res = fl0king01_fMNKC(x0, naik, vthik, nModel1, Mhck1, NL_solve; 
                    vhth=vhth,vhth2=vhth2,lbs=lbs, ubs=ubs,
                    optimizer=optimizer, factor=factor, autodiff=autodiff,
                    is_Jacobian=is_Jacobian, show_trace=show_trace, maxIterKing=maxIterKing,
                    p_tol=p_tol, f_tol=f_tol, g_tol=g_tol, NL_solve_method=NL_solve_method)
  
"""

# nMode = 1

# nMod = 2, Optimization
function fl0king01_fMNKC(x0::AbstractVector{T}, nai::AbstractVector{T}, vhth::AbstractVector{T}, 
    Mhck1::AbstractVector{T}, NL_solve::Symbol; vhth2::AbstractVector{T}=[0.1, 1.0],
    lbs::AbstractVector{T}=[-uhMax, 0.8], ubs::AbstractVector{T}=[uhMax, 1.2],
    optimizer=Dogleg, factor=QR(), autodiff::Symbol=:central,
    is_Jacobian::Bool=true, show_trace::Bool=false, maxIterKing::Int64=200,
    p_tol::Float64=epsT, f_tol::Float64=epsT, g_tol::Float64=epsT,
    NL_solve_method::Symbol=:newton) where {T}

    king01_fMNKC2!(out, x) = king_fMNKC2!(out, x, nai;vhth=vhth,vhth2=vhth2,Mhck1=Mhck1)
    if NL_solve == :LeastSquaresOptim
        if is_Jacobian
            J!(J, x) = king_fMNKC2_g!(J, x, nai;vhth=vhth,vhth2=vhth2,Mhck1=Mhck1)
            nls = LeastSquaresProblem(x=x0, (f!)=king01_fMNKC2!, (g!)=J!, output_length=length(x0), autodiff=autodiff)
        else
            nls = LeastSquaresProblem(x=x0, (f!)=king01_fMNKC2!, output_length=length(x0), autodiff=autodiff)
        end
        res = optimize!(nls, optimizer(factor), iterations=maxIterKing, show_trace=show_trace,
            x_tol=p_tol, f_tol=f_tol, g_tol=g_tol, lower=lbs, upper=ubs)
    elseif NL_solve == :NLsolve
        if is_Jacobian
            Js!(J, x) = king_fMNKC2_g!(J, x, nai;vhth=vhth,vhth2=vhth2,Mhck1=Mhck1)
            nls = OnceDifferentiable(king01_fMNKC2!, Js!, x0, similar(x0))
            if NL_solve_method == :trust_region
                res = nlsolve(nls, x0, method=NL_solve_method, factor=1.0, autoscale=true, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace)
            elseif NL_solve_method == :newton
                res = nlsolve(nls, x0, method=NL_solve_method, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace)
            end
        else
            nls = OnceDifferentiable(king01_fMNKC2!, x0, similar(x0))
            if NL_solve_method == :trust_region
                res = nlsolve(nls, x0, method=NL_solve_method, factor=1.0, autoscale=true, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace, autodiff=:forward)
            elseif NL_solve_method == :newton
                res = nlsolve(nls, x0, method=NL_solve_method, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace, autodiff=:forward)
            end
        end
    elseif NL_solve == :JuMP
        gyhhjjmj
    else
        esfgroifrtg
    end
    return res
end

# [nMod ≥ 2], Optimization
function fl0king01_fMNKC(x0::AbstractVector{T}, nai::AbstractVector{T}, vhth::AbstractVector{T}, nModel1::Int64, 
    Mhck1::AbstractVector{T}, NL_solve::Symbol; vhth2::AbstractVector{T}=[0.1, 1.0], 
    lbs::AbstractVector{T}=[-uhMax, 0.8], ubs::AbstractVector{T}=[uhMax, 1.2],
    optimizer=Dogleg, factor=QR(), autodiff::Symbol=:central,
    is_Jacobian::Bool=true, show_trace::Bool=false, maxIterKing::Int64=200,
    p_tol::Float64=epsT, f_tol::Float64=epsT, g_tol::Float64=epsT,
    NL_solve_method::Symbol=:newton) where {T}

    king01_fMNKC3!(out, x) = king_fMNKC3!(out, x, nai;vhth=vhth,vhth2=vhth2,nModel1=nModel1,Mhck1=Mhck1)
    if NL_solve == :LeastSquaresOptim
        if is_Jacobian
            J!(J, x) = king_fMNKC3_g!(J, x, nai;vhth=vhth,vhth2=vhth2,nModel1=nModel1,Mhck1=Mhck1)
            nls = LeastSquaresProblem(x=x0, (f!)=king01_fMNKC3!, (g!)=J!, output_length=length(x0), autodiff=autodiff)
        else
            nls = LeastSquaresProblem(x=x0, (f!)=king01_fMNKC3!, output_length=length(x0), autodiff=autodiff)
        end
        res = optimize!(nls, optimizer(factor), iterations=maxIterKing, show_trace=show_trace,
            x_tol=p_tol, f_tol=f_tol, g_tol=g_tol, lower=lbs, upper=ubs)
    elseif NL_solve == :NLsolve
        if is_Jacobian
            Js!(J, x) = king_fMNKC3_g!(J, x, nai;vhth=vhth,vhth2=vhth2,nModel1=nModel1,Mhck1=Mhck1)
            nls = OnceDifferentiable(king01_fMNKC3!, Js!, x0, similar(x0))
            if NL_solve_method == :trust_region
                res = nlsolve(nls, x0, method=NL_solve_method, factor=1.0, autoscale=true, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace)
            elseif NL_solve_method == :newton
                res = nlsolve(nls, x0, method=NL_solve_method, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace)
            end
        else
            nls = OnceDifferentiable(king01_fMNKC3!, x0, similar(x0))
            if NL_solve_method == :trust_region
                res = nlsolve(nls, x0, method=NL_solve_method, factor=1.0, autoscale=true, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace, autodiff=:forward)
            elseif NL_solve_method == :newton
                res = nlsolve(nls, x0, method=NL_solve_method, xtol=p_tol, ftol=f_tol,
                    iterations=maxIterKing, show_trace=show_trace, autodiff=:forward)
            end
        end
    elseif NL_solve == :JuMP
        gyhhjjmj
    else
        esfgroifrtg
    end
    return res
end

"""
  Inputs:
    out: = zeros(nModel1,nModel1)
    x: = x(nModel1)
    nh: = nai
    vhth: = vthi
    nModel1 = nModel - 1
    x: = x(nModel1)
    Mhck1: = M̂ⱼₗᵐ*, which is the renormalized general kinetic moments.

  Outputs:
    king_fMNKC!(out, x, nh;vhth=vhth, vhth2=vhth2, Mhck1=Mhck1)

"""

# When `is_nai_const == true`
#      `∑(nai .* vthi^2 = Mh(1,1) = Kh = Th`
# nMode = 1

# nMode = 2
function king_fMNKC2!(out::AbstractVector{T}, x::AbstractVector{T}, nh::AbstractVector{T}; 
    vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0], Mhck1::AbstractVector{T}=[0.1, 1.0]) where{T}

    nh[1] = x[1]
    vhth[1] = x[2]

    vhth2[1] = vhth[1] ^ 2

    nh[2] = Mhck1[1] - nh[1]

    if nh[2] ≥ n0_c
        vhth2[2] = (Mhck1[2] - (nh[1] .* vhth2[1])) / nh[2]
        vhth[2] = vhth2[2]^0.5

        nj = 1
        # # (l,j) = (0,4)
        out[nj] = sum(nh .* vhth2 .^ 2) - Mhck1[nj+2]
    
        nj += 1
        # # (l,j) = (0,6)
        out[nj] = sum(nh .* vhth2 .^ 3) - Mhck1[nj+2]
    else
        nh[2] = 0.0
        vhth2[2] = vhthMin^2
        vhth[2] = vhthMin

        nj = 1
        # # (l,j) = (0,4)
        out[nj] = nh[1] .* vhth2[1] .^ 2 - Mhck1[nj+2]
    
        nj += 1
        # # (l,j) = (0,6)
        out[nj] = nh[1] .* vhth2[1] .^ 3 - Mhck1[nj+2]
    end
end

function king_fMNKC2!(out::AbstractVector{T}, x::AbstractVector{T}, nh::AbstractVector{T}, Mhck1::AbstractVector{T}; 
    vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0]) where{T}

    nh[1] = x[1]
    vhth[1] = x[2]

    vhth2[1] = vhth[1] ^ 2

    nh[2] = Mhck1[1] - nh[1]

    if nh[2] ≥ n0_c
        vhth2[2] = (Mhck1[2] - (nh[1] .* vhth2[1])) / nh[2]
        vhth[2] = vhth2[2]^0.5

        nj = 1
        # # (l,j) = (0,4)
        out[nj] = sum(nh .* vhth2 .^ 2) - Mhck1[nj+2]
    
        nj += 1
        # # (l,j) = (0,6)
        out[nj] = sum(nh .* vhth2 .^ 3) - Mhck1[nj+2]
    else
        nh[2] = 0.0
        vhth2[2] = vhthMin^2
        vhth[2] = vhthMin

        nj = 1
        # # (l,j) = (0,4)
        out[nj] = nh[1] .* vhth2[1] .^ 2 - Mhck1[nj+2]
    
        nj += 1
        # # (l,j) = (0,6)
        out[nj] = nh[1] .* vhth2[1] .^ 3 - Mhck1[nj+2]
    end
    return nh[2], vhth[2]
end

# nMode - 1 ≥ 1
function king_fMNKC3!(out::AbstractVector{T}, x::AbstractVector{T}, nh::AbstractVector{T};
    vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0], nModel1::Int64=1, 
    Mhck1::AbstractVector{T}=[0.1, 1.0]) where{T}

    nh[1:nModel1] = x[1:2:end]
    vhth[1:nModel1] = x[2:2:end]

    vhth2[1:end-1] = vhth[1:end-1] .^ 2

    nh[end] = Mhck1[1] - sum(nh[1:end-1])

    # if nh[end] ≥ n0_c
        Kh9 = (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1]))
        vhth2[end] = Kh9 / nh[end]
        vhth[end] = vhth2[end]^0.5
    # else
    #     nh[end] = 0.0
    #     vhth2[end] = vhthMin^2
    #     vhth[end] = vhthMin
    # end

    nj = 1
    # # (l,j) = (0,4)
    out[nj] = sum(nh .* vhth2 .^ 2) - Mhck1[nj+2]

    nj += 1
    # # (l,j) = (0,6)
    out[nj] = sum(nh .* vhth2 .^ 3) - Mhck1[nj+2]

    for kM in 2:nModel1
        nj += 1
        # (l, j) = (0, 2(nj+1))
        out[nj] = sum(nh .* vhth2 .^ (nj+1)) - Mhck1[nj+2]

        nj += 1
        # (l, j) = (0, 2(nj+1))
        out[nj] = sum(nh .* vhth2 .^ (nj+1)) - Mhck1[nj+2]
    end
end

function king_fMNKC3!(out::AbstractVector{T}, x::AbstractVector{T}, nh::AbstractVector{T}, nModel1::Int64, 
    Mhck1::AbstractVector{T}; vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0]) where{T}

    nh[1:nModel1] = x[1:2:end]
    vhth[1:nModel1] = x[2:2:end]

    vhth2[1:end-1] = vhth[1:end-1] .^ 2

    nh[end] = Mhck1[1] - sum(nh[1:end-1])

    # if nh[end] ≥ n0_c
        vhth2[end] = (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1])) / nh[end]
        vhth[end] = vhth2[end]^0.5
    # else
    #     nh[end] = 0.0
    #     vhth2[end] = vhthMin^2
    #     vhth[end] = vhthMin
    # end

    nj = 1
    # # (l,j) = (0,4)
    out[nj] = sum(nh .* vhth2 .^ 2) - Mhck1[nj+2]

    nj += 1
    # # (l,j) = (0,6)
    out[nj] = sum(nh .* vhth2 .^ 3) - Mhck1[nj+2]

    for kM in 2:nModel1
        nj += 1
        # (l, j) = (0, 2(nj+1))
        out[nj] = sum(nh .* vhth2 .^ (nj+1)) - Mhck1[nj+2]

        nj += 1
        # (l, j) = (0, 2(nj+1))
        out[nj] = sum(nh .* vhth2 .^ (nj+1)) - Mhck1[nj+2]
    end
    return nh[end], vhth[end]
end

"""
  Inputs:
    J: = zeros(nModel1,nModel1)
    x: = x(nModel1)
    nh: = nai
    vhth: = vthi[1:nModel1]
    nModel1 = nModel - 1

  Outputs:
    king_fMNKC2_g!(J, x, nh;vhth=vhth,vhth2=vhth2, Mhck1)
    king_fMNKC3_g!(J, x, nh, nModel1;vhth=vhth,vhth2=vhth2, Mhck1)
"""

# The Jacobian matrix: J = zeros(T,nMod-1,nMod-1)
# nMode = 2
function king_fMNKC2_g!(J::AbstractArray{T,N}, x::AbstractVector{T}, nh::AbstractVector{T}; 
    vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0], Mhck1::AbstractVector{T}=[0.1, 1.0]) where{T,N}

    fill!(J, 0.0)
    nh[1] = x[1]
    vhth[1] = x[1]

    vhth2[1] = vhth[1] ^ 2

    nh[2] = Mhck1[1] - nh[1]

    # if nh[2] ≥ n0_c
        vhth2[2] = (Mhck1[2] - nh[1] .* vhth2[1]) / nh[2]
        vhth[2] = vhth2[2]^0.5

        nj = 1
        # (l, j) = (0, 4)
        ration = (1 + nh[nj] / nh[2])
        J[nj, 1] = vhth2[nj]^2 - vhth2[2]^2 - 2vhth2[nj] * vhth2[2] * ration
        J[nj, 2] = 4nh[nj] * vhth[nj] * (vhth2[nj] - vhth2[2])
    
        nj += 1
        # (l, j) = (0, 6)   # j = 2(nj+1)
        s = 1
        ration = (1 + nh[nj] / nh[2])
        J[nj, s] = vhth2[nj]^3 - vhth2[2]^3 - 3vhth2[nj] * vhth2[2]^2 * ration
        J[nj, s+1] = 6nh[nj] * vhth[nj] * (vhth2[nj]^2 - vhth2[2]^2)
    # else
    #     nh[2] = 0.0
    #     vhth2[2] = vhthMin^2
    #     vhth[2] = vhthMin

    #     nj = 1
    #     # (l, j) = (0, 4)
    #     J[nj, 1] = vhth2[nj]^2
    #     J[nj, 2] = 4nh[nj] * vhth[nj] * vhth2[nj]
    
    #     nj += 1
    #     # (l, j) = (0, 6)   # j = 2(nj+1)
    #     s = 1
    #     J[nj, s] = vhth2[nj]^3
    #     J[nj, s+1] = 6nh[nj] * vhth[nj] * vhth2[nj]^2
    # end
end

function king_fMNKC2_g!(J::AbstractArray{T,N}, x::AbstractVector{T}, nh::AbstractVector{T}, Mhck1::AbstractVector{T}; 
    vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0]) where{T,N}

    fill!(J, 0.0)
    nh[1] = x[1]
    vhth[1] = x[1]

    vhth2[1] = vhth[1] ^ 2

    nh[2] = Mhck1[1] - nh[1]

    # if nh[2] ≥ n0_c
        vhth2[2] = (Mhck1[2] - nh[1] .* vhth2[1]) / nh[2]
        vhth[2] = vhth2[2]^0.5

        nj = 1
        # (l, j) = (0, 4)
        ration = (1 + nh[nj] / nh[2])
        J[nj, 1] = vhth2[nj]^2 - vhth2[2]^2 - 2vhth2[nj] * vhth2[2] * ration
        J[nj, 2] = 4nh[nj] * vhth[nj] * (vhth2[nj] - vhth2[2])
    
        nj += 1
        # (l, j) = (0, 6)   # j = 2(nj+1)
        s = 1
        ration = (1 + nh[nj] / nh[2])
        J[nj, s] = vhth2[nj]^3 - vhth2[2]^3 - 3vhth2[nj] * vhth2[2]^2 * ration
        J[nj, s+1] = 6nh[nj] * vhth[nj] * (vhth2[nj]^2 - vhth2[2]^2)
    # else
    #     nh[2] = 0.0
    #     vhth2[2] = vhthMin^2
    #     vhth[2] = vhthMin

    #     nj = 1
    #     # (l, j) = (0, 4)
    #     J[nj, 1] = vhth2[nj]^2
    #     J[nj, 2] = 4nh[nj] * vhth[nj] * vhth2[nj]
    
    #     nj += 1
    #     # (l, j) = (0, 6)   # j = 2(nj+1)
    #     s = 1
    #     J[nj, s] = vhth2[nj]^3
    #     J[nj, s+1] = 6nh[nj] * vhth[nj] * vhth2[nj]^2
    # end
    return nh[2], vhth[2]
end

# NK ≥ 2
function king_fMNKC3_g!(J::AbstractArray{T,N}, x::AbstractVector{T}, nh::AbstractVector{T}; 
    vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0], nModel1::Int64=2, Mhck1::AbstractVector{T}=[0.1, 1.0]) where{T,N}

    fill!(J, 0.0)
    nh[1:nModel1] = x[1:2:end]
    vhth[1:nModel1] = x[2:2:end]

    vhth2[1:end-1] = vhth[1:end-1] .^ 2

    nh[end] = Mhck1[1] - sum(nh[1:end-1])

    # if nh[end] ≥ n0_c
        vhth2[end] = (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1])) / nh[end]
        vhth[end] = vhth2[end]^0.5
        
        nj = 1
        # (l, j) = (0, 4)   # j = 2(nj+1)
        s = 1
        ration = (1 + nh[nj] / nh[end])
        J[nj, s] = vhth2[nj]^2 - vhth2[end]^2 - 2vhth2[nj] * vhth2[end] * ration
        J[nj, s+1] = 4nh[nj] * vhth[nj] * (vhth2[nj] - vhth2[end])
    
        nj += 1
        # (l, j) = (0, 6)   # j = 2(nj+1)
        s = 1
        ration = (1 + nh[nj] / nh[end])
        J[nj, s] = vhth2[nj]^3 - vhth2[end]^3 - 3vhth2[nj] * vhth2[end]^2 * ration
        J[nj, s+1] = 6nh[nj] * vhth[nj] * (vhth2[nj]^2 - vhth2[end]^2)
        
        # (l, j) = (0, nj)
        for kM in 2:nModel1
            nj += 1
            ration = (1 + nh[nj] / nh[end])
            for s in 1:2
                s2 = 2(s - 1)
                # j = 2(nj + 1)
                j2 = (nj + 1)
                J[nj, s2+1] = vhth2[nj]^j2 - vhth2[end]^j2 - j2 * vhth2[nj] * vhth2[end]^(j2-1) * ration
                J[nj, s2+2] = 2j2 * nh[nj] * vhth[nj] * (vhth2[nj]^(j2-1) - vhth2[end]^(j2-1))
            end
        end
    # else
    #     @show vhth2
    #     @show (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1]))
    #     nh[end] = 0.0
    #     vhth2[end] = vhthMin^2
    #     vhth[end] = vhthMin
        
    #     nj = 1
    #     # (l, j) = (0, 4)   # j = 2(nj+1)
    #     s = 1
    #     J[nj, s] = vhth2[nj]^2
    #     J[nj, s+1] = 4nh[nj] * vhth[nj] * vhth2[nj]
    
    #     nj += 1
    #     # (l, j) = (0, 6)   # j = 2(nj+1)
    #     s = 1
    #     J[nj, s] = vhth2[nj]^3
    #     J[nj, s+1] = 6nh[nj] * vhth[nj] * vhth2[nj]^2
        
    #     # (l, j) = (0, nj)
    #     for kM in 2:nModel1
    #         nj += 1
    #         for s in 1:2
    #             s2 = 2(s - 1)
    #             # j = 2(nj + 1)
    #             j2 = (nj + 1)
    #             J[nj, s2+1] = vhth2[nj]^j2
    #             J[nj, s2+2] = 2j2 * nh[nj] * vhth[nj] * vhth2[nj]^(j2-1)
    #         end
    #     end
    # end
end

function king_fMNKC3_g!(J::AbstractArray{T,N}, x::AbstractVector{T}, nh::AbstractVector{T}, nModel1::Int64, 
    Mhck1::AbstractVector{T}; vhth::AbstractVector{T}=[0.1, 1.0], vhth2::AbstractVector{T}=[0.1, 1.0]) where{T,N}

    fill!(J, 0.0)
    nh[1:nModel1] = x[1:2:end]
    vhth[1:nModel1] = x[2:2:end]

    vhth2[1:end-1] = vhth[1:end-1] .^ 2

    nh[end] = Mhck1[1] - sum(nh[1:end-1])

    # if nh[end] ≥ n0_c
    #     # Kh9 = (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1]))
    #     # if Kh9 ≥ epsT1000
    #     #     vhth2[end] = Kh9 / nh[end]
    #     #     vhth[end] = vhth2[end]^0.5
    #     # else
    #     #     vhth2[end] = vhthMin^2
    #     #     vhth[end] = vhthMin
    #     # end

        vhth2[end] = (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1])) / nh[end]
        vhth[end] = vhth2[end]^0.5
        
        nj = 1
        # (l, j) = (0, 4)   # j = 2(nj+1)
        s = 1
        ration = (1 + nh[nj] / nh[end])
        J[nj, s] = vhth2[nj]^2 - vhth2[end]^2 - 2vhth2[nj] * vhth2[end] * ration
        J[nj, s+1] = 4nh[nj] * vhth[nj] * (vhth2[nj] - vhth2[end])
    
        nj += 1
        # (l, j) = (0, 6)   # j = 2(nj+1)
        s = 1
        ration = (1 + nh[nj] / nh[end])
        J[nj, s] = vhth2[nj]^3 - vhth2[end]^3 - 3vhth2[nj] * vhth2[end]^2 * ration
        J[nj, s+1] = 6nh[nj] * vhth[nj] * (vhth2[nj]^2 - vhth2[end]^2)
        
        # (l, j) = (0, nj)
        for kM in 2:nModel1
            nj += 1
            ration = (1 + nh[nj] / nh[end])
            for s in 1:2
                s2 = 2(s - 1)
                # j = 2(nj + 1)
                j2 = (nj + 1)
                J[nj, s2+1] = vhth2[nj]^j2 - vhth2[end]^j2 - j2 * vhth2[nj] * vhth2[end]^(j2-1) * ration
                J[nj, s2+2] = 2j2 * nh[nj] * vhth[nj] * (vhth2[nj]^(j2-1) - vhth2[end]^(j2-1))
            end
        end
    # else
    #     @show vhth2
    #     @show (Mhck1[2] - sum(nh[1:end-1] .* vhth2[1:end-1]))
    #     nh[end] = 0.0
    #     vhth2[end] = vhthMin^2
    #     vhth[end] = vhthMin
        
    #     nj = 1
    #     # (l, j) = (0, 4)   # j = 2(nj+1)
    #     s = 1
    #     J[nj, s] = vhth2[nj]^2
    #     J[nj, s+1] = 4nh[nj] * vhth[nj] * vhth2[nj]
    
    #     nj += 1
    #     # (l, j) = (0, 6)   # j = 2(nj+1)
    #     s = 1
    #     J[nj, s] = vhth2[nj]^3
    #     J[nj, s+1] = 6nh[nj] * vhth[nj] * vhth2[nj]^2
        
    #     # (l, j) = (0, nj)
    #     for kM in 2:nModel1
    #         nj += 1
    #         for s in 1:2
    #             s2 = 2(s - 1)
    #             # j = 2(nj + 1)
    #             j2 = (nj + 1)
    #             J[nj, s2+1] = vhth2[nj]^j2
    #             J[nj, s2+2] = 2j2 * nh[nj] * vhth[nj] * vhth2[nj]^(j2-1)
    #         end
    #     end
    # end
    return nh[end], vhth[end]
end
