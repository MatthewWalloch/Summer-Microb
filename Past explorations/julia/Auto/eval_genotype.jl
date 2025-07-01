# function for indvidial fitness evaluation 
function eval_genotype(fit_Pop::Vector{Float64},coopPayoff_Pop::Vector{Float64},coopCost_Pop::Vector{Float64},coopEff_Pop::Vector{Float64},sigCost_Pop::Vector{Float64},auto_pro_Rate::Vector{Float64},pro_Rate::Vector{Float64},sig_Th::Vector{Float64},auto_R::Vector{Float64},baseline::Float64,coop_Benefit::Float64,coop_Cost::Float64,sig_Cost::Float64,size_Pop::Int64,lambda::Float64,env_CellDen::Vector{Float64},grid_Size::Int64,base_Volume::Float64,decay_Rate::Float64,median_CellDen::Float64,K::Float64)

    # genotype evaluation
    counter = 0
    while counter<size_Pop
        # select number of mixing genotypes
        mix_Num = sample_ztp(lambda)

        # randomly select mix_Num genotypes
        index_Geno = rand(1:size_Pop,mix_Num)

        pool_CellDen = repeat(env_CellDen,outer=[1,mix_Num])


        # calculate local signal density for mix_Num genotypes
        sig_Concentration = cal_Signal_Concentration(repeat(pro_Rate[index_Geno]',outer=[grid_Size,1]),pool_CellDen,decay_Rate,repeat(auto_R[index_Geno]',outer=[grid_Size,1]),K)./mix_Num
        # individual who turns on cooperation at certain environment
        coop_ON = repeat(sig_Concentration,outer=[1,mix_Num]).>repeat((sig_Th[index_Geno])',outer=[grid_Size,1])

        # benefit for cooperation
        coopPayoff_Pop[index_Geno] .= coop_Benefit.*sum(sum(coop_ON)./mix_Num.*env_CellDen.*base_Volume.> fill(median_CellDen,grid_Size,1).*base_Volume)
        
        # efficiency for cooperation
        coopEff_Pop[index_Geno] .= sum(coop_ON.==repeat(env_CellDen.>median_CellDen,outer=[1,mix_Num]))/Float64(grid_Size)

        # cost for cooperation
        coopCost_Pop[index_Geno] .= coop_Cost.*sum(coop_ON)'

        # cost for signalling
        auto_pro_Rate[index_Geno] = pro_Rate[index_Geno].*(mean(sig_Concentration)./(K+mean(sig_Concentration)).*auto_R[index_Geno])
        sigCost_Pop[index_Geno] = sig_Cost.*(pro_Rate[index_Geno]+auto_pro_Rate[index_Geno])

        # calculate genotype fitness
        fit_Pop[index_Geno] .= baseline+coopPayoff_Pop[index_Geno][1]-coopCost_Pop[index_Geno][1]-sigCost_Pop[index_Geno][1]

        # update counter
        counter += mix_Num
    end

end
