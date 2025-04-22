function [BFIM_tb,BFIM_xyb] = PFIM_MCMT_new(N_BS,K,Nt,Nr,L,A_prob,dA_prob,V_prob,dV_prob,B_prob,R,D,ind_C,C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob)
%   PFIM_MCMT_new Summary of this function goes here
%   Detailed explanation goes here
%   Receive Beamforming C is input


n_prob = size(tan_psens_prob,3);

for i = 1:n_prob
    [F_tb_prob,F_xyb_prob] = FIM_MCMT_new(N_BS,K,Nt,Nr,L,A_prob(:,:,:,i),dA_prob(:,:,:,i),V_prob(:,:,:,i),dV_prob(:,:,:,i),B_prob(:,:,:,:,i),R,D,ind_C,C,tan_psens_prob(:,:,i),cot_psens_prob(:,:,i),p_diff_prob(:,:,:,i));
        
    if i == 1
        BFIM_xyb = weight_prob(i)*F_xyb_prob;
        BFIM_tb  = weight_prob(i)*F_tb_prob;
    else
        BFIM_xyb = BFIM_xyb + weight_prob(i)*F_xyb_prob;
        BFIM_tb  = BFIM_tb  + weight_prob(i)*F_tb_prob;
    end

end





end

