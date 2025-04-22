clc;
clear;

% Note: scale_factor (like scale_R, scale_D) should be carefully selected,
% since cvx can usually report with error with inappropriate
% hyperparameters

dst_gau = @(ps,r) 1/(2*pi*(r/3)^2).*exp(-ps./(2*(r/3)^2));
dst_unf = @(ps,r) (ps./ps)*(1/r_sens);


N_BS = 2;
if N_BS == 2
    p_BS   = [sqrt(3)/2,0;-sqrt(3)/2,0];

end

N_US   = 2;
p_sens = [sqrt(3)/4,3/4;0,3/4];
K = N_US;

%%%%%% Study the impact of number of BS on PCRB
% N_BS = 7;
% % p_BS   = [sqrt(3)/2,-3/2;-sqrt(3)/2,-3/2];
% p_BS   = [sqrt(3)/2,-3/2;0,0];
% if N_BS >= 3
%     % p_BS   = cat(1,p_BS,[0,0]);
%     p_BS   = cat(1,p_BS,[-sqrt(3)/2,-3/2]);
% end
% if N_BS >= 4
%     p_BS   = cat(1,p_BS,[sqrt(3)/2,3/2]);
% end
% if N_BS >= 5
%     p_BS   = cat(1,p_BS,[-sqrt(3)/2,3/2]);
% end
% if N_BS >= 6
%     p_BS   = cat(1,p_BS,[sqrt(3),0]);
% end
% if N_BS >= 7
%     p_BS   = cat(1,p_BS,[-sqrt(3),0]);
% end
%%%%%% End: Study the impact of number of BS on PCRB



n_cent   = 1;
p_sens_p = zeros([K,2,n_cent]);
p_sens_p(:,:,1) = p_sens;
if n_cent > 1
    w_gau    = [0.8,0.2];
    p_sens_p(:,:,2) = [-sqrt(3)/3,3/5;0,3/4];
end

%%%%%% System setup
Nt    = 8;
Nr    = 8;
SNR   = 20;
var_n = 10^(-SNR/10);
TFn   = 8*ones([N_BS,1]);
% {ScalarQ: 1} or {VecterQ: 0}
is_diagQ = 0;   
%%%%%% End: System setup 

%%%%%%% user is random distributed
r0 = 0.1;
r_sens = [r0];
if K >= 2
    r_sens = cat(1,r_sens,[r0]); 
end
if K >= 3
    r_sens = cat(1,r_sens,[r0]); 
end
if K >= 4
    r_sens = cat(1,r_sens,[r0]); 
end
if K >= 5
    r_sens = cat(1,r_sens,[r0]); 
end


n_xy_prob      = 08;
n_prob         = (2*n_xy_prob+1)^2;
delta_p        = r_sens/2/n_xy_prob;
p_sens_prob    = zeros([K,2,n_cent*n_prob]);
pdf_prob       = zeros([K,n_prob]);
weight_k_prob  = zeros([K,n_prob]);

weight_k_prob_cent = zeros([K,n_cent*n_prob]);

for ik = 1:K
    delta_p_x_c  = [-n_xy_prob:1:-1,1:n_xy_prob]/n_xy_prob*r_sens(ik,1);
    delta_p_y_c  = [-n_xy_prob:1:-1,1:n_xy_prob]/n_xy_prob*r_sens(ik,1);
    delta_p_xy_c = zeros([2,4*n_xy_prob]);
    delta_p_xy_c(1,1:2*n_xy_prob)               = delta_p_x_c;
    delta_p_xy_c(2,(1:2*n_xy_prob)+2*n_xy_prob) = delta_p_y_c;

    delta_p_x_nc      = delta_p_x_c'*ones([1,2*n_xy_prob]);
    delta_p_y_nc      = ones([2*n_xy_prob,1])*delta_p_y_c;
    delta_p_xy_nc_mat = cat(3,delta_p_x_nc,delta_p_y_nc);
    delta_p_xy_nc     = reshape(delta_p_xy_nc_mat,[4*n_xy_prob*n_xy_prob,2]);
    delta_p_xy_nc     = delta_p_xy_nc.';

    delta_p_xy  = cat(2,delta_p_xy_c,delta_p_xy_nc);
    delta_p_xy0 = cat(2,delta_p_xy,[0;0]);
    if n_cent > 1
        delta_p_xy0 = cat(2,delta_p_xy0,delta_p_xy0);
    end

    delta_p_sdis        = sum(abs(delta_p_xy0).^2,1);
    if n_cent == 1
        weight_k_prob(ik,:) = dst_gau(delta_p_sdis,r_sens(ik,1));
        weight_k_prob_cent(ik,:) = dst_gau(delta_p_sdis,r_sens(ik,1));      % gaussian distributed users
    elseif n_cent > 1
        weight_cent_tmp = zeros([1,n_prob]);
        for icent = 1:n_cent
            icent_rg = ((icent-1)*n_prob+1):(icent*n_prob);
            weight_cent_tmp = w_gau(icent)*dst_gau(delta_p_sdis(:,icent_rg),r_sens(ik,1));      % uniform distributed users
            weight_k_prob_cent(ik,icent_rg) = weight_cent_tmp;
        end
    end
    
    if n_cent == 1
        p_xy     = (p_sens(ik,:).')*ones([1,n_prob]) + delta_p_xy0;
    elseif n_cent > 1
        p_sens_icent = [];
        for icent = 1:n_cent
            p_sens_icent = [p_sens_icent,(p_sens_p(ik,:,icent).')*ones([1,n_prob])];
        end
        p_xy = p_sens_icent + delta_p_xy0;
    end

    if n_cent == 1
        p_xy_ten = reshape(p_xy,[1,2,n_prob]);
    else
        p_xy_ten = reshape(p_xy,[1,2,n_cent*n_prob]);
    end
    p_sens_prob(ik,:,:) = p_xy_ten;

end

if n_cent == 1
    weight_prob = prod(weight_k_prob,1);
    weight_prob = weight_prob ./ (sum(weight_prob,2));
else
    weight_prob = prod(weight_k_prob_cent,1);
    weight_prob = weight_prob ./ (sum(weight_prob,2));
end

n_prob_ss = n_prob;
if n_cent > 1
    n_prob = n_prob*n_cent;
end

dis_psens_prob = zeros([K,N_BS,n_prob]);
tan_psens_prob = zeros([K,N_BS,n_prob]);
cot_psens_prob = zeros([K,N_BS,n_prob]);
p_diff_prob    = zeros([K,2,N_BS,n_prob]);
for ip = 1:n_prob
    for iN = 1:N_BS
        p_diff_prob(:,:,iN,ip)  = p_sens_prob(:,:,ip) - ones(K,1)*p_BS(iN,:);
        dis_psens_prob(:,iN,ip) = sqrt(sum(p_diff_prob(:,:,iN,ip).^2,2));
        tan_psens_prob(:,iN,ip) = p_diff_prob(:,1,iN,ip) ./ p_diff_prob(:,2,iN,ip);
        cot_psens_prob(:,iN,ip) = p_diff_prob(:,2,iN,ip) ./ p_diff_prob(:,1,iN,ip);
    end
end
theta_psens_prob = atan(tan_psens_prob); % (K , N_BS , n_prob)

rc_sens_prob = 1;
pl_sens_prob = 1./dis_psens_prob;    % （K , N_BS , n_prob)
rcs_sc_prob  = 0.1*ones([N_BS,N_BS,K,n_prob]);
prcs_prob    = zeros([N_BS,N_BS,K,n_prob]);
for ip = 1:n_prob
    for iN_rBS = 1:N_BS
        for iN_tBS = 1:N_BS
            for k = 1:K
                prcs_prob(iN_rBS,iN_tBS,k,ip) = pl_sens_prob(k,iN_rBS,ip)*pl_sens_prob(k,iN_tBS,ip)*rcs_sc_prob(iN_rBS,iN_tBS,k,ip);
            end
        end
    end
end

%%%%%%% End: user is random distributed

%%%%%% steering vector/matrix derivation
H_sens_prob    = zeros([Nr,Nt,N_BS,N_BS,n_prob]);
A_sens_prob    = zeros([Nr,K,N_BS,n_prob]);
V_sens_prob    = zeros([Nt,K,N_BS,n_prob]);
de_A_sens_prob = zeros([Nr,K,N_BS,n_prob]);
de_V_sens_prob = zeros([Nt,K,N_BS,n_prob]);
for ip = 1:n_prob
    for iN = 1:N_BS
        %%%%%% 1st antenna as reference point
        A_sens_prob(:,:,iN,ip)    = exp(-1j*pi*(0:(Nr-1))'*sin(theta_psens_prob(:,iN,ip)'));
        V_sens_prob(:,:,iN,ip)    = exp( 1j*pi*(0:(Nt-1))'*sin(theta_psens_prob(:,iN,ip)'));
        de_A_sens_prob(:,:,iN,ip) = (-1j*pi*(0:(Nr-1))'*cos(theta_psens_prob(:,iN,ip)')) .* A_sens_prob(:,:,iN,ip);
        de_V_sens_prob(:,:,iN,ip) = ( 1j*pi*(0:(Nt-1))'*cos(theta_psens_prob(:,iN,ip)')) .* V_sens_prob(:,:,iN,ip);
    end
    for iN = 1:N_BS
        for iN2 = 1:N_BS
            H_sens_prob(:,:,iN,iN2,ip) = A_sens_prob(:,:,iN,ip)*diag(reshape(prcs_prob(iN,iN2,:,ip),[K,1]))*V_sens_prob(:,:,iN2,ip).';
        end
    end
end

B_block_prob = zeros([K,N_BS*K,N_BS,n_prob]);
Bnu_prob     = zeros([K,K,N_BS,N_BS,n_prob]);
for ip = 1:n_prob
    for in = 1:N_BS
        Bn_block_prob = [];
        for iu = 1:N_BS
            Bnu_prob(:,:,in,iu,ip) = diag(reshape(prcs_prob(in,iu,:,ip),[K,1]));
            Bn_block_prob          = [Bn_block_prob,Bnu_prob(:,:,in,iu,ip)];
        end 
        B_block_prob(:,:,in,ip) = Bn_block_prob;
    end
end
VV_diagT_prob = zeros([N_BS*K,N_BS*Nt,n_prob]);
VV_diagC_prob = zeros([N_BS*Nt,N_BS*K,n_prob]);
for ip = 1:n_prob
    VV = [];
    for i = 1:N_BS
        VV = [VV;V_sens_prob(:,:,i,ip)];
    end
    VV_diagT_prob(:,:,ip) = kron(ones(N_BS,1),VV.')    .*kron(eye(N_BS),ones([K,Nt]));
    VV_diagC_prob(:,:,ip) = kron(ones(1,N_BS),conj(VV)).*kron(eye(N_BS),ones([Nt,K]));
end
%%%%%% End: steering vector/matrix derivation

%%%%%% projection matrix derivation
%%%%%%%%%%%%% Scheme Selection %%%%%%%%%%%%%%%%%%
ind_C = 0; % 0: no projection, 1: lossless projection, 2: DFT projection, 3: K dimensional DFT(theta_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ind_C == 1
    L = 2*K; %L=K+1;
    C = zeros([Nr,L,N_BS]);
    for iN = 1:N_BS
        PidA_n = (eye(Nr)-A_sens_prob(:,:,iN,n_prob)*inv(A_sens_prob(:,:,iN,n_prob)'*A_sens_prob(:,:,iN,n_prob))*A_sens_prob(:,:,iN,n_prob)')*de_A_sens_prob(:,:,iN,n_prob);
        C_tmp = [A_sens_prob(:,:,iN,n_prob),PidA_n];
        [UCCH,eCCH,VCCH] = svd(C_tmp*C_tmp');
        C(:,:,iN) = UCCH(:,1:L);

    end

elseif ind_C == 2
    L = 2*K; 
    C = zeros([Nr,L,N_BS]);
    
    tmpaa1 = zeros([L,Nr,N_BS]);
    tmpaa1_max = zeros([N_BS,L]);
    for iN = 1:N_BS
        DFT_mtx = 1/sqrt(Nr)*fft(eye(Nr));
        % Cc_DFT  = randperm(Nr);
        PidA_n = (eye(Nr)-A_sens_prob(:,:,iN,n_prob)*inv(A_sens_prob(:,:,iN,n_prob)'*A_sens_prob(:,:,iN,n_prob))*A_sens_prob(:,:,iN,n_prob)')*de_A_sens_prob(:,:,iN,n_prob);
        C_tmp  = [A_sens_prob(:,:,iN,n_prob),PidA_n];
        C_tmp1 = C_tmp ./ (ones(Nr,1)*sqrt(sum(abs(C_tmp).^2,1)));

        % tmpaa1(:,:,iN) = abs((A_sens_prob(:,:,iN,ip)./(ones(Nr,1)*sqrt(sum(abs(A_sens_prob(:,:,iN,ip)).^2,1))))'*DFT_mtx);
        tmpaa1(:,:,iN) = abs(C_tmp1'*DFT_mtx);
        for im = 1:L
            [~,tmpaa1_max(iN,im)] = max(tmpaa1(im,:,iN));
        end

        if iN == 1
            % C(:,:,iN) = DFT_mtx(:,[5,6,7,8]);
            % C(:,:,iN) = DFT_mtx(:,[7,8,9,10]);    % Nr=10;
            % C(:,:,iN) = DFT_mtx(:,[12,13,14,16]);   % Nr=20;
            % C(:,:,iN) = DFT_mtx(:,[24,20,23,19]);   % Nr=30;
            C(:,:,iN) = DFT_mtx(:,[31,26,32,25]);   % Nr=40;
            C(:,:,iN) = DFT_mtx(:,[39,32,33,38]);   % Nr=50;
            C(:,:,iN) = DFT_mtx(:,[46,38,45,39]);   % Nr=60;
            C(:,:,iN) = DFT_mtx(:,[54,45,53,44]);   % Nr=70;
        elseif iN == 2
            % C(:,:,iN) = DFT_mtx(:,[3,4,5,6]);
            % C(:,:,iN) = DFT_mtx(:,[4,5,6,7]);    % Nr=10;
            % C(:,:,iN) = DFT_mtx(:,[8,9,10,11]);   % Nr=20;
            % C(:,:,iN) = DFT_mtx(:,[14,12,13,11]);   % Nr=30;
            C(:,:,iN) = DFT_mtx(:,[16,17,18,19]);   % Nr=40;
            C(:,:,iN) = DFT_mtx(:,[23,20,22,19]);   % Nr=50;
            C(:,:,iN) = DFT_mtx(:,[27,24,26,23]);   % Nr=60;
            C(:,:,iN) = DFT_mtx(:,[31,27,32,28]);   % Nr=70;
        end
        test = C(:,:,iN)'*C(:,:,iN);
    end
elseif ind_C == 3
    % L = 2*K; 
    L = 2*K-1;
    % L = 2*K+1;
    C = zeros(Nr,L,N_BS);
    for iN = 1:N_BS
        DFT_mtx = 1/sqrt(Nr)*fft(eye(Nr));
        Cc_DFT  = randperm(Nr);
        if iN == 1
            C(:,:,iN) = DFT_mtx(:,[8,6,7]);
        elseif iN == 2
            C(:,:,iN) = DFT_mtx(:,[4,5,2]);
        end
        test = C(:,:,iN)'*C(:,:,iN);
    end

elseif ind_C == 0
    L = Nr;
    C = zeros(Nr,L,N_BS);
    for iN = 1:N_BS
        DFT_mtx = 1/sqrt(Nr)*fft(eye(Nr));
        C(:,:,iN) = DFT_mtx;
        C(:,:,iN) = eye(Nr);
    end
elseif ind_C == 4
    L = Nr;
    C = zeros([Nr,L,N_BS]);

    H_tmp  = zeros([Nt*N_BS,Nr,N_BS]);
    HH_tmp = zeros([Nr,Nr,N_BS]);
    for iN = 1:N_BS
        H_tmp(:,:,iN)  = VV_diagC_prob(:,:,n_prob)*B_block_prob(:,:,i,n_prob)'*A_sens_prob(:,:,i,n_prob)';
        HH_tmp(:,:,iN) = H_tmp(:,:,iN)'*H_tmp(:,:,iN);
        [UHH,eHH,VHH]  = svd(HH_tmp(:,:,iN));
        C(:,:,iN)      = UHH;
    end
end
%%%%%% End: projection matrix derivation

tic;
%%%%%% Initialization of R
%%%%%% optimize transmit covariance: R with no quantization
I_F    = eye(K*N_BS+2*K*N_BS*N_BS);
I_Fxyb = eye(2*K   +2*K*N_BS*N_BS);
I_R    = eye(Nt*N_BS);
I_BS   = eye(N_BS);
I_subR = zeros([Nt*N_BS,Nt,N_BS]);
for i = 1:N_BS
    I_subR(:,:,i) = kron(I_BS(:,i),eye(Nt));
end
I_Fxybp1 = eye(2*K+2*K*N_BS*N_BS+1);


% test prior
p0 = 1/(r0^2);
Fpr = zeros([2*K+2*K*N_BS*N_BS,2*K+2*K*N_BS*N_BS]);
Fpr(1:2*K,1:2*K) = p0*eye(2*K); %Fpr(3:4,3:4) = 100*eye(2);
Fpr((2*K+1):(2*K+2*K*N_BS*N_BS),(2*K+1):(2*K+2*K*N_BS*N_BS)) = 100000000*eye(2*K*N_BS*N_BS);

eps_R = 100;
rho1  = 1e-6;
scale_R_ini = 1e-6;

IND_MAT = zeros([2*K+2*K*N_BS*N_BS,2*K]);
IND_MAT(1:2*K,1:2*K) = eye(2*K);

% When study PCRB vs N_BS, 1.(N_BS>=4,scale_R_ini=1e-4); 2.(N_BS=3,scale_R_ini=1e-3);
%%%%%%%% cvx start
cvx_precision best
cvx_begin sdp
    variable Rini(Nt,Nt,N_BS) hermitian semidefinite complex
    variable t_ini(2*K)

    % R_ini_opt = I_subR(:,:,1)*Rini(:,:,1)*I_subR(:,:,1)' + I_subR(:,:,2)*Rini(:,:,2)*I_subR(:,:,2)';
    for i = 1:N_BS
        if i == 1
            R_ini_opt = I_subR(:,:,i)*Rini(:,:,i)*I_subR(:,:,i)';
        else
            R_ini_opt = R_ini_opt + I_subR(:,:,i)*Rini(:,:,i)*I_subR(:,:,i)';
        end
    end

    D_ini  = eye(L*N_BS);
    if is_diagQ == 0
        [~,F1] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                           R_ini_opt,D_ini,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
    elseif is_diagQ == 1
        [~,F1] = PFIM_MCMT_new(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                           R_ini_opt,D_ini,ind_C,C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
    end

    minimize(sum(t_ini))
    subject to
        for i2K = 1:(2*K)
            [scale_R_ini*(F1+Fpr),I_Fxyb(:,i2K);I_Fxyb(i2K,:),t_ini(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
        end

        for iNB = 1:N_BS
            real(trace(Rini(:,:,iNB)))*1000 <= 8/var_n*1000;       % change Nt
        end
cvx_end
R_ini = R_ini_opt;
%%%%%% End: Initialization of R (Optimization of R without quantization)
toc;


%%%%% Derivation of Jn      
if TFn(1) > 1
    R0 = 1/var_n*eye(Nt*N_BS);
else
    R0 = R_ini;
end

%%%%% End: Derivation of Jn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
J           = zeros([L,L,N_BS]);
J_org_prob  = zeros([Nr,Nr,N_BS,n_prob]);
G_prob      = zeros([Nt*N_BS,L,N_BS,n_prob]);
G_org_prob  = zeros([Nt*N_BS,Nr,N_BS,n_prob]);
J_test      = zeros([L,L,N_BS]);
Jn_prob     = zeros([L,L,N_BS,n_prob]);
I_RBS       = zeros([N_BS*Nt,Nt,N_BS]);
for i = 1:N_BS
    I_RBS(:,:,i) = kron(I_BS(:,i),eye(Nt));
end
for i = 1:N_BS
    for ip = 1:n_prob
        G_org_prob(:,:,i,ip) = VV_diagC_prob(:,:,ip)*B_block_prob(:,:,i,ip)'*A_sens_prob(:,:,i,ip)';
        J_org_prob(:,:,i)    = G_org_prob(:,:,i,ip)'*R0*G_org_prob(:,:,i,ip);
        G_prob(:,:,i,ip)     = G_org_prob(:,:,i,ip)*C(:,:,i);
        Jn_prob(:,:,i,ip)    = G_prob(:,:,i,ip)'*R0*G_prob(:,:,i,ip);
        if ip == 1
            Jns = weight_prob(ip)*Jn_prob(:,:,i,ip);
        else
            Jns = Jns + weight_prob(ip)*Jn_prob(:,:,i,ip);
        end
    end
    J(:,:,i)     = Jns;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Derivation of Tn (inverse of the noise covariance)
C_mat  = reshape(C,[Nr,L*N_BS]);
T      = (C_mat'*C_mat).*kron(eye(N_BS),ones([L,L]));
I_subL = zeros(N_BS*L,L,N_BS);
T_ten  = zeros(L,L,N_BS);
for i = 1:N_BS
    I_subL(:,:,i) = kron(I_BS(:,i),eye(L));
end
D0    = zeros(N_BS*L,N_BS*L); % inverse of Tn: D0 = Tn^{-1} = (sigma_n * C^H * C)^{-1}
for i = 1:N_BS
    rg           = ((i-1)*L+1):(i*L);
    Tn           = I_subL(:,:,i).'*T*I_subL(:,:,i);
    D0(rg,rg)    = eye(L)/Tn;
    T_ten(:,:,i) = Tn;
    eig_test     = eig(Tn);
end
%%%%%%% End: Derivation of Tn

if is_diagQ == 1 || SNR >= 0        
%%%%%% Derivation of the intialization of On (Dn): Identical Uniform On
    D_ini      = zeros([L,L,N_BS]);
    JpT        = J + T_ten;
    Qvar_n     = zeros([L,L,N_BS]);
    for i = 1:N_BS
        tmp_term_On   = 2^(TFn(i)/L);
        tmp_On        = tmp_term_On-1;
        JpTn          = JpT(:,:,i);
        if is_diagQ == 1
            JpTn = eye(L).*JpT(:,:,i);
        end
        [eigv_n,eig_vn,~] = svd(JpTn);
        eig_n         = real(diag(eig_vn));
        eig_On        = eig_n./tmp_On;
        %%%%%% Bisection search for On (Uniform)
        eig_max       = max(eig_On);
        eig_min       = min(eig_On);
        eig_sub       = eig_max - eig_min;
        while(eig_sub >= 1e-2)
            eig_ave   = (eig_max+eig_min)/2;
            Qvar_ave  = eigv_n*diag(eig_ave)*eigv_n';
            TFn_ave   = log2(real(det(JpTn+Qvar_ave))) - log2(real(det(Qvar_ave)));
            if TFn_ave > TFn(i)
                eig_min = eig_ave;
            else
                eig_max = eig_ave;
            end
            eig_sub = eig_max - eig_min;
        end
        eig_fin  = 1*eig_max;
        %%%%%% Tn is scalared identity matrix, initialize On to be also
        %%%%%% scalared identity matrix (if)
        Qvar_n(:,:,i) = eig_fin*eye(L);
        D_ini(:,:,i)  = 1/(1+eig_fin)*eye(L);
    end
    D_ini_diag = kron(ones([N_BS,1]),reshape(D_ini,[L,L*N_BS]));
    D_ini_diag = D_ini_diag.*kron(eye(N_BS),ones([L,L]));
%%%%%%% End: Derivation of the intialization of On (Dn): Uniform On
else
%%%%%% Derivation of the intialization of On (Dn): Non-identical On
    D_ini      = zeros([L,L,N_BS]);
    JpT        = J + T_ten;
    Qvar_n     = zeros([L,L,N_BS]);
    for i = 1:N_BS
        tmp_term_On   = 2^(TFn(i)/L);
        tmp_On        = tmp_term_On-1;
        JpTn          = JpT(:,:,i);
        [eigv_n,eig_vn,~] = svd(JpTn);
        eig_n         = real(diag(eig_vn));
        eig_On        = eig_n./tmp_On;
        %%%%%% Tn is scalared identity matrix, initialize On to be also
        %%%%%% scalared identity matrix (if)
        Qvar_n(:,:,i) = eigv_n*diag(5*eig_On)*eigv_n';
        D_ini(:,:,i)  = eye(L)/(T_ten(:,:,i)+Qvar_n(:,:,i));
    end
    D_ini_diag = kron(ones([N_BS,1]),reshape(D_ini,[L,L*N_BS]));
    D_ini_diag = D_ini_diag.*kron(eye(N_BS),ones([L,L]));
%%%%%% End: Derivation of the intialization of On (Dn): Non-identical On
end



I_F    = eye(K*N_BS+2*K*N_BS*N_BS);
I_Fxyb = eye(2*K+2*K*N_BS*N_BS);
I_R    = eye(Nt*N_BS);
I_subR = zeros([Nt*N_BS,Nt,N_BS]);
for i = 1:N_BS
    I_subR(:,:,i) = kron(I_BS(:,i),eye(Nt));
end

iter   = 20;
O      = zeros([L,L,N_BS]);
Mat_ZL = zeros([L,L]);

t_record1 = [];
it_RQ     = [];
t_record  = [];

%%%%%% Optimize either O (1) or R (2)
is_OptEOR = 0;
if is_OptEOR == 1

    D_update = D_ini_diag;
    R_update = R0;
    R_new    = R_update;
    eps      = 100;

    while(eps > 1e-3)
    
        D_past     = D_update;
        D_current  = D_past;
        tem_SCA1   = zeros([N_BS,1]);
        tem_SCA2   = zeros([L,L,N_BS]);
        Dn_current = zeros([L,L,N_BS]);
        Sigma      = zeros([L,L,N_BS]);
        for i = 1:N_BS
            Dn_current(:,:,i) = I_subL(:,:,i).'*D_current*I_subL(:,:,i);
            Sigma(:,:,i)      = eye(L) + J(:,:,i)^(0.5)*Dn_current(:,:,i)*(J(:,:,i)^(0.5));
            tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
            tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
        end

        scale_D = 2e-4; % Nr=8, 2e-5 Nr=32, 1e-2

        T_i = zeros([L,L,N_BS]);

        cvx_precision best
        cvx_begin SDP
            variable D_opt(L,L,N_BS) hermitian semidefinite complex
            variable tQ(2*K)
    
            for i = 1:N_BS
                if i == 1
                    D_mat = I_subL(:,:,i)*(D_opt(:,:,i))*I_subL(:,:,i)';
                else
                    D_mat = D_mat + I_subL(:,:,i)*(D_opt(:,:,i))*I_subL(:,:,i)';
                end
                T_i(:,:,i) = I_subL(:,:,i).'*D0*I_subL(:,:,i);
            end
    
    
            [~,F2] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                   R_new,D_mat,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
                

            minimize(sum(tQ))
            subject to
                for i2K = 1:(2*K)
                    [scale_D*(F2+Fpr),I_Fxyb(:,i2K);I_Fxyb(i2K,:),tQ(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                end
    
                %%%%%% Noise variance is positive semidefinite
                for i = 1:N_BS
                    (T_i(:,:,i) - D_opt(:,:,i))*10000 >= 0;
                end

                %%%%%% Lemma 1 approximation
                for i = 1:N_BS
                    (tem_SCA1(i) + real(trace(Sigma(:,:,i)\(eye(L)+J(:,:,i)^(0.5)*D_opt(:,:,i)*(J(:,:,i)^(0.5)))))/log(2) ...
                    - log_det(eye(L)-T_ten(:,:,i)^(0.5)*D_opt(:,:,i)*T_ten(:,:,i)^(0.5))/log(2) - L/log(2))*0.003 <= TFn(i)*0.003;
                end
        cvx_end

        for i = 1:N_BS
            Dn_current(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
            Sigma(:,:,i)      = eye(L) + J(:,:,i)^(0.5)*Dn_current(:,:,i)*(J(:,:,i)^(0.5));
            tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
            tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
        end
        disp(tmpu);

    end
    D_new = D_update;
    
elseif is_OptEOR == 2
       
    it_OptEOR = 16;
    
    R_EOR    = R0;
    D_update = D_ini_diag;

    for itR = 1:it_OptEOR

        % D_update = D_ini_diag;
        R_update = R_EOR;
    
        [~,PFu]  = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                     R_update,D_update,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
        PCRBu    = inv(PFu);
        mse_record0 = trace(PCRBu(1:(2*K),1:(2*K)));
        disp(mse_record0);
    
        Dn_update = zeros(L,L,N_BS);
        for i = 1:N_BS
            Dn_update(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
        end
        eps_R = 100;
        rho1  = 6e-3;
        scale_R = 1e-3; %1e-3 %Nr=32,ind_C=2,choice=1,2e-3 %reference central: 8e-4
    
        %%% opt Rn

        while(eps_R > 1e-2)

            if itR == 1
                break;
            end

            R_past       = R_update;
            SigmaR       = zeros([L,L,N_BS]);
            GRG_prob     = zeros([L,L,N_BS,n_prob]);
            tmp_MIs      = zeros([N_BS,1]);
            On_update    = zeros([L,L,N_BS]);
            invDn_update = zeros([L,L,N_BS]);
            for i = 1:N_BS
                On_update(:,:,i)    = eye(L)/(Dn_update(:,:,i))-T_ten(:,:,i);
                invDn_update(:,:,i) = eye(L)/(Dn_update(:,:,i));
                    
                for ip = 1:n_prob
                    GRG_prob(:,:,i,ip) = G_prob(:,:,i,ip)'*R_past*G_prob(:,:,i,ip);
                    SigmaR(:,:,i)      = SigmaR(:,:,i) + weight_prob(ip)*GRG_prob(:,:,i,ip);
                end
                SigmaR(:,:,i) = SigmaR(:,:,i) + invDn_update(:,:,i);
                tmp1_MI       = log2(real(det(SigmaR(:,:,i))));
                [~,On_svd,~]  = svd(On_update(:,:,i));
            end
        
            cvx_precision best
            cvx_begin sdp
                variable Rnn(Nt,Nt,N_BS) hermitian semidefinite complex
                variable t(2*K)
                    
                for i = 1:N_BS
                    if i == 1
                        R1 = I_subR(:,:,i)*Rnn(:,:,i)*I_subR(:,:,i)';
                    else
                        R1 = R1 + I_subR(:,:,i)*Rnn(:,:,i)*I_subR(:,:,i)';
                    end
                end
    
                J_optVar = [];
                for i = 1:N_BS
                    for ip = 1:n_prob
                        if ip == 1
                            J_optVar_iBS = weight_prob(ip)*G_prob(:,:,i,ip)'*R1*G_prob(:,:,i,ip);
                        else
                            J_optVar_iBS = J_optVar_iBS + weight_prob(ip)*G_prob(:,:,i,ip)'*R1*G_prob(:,:,i,ip);
                        end
                    end
                    J_optVar = cat(3,J_optVar,J_optVar_iBS);
                end
    
                [~,FR] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                       R1,D_update,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
    
                minimize(sum(t))
                subject to
                    for i2K = 1:(2*K)
                        [scale_R*(FR+Fpr),I_Fxyb(:,i2K);I_Fxyb(i2K,:),t(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                    end
                    for iNB = 1:N_BS
                        real(trace(Rnn(:,:,iNB)))*100 <= 8/var_n*100;      % Nt changes
                        (real(tmp_MIs(iNB)) + real(trace(SigmaR(:,:,iNB)\(J_optVar(:,:,iNB)+inv(Dn_update(:,:,iNB)))))/log(2) - L/log(2))*0.01 <= TFn(iNB)*0.01; %0.001
                    end
            cvx_end
    
            R_update = R1;
            eps_R    = norm(R_update-R_past,'fro')/(N_BS*Nt*N_BS*Nt)/norm(R_past,'fro');
        
            tmp_cons_fr = zeros([N_BS,1]);
            for i = 1:N_BS
                tmp_cons_fr(i) = (real(tmp_MIs(iNB)) + real(trace(SigmaR(:,:,iNB)\(J_optVar(:,:,iNB)+inv(Dn_update(:,:,iNB)))))/log(2) - L/log(2));
            end
            disp(tmp_cons_fr);
        
            GRG_prob1 = zeros([L,L,N_BS,n_prob]);
            SigmaR1   = zeros([L,L,N_BS]);
            J_update  = zeros([L,L,N_BS]);
            for i = 1:N_BS
                On_update(:,:,i)    = eye(L)/(Dn_update(:,:,i))-T_ten(:,:,i);
                invDn_update(:,:,i) = eye(L)/(Dn_update(:,:,i));
                    
                for ip = 1:n_prob
                    GRG_prob1(:,:,i,ip) = G_prob(:,:,i,ip)'*R_update*G_prob(:,:,i,ip);
                    SigmaR1(:,:,i)      = SigmaR1(:,:,i) + weight_prob(ip)*GRG_prob1(:,:,i,ip);
                end
                J_update(:,:,i) = SigmaR1(:,:,i);
                SigmaR1(:,:,i)  = SigmaR1(:,:,i) + invDn_update(:,:,i);
                tmp1_MI1        = log2(real(det(SigmaR1(:,:,i))));
                [~,On_svd1,~]   = svd(On_update(:,:,i));
                tmp2_MI1        = log2(real(det(On_update(:,:,i))));
        
            end
                
        end
        R_new = R_update;
        R_EOR = R_update;
                    
        if itR >= 2
        for i = 1:N_BS
            J(:,:,i) = J_update(:,:,i);
        end 
        end
    
        eps   = 100;
        
        if itR >= 1
        %%% opt On
        while(eps > 1e-2)

            D_past     = D_update;
            D_current  = D_past;
            tem_SCA1   = zeros([N_BS,1]);
            tem_SCA2   = zeros([L,L,N_BS]);
            Dn_current = zeros([L,L,N_BS]);
            Sigma      = zeros([L,L,N_BS]);

            for i = 1:N_BS
                Dn_current(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
                Sigma(:,:,i)      = eye(L) + J(:,:,i)^(0.5)*Dn_current(:,:,i)*(J(:,:,i)^(0.5));
                tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
                tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);

            end
            disp(tmpu);

            


            scale_D = 1e-3;
            cvx_begin

                variable d_opt(L*N_BS,1)
                variable tQv(2*K)

                D_mat = eye(L*N_BS).*(d_opt*ones([1,L*N_BS]));
               
                [~,F2] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                   R_new,D_mat,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);


                D_opt = [];
                for i = 1:N_BS
                    if i == 1
                        D_opt = I_subL(:,:,i).'*D_mat*I_subL(:,:,i);
                    else
                        D_opt = cat(3,D_opt,I_subL(:,:,i)'*D_mat*I_subL(:,:,i));
                    end
                end

                T_i = diag(D0);

                minimize(sum(tQv))
                subject to
                    for i2K = 1:(2*K)
                        [scale_D*(F2+Fpr),I_Fxyb(:,i2K);I_Fxyb(i2K,:),tQv(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                    end

                    real(T_i - d_opt)*1000 >= 0;

                    %%%%%% Lemma 1 approximation
                    for i = 1:N_BS
                        (tem_SCA1(i) + real(trace(Sigma(:,:,i)\(eye(L)+J(:,:,i)^(0.5)*D_opt(:,:,i)*(J(:,:,i)^(0.5)))))/log(2) ...
                        - log_det(eye(L)-T_ten(:,:,i)^(0.5)*D_opt(:,:,i)*T_ten(:,:,i)^(0.5))/log(2) - L/log(2))*0.003 <= TFn(i)*0.003;
                    end

                    real(d_opt) >= 0;

            cvx_end

            CRBF2    = inv(F2+Fpr);
            disp(trace(CRBF2(1:(2*K),1:(2*K))));

            D_update = D_mat.*kron(eye(N_BS),ones([L,L]));
            eps      = norm(D_update-D_past,'fro')/(N_BS*L)/norm(D_past,'fro');

            for i = 1:N_BS
                Dn_current(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
                Sigma(:,:,i)      = eye(L) + J(:,:,i)^(0.5)*Dn_current(:,:,i)*(J(:,:,i)^(0.5));
                tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
                tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
            end
            disp(tmpu);

        end
        D_new = D_update;
        end

        disp(itR);
    end

end

tic;
%%%%%% Optimization (Generalized N BSs scenario)
if is_diagQ == 0        %%%% Vector Quantization
    for it = 1:iter
        if it == 1
            D_update  = D_ini_diag;
            R_update  = R0;
        else
            D_update  = D_new;
            R_update  = R_new;
        end
        
        
        %%%%%% optimize quantization noise covariance: O (here optimize D = inv(T+O) )
        R_new = zeros([N_BS*Nt,N_BS*Nt]) + R_update;
        eps   = 100;
    
        while(eps > 1e-2)
            if it == 1
                if K <= 3
                    break;    % K <= 3;
                end
            end
    
            D_past     = D_update;
            D_current  = D_past;
            tem_SCA1   = zeros([N_BS,1]);
            tem_SCA2   = zeros([L,L,N_BS]);
            Dn_current = zeros([L,L,N_BS]);
            Sigma      = zeros([L,L,N_BS]);
            for i = 1:N_BS
                Dn_current(:,:,i) = I_subL(:,:,i).'*D_current*I_subL(:,:,i);
                Sigma(:,:,i)      = eye(L) + J(:,:,i)^(0.5)*Dn_current(:,:,i)*(J(:,:,i)^(0.5));
                tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
                tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
                
            end

            scale_D = 2e-4; % Nr=8, 2e-5 Nr=32, 1e-2

            T_i = zeros([L,L,N_BS]);

            cvx_precision best
            cvx_begin SDP
                variable D_opt(L,L,N_BS) hermitian semidefinite complex
                variable tQ(2*K)
    
                for i = 1:N_BS
                    if i == 1
                        D_mat = I_subL(:,:,i)*(D_opt(:,:,i))*I_subL(:,:,i)';
                    else
                        D_mat = D_mat + I_subL(:,:,i)*(D_opt(:,:,i))*I_subL(:,:,i)';
                    end
                    T_i(:,:,i) = I_subL(:,:,i).'*D0*I_subL(:,:,i);
                end
    
    
                [~,F2] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                   R_new,D_mat,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
                

                minimize(sum(tQ))
                subject to
                    for i2K = 1:(2*K)
                        [scale_D*(F2+Fpr),I_Fxyb(:,i2K);I_Fxyb(i2K,:),tQ(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                    end
    
                    %%%%%% Noise variance is positive semidefinite
                    for i = 1:N_BS
                        (T_i(:,:,i) - D_opt(:,:,i))*10000 >= 0;
                    end

                    %%%%%% Lemma 1 approximation
                    for i = 1:N_BS
                        (tem_SCA1(i) + real(trace(Sigma(:,:,i)\(eye(L)+J(:,:,i)^(0.5)*D_opt(:,:,i)*(J(:,:,i)^(0.5)))))/log(2) ...
                        - log_det(eye(L)-T_ten(:,:,i)^(0.5)*D_opt(:,:,i)*T_ten(:,:,i)^(0.5))/log(2) - L/log(2))*0.003 <= TFn(i)*0.003;
                    end
            cvx_end

            D_update = D_mat.*kron(eye(N_BS),ones([L,L]));
            eps      = norm(D_update-D_past,'fro')/(N_BS*L)/norm(D_past,'fro');
           
            for i = 1:N_BS
                Dn_current(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
                Sigma(:,:,i)      = eye(L) + J(:,:,i)^(0.5)*Dn_current(:,:,i)*(J(:,:,i)^(0.5));
                tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
                tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
    
            end
            disp(tmpu);

        end
        D_new = D_update;
        %%%%%% End: optimize quantization noise covariance: O ( here optimize D = inv(T+O) )
        toc;
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%% optimize transmit covariance: R
        if it == iter
        % if it ~= 0
            break;
        end
        R_new = R_update;
        Dn_update = zeros(L,L,N_BS);
        for i = 1:N_BS
            Dn_update(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
        end
        eps_R = 100;
        rho1  = 6e-3;
        scale_R = 1e-3; %1e-3 %Nr=32,ind_C=2,choice=1,2e-3 %reference central: 8e-4


        itR = 1;
        while(eps_R > 1e-1)% || itR <= 2)              % Note: 1.                 (most cases, esp_R>1e-1)
            R_past       = R_update;
            SigmaR       = zeros([L,L,N_BS]);
            GRG_prob     = zeros([L,L,N_BS,n_prob]);
            tmp_MIs      = zeros([N_BS,1]);
            On_update    = zeros([L,L,N_BS]);
            invDn_update = zeros([L,L,N_BS]);
            for i = 1:N_BS
                On_update(:,:,i)    = eye(L)/(Dn_update(:,:,i))-T_ten(:,:,i);
                invDn_update(:,:,i) = eye(L)/(Dn_update(:,:,i));
                
                for ip = 1:n_prob
                    GRG_prob(:,:,i,ip) = G_prob(:,:,i,ip)'*R_past*G_prob(:,:,i,ip);
                    SigmaR(:,:,i)      = SigmaR(:,:,i) + weight_prob(ip)*GRG_prob(:,:,i,ip);
                end
                SigmaR(:,:,i) = SigmaR(:,:,i) + invDn_update(:,:,i);
                [~,On_svd,~]  = svd(On_update(:,:,i));

            end
    
            cvx_precision best
            cvx_begin sdp
                variable Rnn(Nt,Nt,N_BS) hermitian semidefinite complex
                variable t(2*K)
                
                for i = 1:N_BS
                    if i == 1
                        R1 = I_subR(:,:,i)*Rnn(:,:,i)*I_subR(:,:,i)';
                    else
                        R1 = R1 + I_subR(:,:,i)*Rnn(:,:,i)*I_subR(:,:,i)';
                    end
                end

                J_optVar = [];
                for i = 1:N_BS
                    for ip = 1:n_prob
                        if ip == 1
                            J_optVar_iBS = weight_prob(ip)*G_prob(:,:,i,ip)'*R1*G_prob(:,:,i,ip);
                        else
                            J_optVar_iBS = J_optVar_iBS + weight_prob(ip)*G_prob(:,:,i,ip)'*R1*G_prob(:,:,i,ip);
                        end
                    end
                    J_optVar = cat(3,J_optVar,J_optVar_iBS);
                end

                [~,FR] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                   R1,D_update,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);

                minimize(sum(t))
                subject to
                    for i2K = 1:(2*K)
                        [scale_R*(FR+Fpr),I_Fxyb(:,i2K);I_Fxyb(i2K,:),t(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                    end
                    for iNB = 1:N_BS
                        real(trace(Rnn(:,:,iNB)))*100 <= 8/var_n*100;      % Nt changes
                        
                        (real(tmp_MIs(iNB)) + real(trace(SigmaR(:,:,iNB)\(J_optVar(:,:,iNB)+inv(Dn_update(:,:,iNB)))))/log(2) - L/log(2))*0.01 <= TFn(iNB)*0.01; %0.001
                    end
            cvx_end

            
            R_update = R1;
            eps_R    = norm(R_update-R_past,'fro')/(N_BS*Nt*N_BS*Nt)/norm(R_past,'fro');
            
            tmp_cons_fr = zeros([N_BS,1]);
            for i = 1:N_BS
                tmp_cons_fr(i) = (real(tmp_MIs(iNB)) + real(trace(SigmaR(:,:,iNB)\(J_optVar(:,:,iNB)+inv(Dn_update(:,:,iNB)))))/log(2) - L/log(2));
            end
            disp(tmp_cons_fr);
    
            GRG_prob1 = zeros([L,L,N_BS,n_prob]);
            SigmaR1   = zeros([L,L,N_BS]);
            J_update  = zeros([L,L,N_BS]);
            for i = 1:N_BS
                On_update(:,:,i)    = eye(L)/(Dn_update(:,:,i))-T_ten(:,:,i);
                invDn_update(:,:,i) = eye(L)/(Dn_update(:,:,i));
                
                for ip = 1:n_prob
                    GRG_prob1(:,:,i,ip) = G_prob(:,:,i,ip)'*R_update*G_prob(:,:,i,ip);
                    SigmaR1(:,:,i)      = SigmaR1(:,:,i) + weight_prob(ip)*GRG_prob1(:,:,i,ip);
                end
                J_update(:,:,i) = SigmaR1(:,:,i);
                SigmaR1(:,:,i)  = SigmaR1(:,:,i) + invDn_update(:,:,i);
                tmp1_MI1        = log2(real(det(SigmaR1(:,:,i))));
                [~,On_svd1,~]   = svd(On_update(:,:,i));
                tmp2_MI1        = log2(real(det(On_update(:,:,i))));

            end
            disp(tmp_MIs_1);

            itR = itR + 1;
            
        end
        R_new = R_update;
        toc;
        
        disp(it);
        
        for i = 1:N_BS
            J(:,:,i) = J_update(:,:,i);
        end 
    
    end

elseif is_diagQ ~= 0    %%%% Scalar Quantization

    for it = 1:iter
        if it == 1
            D_update  = D_ini_diag;
            R_update  = R0;
        else
            D_update  = D_new;
            R_update  = R_new;
        end
        
        
        %%%%%% optimize quantization noise covariance: O (here optimize D = inv(T+O) )
        R_new = R_update;
        eps   = 100;
        
        J_diag = [];
        for i = 1:N_BS
            if i == 1
                J_diag = eye(L).*J(:,:,i);
            else
                J_diag = cat(3,J_diag,eye(L).*J(:,:,i));
            end
        end

        while(eps > 1e-2)
            if it == 1
                % break;
            end
    
            D_past     = D_update;
            D_current  = D_past;
            tem_SCA1   = zeros([N_BS,1]);
            tem_SCA2   = zeros([L,L,N_BS]);
            Dn_current = zeros([L,L,N_BS]);
            Sigma      = zeros([L,L,N_BS]);

            for i = 1:N_BS
                Dn_current(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
                Sigma(:,:,i)      = eye(L) + J_diag(:,:,i)^(0.5)*Dn_current(:,:,i)*(J_diag(:,:,i)^(0.5));
                tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
                tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
   

            end
            disp(tmpu);
            

            scale_D = 1e-3;
            cvx_begin
    
                variable d_opt(L*N_BS,1)
                variable tQv(2*K)
    
                D_mat = eye(L*N_BS).*(d_opt*ones([1,L*N_BS]));
    
                [~,F2] = PFIM_MCMT_new(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                   R_new,D_mat,ind_C,C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
                

                D_opt = [];
                for i = 1:N_BS
                    if i == 1
                        D_opt = I_subL(:,:,i).'*D_mat*I_subL(:,:,i);
                    else
                        D_opt = cat(3,D_opt,I_subL(:,:,i)'*D_mat*I_subL(:,:,i));
                    end
                end
                      
                T_i = diag(D0);
        
                minimize(sum(tQv))
                subject to
                    for i2K = 1:(2*K)
                        [scale_D*F2,I_Fxyb(:,i2K);I_Fxyb(i2K,:),tQv(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                    end

                    real(T_i - d_opt)*1000 >= 0;

                    %%%%%% Lemma 1 approximation
                    for i = 1:N_BS
                        (tem_SCA1(i) + real(trace(Sigma(:,:,i)\(eye(L)+J_diag(:,:,i)^(0.5)*D_opt(:,:,i)*(J_diag(:,:,i)^(0.5)))))/log(2) ...
                        - log_det(eye(L)-T_ten(:,:,i)^(0.5)*D_opt(:,:,i)*T_ten(:,:,i)^(0.5))/log(2) - L/log(2))*0.003 <= TFn(i)*0.003;
                    end

                    d_opt >= 0;
    
            cvx_end
    
            CRBF2    = inv(F2);
            disp(trace(CRBF2(1:(2*K),1:(2*K))));

            D_update = D_mat.*kron(eye(N_BS),ones([L,L]));
            eps      = norm(D_update-D_past,'fro')/(N_BS*L)/norm(D_past,'fro');
            [~,PFu]  = PFIM_MCMT_new(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                 R_new,D_update,ind_C,C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
            CRBu     = inv(PFu);
            disp(trace(CRBu(1:(2*K),1:(2*K))));

    
            for i = 1:N_BS
                Dn_current(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
                Sigma(:,:,i)      = eye(L) + J_diag(:,:,i)^(0.5)*Dn_current(:,:,i)*(J_diag(:,:,i)^(0.5));
                tem_SCA1(i)       = log2(real(det(Sigma(:,:,i))));
                tem_SCA2(:,:,i)   = eye(L)/Sigma(:,:,i);
    
            end
            disp(tmpu);
    
        end
        D_new = D_update;
        toc;


        %%%%%% optimize transmit covariance: R
        if it == iter
            break;
        end
        Dn_update = zeros(L,L,N_BS);
        for i = 1:N_BS
            Dn_update(:,:,i) = I_subL(:,:,i).'*D_update*I_subL(:,:,i);
        end
        eps_R   = 100;
        scale_R = 1e-3; %1e-3 %Nr=32,ind_C=2,choice=1,2e-3 %reference central: 8e-4
        
        while(eps_R > 1e-1)
            R_past       = R_update;
            GRG_prob     = zeros([L,L,N_BS,n_prob]);
            SigmaR       = zeros([L,L,N_BS]);
            tmp_MIs      = zeros([N_BS,1]);
            On_update    = zeros([L,L,N_BS]);
            invDn_update = zeros([L,L,N_BS]);
            diag_det     = zeros([L,N_BS]);
            for i = 1:N_BS
                On_update(:,:,i)    = eye(L)/(Dn_update(:,:,i))-T_ten(:,:,i);
                invDn_update(:,:,i) = eye(L)/(Dn_update(:,:,i));

                for ip = 1:n_prob
                    GRG_prob(:,:,i,ip) = G_prob(:,:,i,ip)'*R_past*G_prob(:,:,i,ip);
                    SigmaR(:,:,i)      = SigmaR(:,:,i) + weight_prob(ip)*GRG_prob(:,:,i,ip);
                end
                SigmaR(:,:,i)    = eye(L).*SigmaR(:,:,i) + invDn_update(:,:,i);
                [~,On_svd,~]     = svd(On_update(:,:,i));
                
                diag_det(:,i)    = diag(SigmaR(:,:,i))./diag(On_update(:,:,i));
                tmp_MIs(i)       = sum(log2(real(diag_det(:,i))));
              
            end
    
            cvx_precision best
            cvx_begin sdp
                variable Rnn(Nt,Nt,N_BS) hermitian semidefinite complex
                variable t(2*K)
                
                for i = 1:N_BS
                    if i == 1
                        R1 = I_subR(:,:,i)*Rnn(:,:,i)*I_subR(:,:,i)';
                    else
                        R1 = R1 + I_subR(:,:,i)*Rnn(:,:,i)*I_subR(:,:,i)';
                    end
                end

                [~,FR] = PFIM_MCMT_new(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                                   R1,D_update,ind_C,C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);

                J_optVar = [];
                for i = 1:N_BS
                    for ip = 1:n_prob
                        if ip == 1
                            J_optVar_iBS = weight_prob(ip)*G_prob(:,:,i,ip)'*R1*G_prob(:,:,i,ip);
                        else
                            J_optVar_iBS = J_optVar_iBS + weight_prob(ip)*G_prob(:,:,i,ip)'*R1*G_prob(:,:,i,ip);
                        end
                    end
                    J_optVar = cat(3,J_optVar,J_optVar_iBS);
                end


                minimize(sum(t))
                subject to
                    for i2K = 1:(2*K)
                        [scale_R*FR,I_Fxyb(:,i2K);I_Fxyb(i2K,:),t(i2K)] == semidefinite(2*K+2*K*N_BS*N_BS+1);
                    end

                    for iNB = 1:N_BS
                        real(trace(Rnn(:,:,iNB)))*100 <= Nt/var_n*100;
                        (real(tmp_MIs(iNB)) + real(trace(SigmaR(:,:,iNB)\(eye(L).*J_optVar(:,:,iNB)+inv(Dn_update(:,:,iNB)))))/log(2) - L/log(2))*0.01 <= TFn(iNB)*0.01; %0.001

                    end
            cvx_end

            R_update = R1;
            eps_R    = norm(R_update-R_past,'fro')/(N_BS*Nt*N_BS*Nt)/norm(R_past,'fro');
    
            tmp_cons_fr = zeros([N_BS,1]);
            for i = 1:N_BS
                tmp_cons_fr(i) = (real(tmp_MIs(i)) + real(trace(SigmaR(:,:,i)\(eye(L).*J_optVar(:,:,i)+inv(Dn_update(:,:,i)))))/log(2) - L/log(2));
            end
            disp(tmp_cons_fr);
    
            GRG_prob1 = zeros([L,L,N_BS,n_prob]);
            SigmaR1   = zeros([L,L,N_BS]);
            J_update  = zeros([L,L,N_BS]);
            for i = 1:N_BS
                On_update(:,:,i)    = eye(L)/(Dn_update(:,:,i))-T_ten(:,:,i);
                invDn_update(:,:,i) = eye(L)/(Dn_update(:,:,i));
                
                for ip = 1:n_prob
                    GRG_prob1(:,:,i,ip) = G_prob(:,:,i,ip)'*R_update*G_prob(:,:,i,ip);
                    SigmaR1(:,:,i)      = SigmaR1(:,:,i) + weight_prob(ip)*GRG_prob1(:,:,i,ip);
                end
                J_update(:,:,i) = SigmaR1(:,:,i);
                SigmaR1(:,:,i)  = eye(L).*SigmaR1(:,:,i) + invDn_update(:,:,i);
                tmp1_MI1        = log2(real(det(SigmaR1(:,:,i))));
                [~,On_svd1,~]   = svd(On_update(:,:,i));
                tmp2_MI1        = log2(real(det(On_update(:,:,i))));
    
            end
            disp(tmp_MIs);
            
        end
        R_new = R_update;
        toc;
        
        disp(it);
    
        for i = 1:N_BS
            J(:,:,i) = J_update(:,:,i);
        end
    


    end
end
%%%%%%%%%% End: Joint Optimization for N BS scenario
toc;

if is_diagQ == 0
    [~,PFo] = PFIM_MCMT(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                       R_update,D_update,ind_C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
elseif is_diagQ == 1
    [~,PFo] = PFIM_MCMT_new(N_BS,K,Nt,Nr,L,A_sens_prob,de_A_sens_prob,V_sens_prob,de_V_sens_prob,Bnu_prob,...
                       R_update,D_update,ind_C,C,tan_psens_prob,cot_psens_prob,p_diff_prob,weight_prob);
end
PCRBo      = inv(PFo);
mse_final  = trace(PCRBo(1:(2*K),1:(2*K)));
disp(mse_final);

disp('done.')