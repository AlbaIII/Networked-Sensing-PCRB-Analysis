function [FIM_full, FIM_xyb] = FIM_MCMT(N_BS,K,Nt,Nr,L,A,dA,V,dV,B,R,D,ind_C,tan_psens,cot_psens,p_diff)
%CRB_MCMT Summary of this function goes here
%   Detailed explanation goes here

AA  = [];
dAA = [];
for iN = 1:N_BS
    AA  = [AA,A(:,:,iN)];
    dAA = [dAA,dA(:,:,iN)];
end

CAA  = [];
CdAA = [];
if ind_C == 0
    % L     = Nr;
    % block   = ones(Nr,K);
    DFT_mtx = 1/sqrt(Nr)*fft(eye(Nr));
    C       = zeros(Nr,L,N_BS);
    CA      = zeros(L,K,N_BS);
    CdA     = zeros(L,K,N_BS);
    for iN = 1:N_BS
        C(:,:,iN)   = eye(L);
        CA(:,:,iN)  = C(:,:,iN)'*A(:,:,iN);
        CdA(:,:,iN) = C(:,:,iN)'*dA(:,:,iN);
        CAA         = [CAA,CA(:,:,iN)];
        CdAA        = [CdAA,CdA(:,:,iN)];
    end
elseif ind_C == 1
    % L   = 2*K;
    C   = zeros(Nr,L,N_BS);
    CA  = zeros(L,K,N_BS);
    CdA = zeros(L,K,N_BS);
    % block = ones(L,K);
    for iN = 1:N_BS
        pm_A             = (eye(Nr)-A(:,:,iN)/(A(:,:,iN)'*A(:,:,iN))*A(:,:,iN)')*dA(:,:,iN);
        C_tmp            = [A(:,:,iN),pm_A];
        [UCCH,eCCH,VCCH] = svd(C_tmp*C_tmp');
        C(:,:,iN) = UCCH(:,1:L);

        CA(:,:,iN)  = C(:,:,iN)'*A(:,:,iN);
        CdA(:,:,iN) = C(:,:,iN)'*dA(:,:,iN);
        CAA         = [CAA,CA(:,:,iN)];
        CdAA        = [CdAA,CdA(:,:,iN)];
    end
elseif ind_C == 2
    % L   = 2*K; 
    C   = zeros(Nr,L,N_BS);
    CA  = zeros(L,K,N_BS);
    CdA = zeros(L,K,N_BS);
    % block = ones(L,K);

    tmpaa1 = zeros([L,Nr,N_BS]);
    tmpaa1_max = zeros([N_BS,L]);
    for iN = 1:N_BS

        DFT_mtx = 1/sqrt(Nr)*fft(eye(Nr));
        PidA_n = (eye(Nr)-A(:,:,iN)*inv(A(:,:,iN)'*A(:,:,iN))*A(:,:,iN)')*dA(:,:,iN);
        C_tmp  = [A(:,:,iN),PidA_n];
        C_tmp1 = C_tmp ./ (ones(Nr,1)*sqrt(sum(abs(C_tmp).^2,1)));

        tmpaa1(:,:,iN) = abs(C_tmp1'*DFT_mtx);
        for im = 1:L
            [~,tmpaa1_max(iN,im)] = max(tmpaa1(im,:,iN));
        end

        if iN == 1
            % C(:,:,iN) = DFT_mtx(:,[5,6,7,8]);
            % C(:,:,iN) = DFT_mtx(:,[7,8,9,10]);    % Nr=10;
            % C(:,:,iN) = DFT_mtx(:,[12,13,14,16]);   % Nr=20;
            % C(:,:,iN) = DFT_mtx(:,[15,16,17,18]);   % Nr=20;
            C(:,:,iN) = DFT_mtx(:,[24,20,23,19]);   % Nr=30;
            C(:,:,iN) = DFT_mtx(:,[31,26,32,25]);   % Nr=40;
            C(:,:,iN) = DFT_mtx(:,[39,32,33,38]);   % Nr=50;
            C(:,:,iN) = DFT_mtx(:,[46,38,45,39]);   % Nr=60;
            C(:,:,iN) = DFT_mtx(:,[54,45,53,44]);   % Nr=70;
        elseif iN == 2
            % C(:,:,iN) = DFT_mtx(:,[3,4,5,6]);
            % C(:,:,iN) = DFT_mtx(:,[4,5,6,7]);    % Nr=10;
            % C(:,:,iN) = DFT_mtx(:,[8,9,10,11]);   % Nr=20;
            % C(:,:,iN) = DFT_mtx(:,[15,12,13,14]);   % Nr=20;
            C(:,:,iN) = DFT_mtx(:,[14,12,13,11]);   % Nr=30;
            C(:,:,iN) = DFT_mtx(:,[16,17,18,19]);   % Nr=40;
            C(:,:,iN) = DFT_mtx(:,[23,20,22,19]);   % Nr=50;
            C(:,:,iN) = DFT_mtx(:,[27,24,26,23]);   % Nr=60;
            C(:,:,iN) = DFT_mtx(:,[31,27,32,28]);   % Nr=70;
        end
        test = C(:,:,iN)'*C(:,:,iN);

        CA(:,:,iN)  = C(:,:,iN)'*A(:,:,iN);
        CdA(:,:,iN) = C(:,:,iN)'*dA(:,:,iN);
        CAA  = [CAA,CA(:,:,iN)];
        CdAA = [CdAA,CdA(:,:,iN)];
    end
elseif ind_C == 3   
    % L   = 2*K-1;
    C   = zeros(Nr,L,N_BS);
    CA  = zeros(L,K,N_BS);
    CdA = zeros(L,K,N_BS);
    for iN = 1:N_BS
        DFT_mtx = 1/sqrt(Nr)*fft(eye(Nr));
        if iN == 1
            C(:,:,iN) = DFT_mtx(:,[8,6,7]);
        elseif iN == 2
            C(:,:,iN) = DFT_mtx(:,[4,5,2]);
        end

        CA(:,:,iN)  = C(:,:,iN)'*A(:,:,iN);
        CdA(:,:,iN) = C(:,:,iN)'*dA(:,:,iN);
        CAA  = [CAA,CA(:,:,iN)];
        CdAA = [CdAA,CdA(:,:,iN)];
    end
end

block     = ones([L,K]);
ones_diag = kron(eye(N_BS),block);
AA_diag   = kron(ones(N_BS,1),CAA) .*ones_diag;
dAA_diag  = kron(ones(N_BS,1),CdAA).*ones_diag;

AHA   = AA_diag' *D*AA_diag;
dAHA  = dAA_diag'*D*AA_diag;
AHdA  = AA_diag' *D*dAA_diag;
dAHdA = dAA_diag'*D*dAA_diag;

ID = zeros(Nt*N_BS,Nt,N_BS);
I  = eye(N_BS);
for iN = 1:N_BS
    ID(:,:,iN) = kron(I(:,iN),eye(Nt));
end
IAQA = zeros(K*N_BS,K,N_BS);
for iN = 1:N_BS
    IAQA(:,:,iN) = kron(I(:,iN),eye(K));
end

%%%%%%%% FIM derivation
F11 = [];
F12 = [];
F22 = [];
% Calculate F(theta,theta) 
for iN1 = 1:N_BS
    F11r = [];
    for iN2 = 1:N_BS
        tmp1 = zeros(K,K);
        if iN1 == iN2
            for iN3 = 1:N_BS
                tmp1 = tmp1 + ...
                    (IAQA(:,:,iN3)'*AHA*IAQA(:,:,iN3)).*(conj(B(:,:,iN3,iN1))*dV(:,:,iN1)'*ID(:,:,iN1)'*conj(R)*ID(:,:,iN2)*dV(:,:,iN2)*B(:,:,iN3,iN2)) + ...
                    (IAQA(:,:,iN1)'*dAHA*IAQA(:,:,iN1)).*(conj(B(:,:,iN1,iN3))*V(:,:,iN3)'*ID(:,:,iN3)'*conj(R)*ID(:,:,iN2)*dV(:,:,iN2)*B(:,:,iN1,iN2)) + ...
                    (IAQA(:,:,iN2)'*AHdA*IAQA(:,:,iN2)).*(conj(B(:,:,iN2,iN1))*dV(:,:,iN1)'*ID(:,:,iN1)'*conj(R)*ID(:,:,iN3)*V(:,:,iN3)*B(:,:,iN2,iN3));
            
            end
            tmp2 = zeros(K,K);
            for iu = 1:N_BS
                for iv = 1:N_BS
                    tmp2 = tmp2 + conj(B(:,:,iN1,iu))*V(:,:,iu)'*ID(:,:,iu)'*conj(R)*ID(:,:,iv)*V(:,:,iv)*B(:,:,iN2,iv);
                end
            end
            tmp2 = tmp2.*(IAQA(:,:,iN1)'*dAHdA*IAQA(:,:,iN1));   %iN1 == iN2
            tmp1 = tmp1 + tmp2;
        else
            for iN3 = 1:N_BS
                tmp1 = tmp1 + ...
                    (IAQA(:,:,iN3)'*AHA*IAQA(:,:,iN3)).*(conj(B(:,:,iN3,iN1))*dV(:,:,iN1)'*ID(:,:,iN1)'*conj(R)*ID(:,:,iN2)*dV(:,:,iN2)*B(:,:,iN3,iN2)) + ...
                    (IAQA(:,:,iN1)'*dAHA*IAQA(:,:,iN1)).*(conj(B(:,:,iN1,iN3))*V(:,:,iN3)'*ID(:,:,iN3)'*conj(R)*ID(:,:,iN2)*dV(:,:,iN2)*B(:,:,iN1,iN2)) + ...
                    (IAQA(:,:,iN2)'*AHdA*IAQA(:,:,iN2)).*(conj(B(:,:,iN2,iN1))*dV(:,:,iN1)'*ID(:,:,iN1)'*conj(R)*ID(:,:,iN3)*V(:,:,iN3)*B(:,:,iN2,iN3));
            
            end
        end
        F11r = [F11r,tmp1];
    end
    F11 = [F11;F11r];
end

% Calculate F(theta,b)
for iNth = 1:N_BS
    F12r = [];
    for iNb1 = 1:N_BS
        for iNb2 = 1:N_BS
            tmp3 = zeros(K,K);
            if iNth ~= iNb1
                tmp3 = tmp3 + (IAQA(:,:,iNb1)'*AHA*IAQA(:,:,iNb1)).*(conj(B(:,:,iNb1,iNth))*dV(:,:,iNth)'*ID(:,:,iNth)'*conj(R)*ID(:,:,iNb2)*V(:,:,iNb2));

            else
                for iN4 = 1:N_BS
                    tmp3 = tmp3 + ...
                        (IAQA(:,:,iNth)'*dAHA*IAQA(:,:,iNth)).*(conj(B(:,:,iNth,iN4))*V(:,:,iN4)'*ID(:,:,iN4)'*conj(R)*ID(:,:,iNb2)*V(:,:,iNb2));
                        
                end
                tmp3 = tmp3 + (IAQA(:,:,iNth)'*AHA*IAQA(:,:,iNth)).*(conj(B(:,:,iNth,iNth))*dV(:,:,iNth)'*ID(:,:,iNth)'*conj(R)*ID(:,:,iNb2)*V(:,:,iNb2));
            end
            F12r = [F12r,tmp3];
        end
    end
    F12 = [F12;F12r];
end

% Calculate F(b,b)
for iNb1 = 1:N_BS
    for iNb2 = 1:N_BS
        F22_r = [];
        for iNb3 = 1:N_BS
            for iNb4 = 1:N_BS
                if iNb1 ~= iNb3
                    tmp4 = zeros(K,K);
                else
                    tmp4 = (IAQA(:,:,iNb1)'*AHA*IAQA(:,:,iNb3)) .* (V(:,:,iNb2)'*ID(:,:,iNb2)'*conj(R)*ID(:,:,iNb4)*V(:,:,iNb4)); % iNb1 == iNb3
                    
                end
                F22_r = [F22_r,tmp4];
            end
        end
        F22 = [F22;F22_r];
    end
end
FIM_full = 2*[real(F11),    real(F12),   -imag(F12);...
              real(F12).',  real(F22),   -imag(F22);...
             -imag(F12).', -imag(F22).',  real(F22)];

T_tq = zeros(2*K,N_BS*K);
for iN = 1:N_BS
    T_tqn = zeros(2*K,K);
    for k = 1:K
        T_txn =  1./((1+tan_psens(k,iN).^2).*p_diff(k,2,iN));
        T_tyn = -1./((1+cot_psens(k,iN).^2).*p_diff(k,1,iN));
        K_rg  = ((k-1)*2+1):(k*2);
        T_tqn(K_rg,k) = [T_txn;T_tyn];
    end
    T_tq(:,((iN-1)*K+1):(iN*K)) = T_tqn;
end
T = [T_tq,zeros([2*K,2*N_BS*N_BS*K]);zeros([2*N_BS*N_BS*K,K*N_BS]),eye(2*N_BS*N_BS*K)];

FIM_xyb = T*FIM_full*T.';



end

