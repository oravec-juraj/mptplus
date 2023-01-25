function [iMPC,eMPC,secOutputs] = TMPCController(model,N,option)
% -------------------------------------------------------------------------
%          [iMPC,eMPC,secOutputs] = TMPCController(model,N,option)
% -------------------------------------------------------------------------
%  Function computes a Tube MPC policy based on given uncertain model,
%  prediction horizon, and selected options. For more details we refere
%  users to read:
%  Tube MPC Extension of MPT: Experimental Analysis by Juraj Holaza, 
%  Erika Pavlovi\v{c}ov\'{a}, Juraj Oravec. Proceedings of the 24th 
%  International Conference on Process Control, June 2023
% ------------------------------------------------------------------------
% INPUTS:
%      model                - MPT3 object with ULTISystem 
%      N                    - Prediction horizon
% OPTIONS:
%      eMPC                 - {0/1} defines if the explicit MPC should be
%                             calculated (0 - default value)
%      adjList              - {0/1} computes adjecency list (default = 0)
%      Tube_tol             - allowed tollerace (default = 1e-4 )
%      Tube_MaxIter         - maximal number of itterations (default = 1e3)
%      solType              - {0/1} defines the retuned optimized values 
%                             (the default value is 0):
%                             0: iMPC/eMPC returns TMPC outputs: 
%                                           [u_opt, x_opt] 
%                             1: iMPC/eMPC returns the transformed TMPC 
%                                outputs: 
%                                     Utube = u_opt + K*( x - x_opt ),
%                                where [u_opt, x_opt] are the original TMPC
%                                outputs, K is the LQR for the given 
%                                model (given as u = K*x), and x is the  
%                                state measurement.
%
% OUTPUTS:
%      iMPC                 - online MPC policy (MPT3 object)
%      eMPC                 - explicit MPC policy (MPT3 object)
%                             (computed only if requested by 'option.eMPC')
%      secOutputs.Tube      - approximated minimal robust invariant set
%      secOutputs.problem   - defined optimization problem (opt object)
%      secOutputs.Xc        - state constraints robustification
%      secOutputs.Uc        - input constraints robustification
%      secOutputs.K         - LQR controller (given as u = K*x)
%      secOutputs.mpsol     - eMPC as MPT2 object (retuend by YALMIP)
%      secOutputs.DIAGNOSTIC- further outputs from the YALMIP
%      secOutputs.Z         - further outputs from the YALMIP
%      secOutputs.HPWF      - further outputs from the YALMIP
%      secOutputs.ZPWF      - further outputs from the YALMIP
% ------------------------------------------------------------------------
% % EXAMPLE:
% % LTI system
% model = ULTISystem('A', [1, 1; 0, 1], 'B', [0.5; 1], 'E', [1, 0; 0, 1]);
% model.u.min = [-1]; 
% model.u.max = [ 1];
% model.x.min = [-5; -5]; 
% model.x.max = [ 5;  5];
% model.d.min = [-0.1; -0.1]; 
% model.d.max = [ 0.1;  0.1];
% % Penalty functions
% model.x.penalty = QuadFunction(diag([1, 1]));
% model.u.penalty = QuadFunction(diag([0.01]));
% % Terminal set
% model.x.with('terminalSet');
% model.x.terminalSet = model.LQRSet;
% % Terminal penalty
% model.x.with('terminalPenalty');
% model.x.terminalPenalty = model.LQRPenalty;
% % Prediction horizon
% N = 2;
% % Enable construction of the explicit controller
% option = {'eMPC',1};
% % TMPC controller construction
% [iMPC,eMPC,sol] = TMPCController(model,N,option)
% % TMPC evaluation
% x0 = [-5; 2]; % Initial condition
% u_implicit = iMPC.optimizer(x0) % Implicit MPC evaluation
% u_explicit = eMPC.evaluate(x0)  % Explicit MPC evaluation
% [u, feasible, openloop] = eMPC.evaluate(x0) % Enriched output
% ------------------------------------------------------------------------

%% Checks / initial setup
% Perform checks:
if nargin < 2, error('Ups, not enough inputs were selected'); end
% We do not support parametric uncertainties
if iscell(model.A) || iscell(model.B)
    error('Implemented Tube MPC design method does not support parametric uncertainties (LPV systems), yet.')
end


% set options/default values
ip = inputParser;
ip.addParameter('eMPC', 0);
ip.addParameter('adjList', 0);
ip.addParameter('Tube_tol', 1e-4);
ip.addParameter('Tube_MaxIter', 1e3);
ip.addParameter('solType', 1);
if nargin == 3, ip.parse(option{:}); else ip.parse(); end
options = ip.Results;


%% TMPC setup

% Bounded additive uncertainty 
W = Polyhedron('lb',model.d.min,'ub',model.d.max);
W = model.E*W;
W = Polyhedron('lb',-max(abs(W.V)),'ub',max(abs(W.V)));

% Constraints handling
Sx = Polyhedron('lb',model.x.min,'ub',model.x.max);
Su = Polyhedron('lb',model.u.min,'ub',model.u.max);

% Robust positive invariant set
[Klqr, ~] = dlqr(model.A, model.B, model.x.penalty.weight, model.u.penalty.weight); % LQR controller
K = -Klqr;
Tube = approx_mRPIS(model.A,model.B,K,W,options.Tube_tol,options.Tube_MaxIter,options.Tube_approx);
if Tube.isEmptySet
    error('Cannot continue as the minimal positive invariant set is empty!')
end

% State constraints robustification
Xc = Sx - Tube;
if Xc.isEmptySet
    error('Cannot continue as the state constraint set is empty!')
end
% Input constraints robustification
Uc = Su - ( K*Tube );
if Uc.isEmptySet
    error('Cannot continue as the input constraint set is empty!')
end

%% Construction of the Tube eMPC

% Define decision variables - YALMIP variables
X0 = sdpvar(model.nx,1);
u = sdpvar(model.nu,N,'full');
x = sdpvar(model.nx,N+1,'full');
if isprop(model.x,'reference'), xref = sdpvar(model.nx,1,'full'); end
if isprop(model.u,'deltaPenalty'), uprev = sdpvar(model.nu,1,'full'); end   

% Initialize constraints and objective function
con = [];
obj = 0;

% Constraints for OPTIMIZER
con = [con, Sx.A*X0 <= Sx.b ];
con = [con, Tube.A*( X0 - x(:,1) ) <= Tube.b ];

for k = 1:N
    % Objective - state
    if isprop(model.x,'reference')
        error('MPTplus does not support Tube MPC for reference tracking problem, yet.')
    else
        obj = obj + x(:,k)'*model.x.penalty.weight*x(:,k);
    end
    % Objective - inputs
    if isprop(model.u,'deltaPenalty')
        error('MPTplus does not support Tube MPC with delta-u penalization, yet.')
    else
        obj = obj + u(:,k)'*model.u.penalty.weight*u(:,k);
    end

    % Nominal LTI model
    con = [con, x(:,k+1) == model.A*x(:,k) + model.B*u(:,k) ];

    % Input and state constraints
    con = [con, Xc.A*x(:,k) <= Xc.b ];
    con = [con, Uc.A*u(:,k) <= Uc.b ];
end

% Terminal penalty
if isprop(model.x,'terminalPenalty')
    obj = obj + x(:,end)'*model.x.terminalPenalty.weight*x(:,end);
end

% Terminal constraint
if isprop(model.x,'terminalSet')
    con = [con, model.x.terminalSet.A*x(:,end) <= model.x.terminalSet.b ];
end

% Optimization settings - silent verbose mode
opt_settings = sdpsettings('verbose',0);

p_required = X0;

if options.solType
    % In this case all controllers return: Utube = u_opt + K*( x - x_opt )
    p_return = u(:,1) + K*( X0 - x(:,1)); 
    problem = Opt(con, obj, p_required, p_return);
    
    % Implicit (non-explicit) MPC: iMPC
    iMPC = MPCController(model,1);
    iMPC.construct();
    iMPC.optimizer = optimizer(con,obj,opt_settings,p_required,p_return); 

    % Explicit MPC: eMPC
    if options.eMPC || options.adjList
        % Note that here the Utube outout will be corrected inside of 
        % function ConvertSolutionToMPT3, hence we use the [u_opt, x_opt]
        p_return = [u(:,1); x(:,1)];
        [mpsol, DIAGNOSTIC,Z,HPWF,ZPWF] = solvemp(con, obj, opt_settings, p_required, p_return);
        if cellfun(@isempty,mpsol), error('Ups, the problem is infeasible!');end
        eMPC = ConvertSolutionToMPT3(model,mpsol,K,options);
    else
        eMPC = [];
    end
else
    % In this case all controllers return: [u_opt, x_opt]
    % Note that the true closed-loop input is defined as: 
    % Utube = u_opt + K*( x - x_opt ))
    p_return = [u(:,1); x(:,1)];
    problem = Opt(con, obj, p_required, p_return);

    % iMPC
    iMPC = MPCController(model,N);
    iMPC.optimizer = optimizer(con,obj,opt_settings,p_required,p_return);
    
    % eMPC
    if options.eMPC || options.adjList
        [mpsol, DIAGNOSTIC,Z,HPWF,ZPWF] = solvemp(con, obj, opt_settings, p_required, p_return);
        if cellfun(@isempty,mpsol), error('Ups, the problem is infeasible!');end
        eMPC = ConvertSolutionToMPT3(model,mpsol,K,options);
    else
        eMPC = [];
    end
end

%% Secondary outputs:
secOutputs.problem = problem;
secOutputs.Tube = Tube;
secOutputs.Xc = Xc;
secOutputs.Uc = Uc;
secOutputs.K = K;
if options.eMPC || options.adjList
    secOutputs.mpsol = mpsol;
    secOutputs.DIAGNOSTIC = DIAGNOSTIC;
    secOutputs.Z = Z;
    secOutputs.HPWF = HPWF;
    secOutputs.ZPWF = ZPWF;
end
end


function [eMPC] = ConvertSolutionToMPT3(model,mpsol,K,option)
% Fucntion converts solution from YALMIP (MPT2 object) into the MPT3
% framework. 
% The final controller eMPC returns the true closed-loop inputs
%                   Utube = u_opt + K*( x - x_opt )
% where [u_opt, x_opt] is solution from mpsol, K is the LQR for the model,
% and x is the state measurement.
%
% model  - MPT3 object
% mpsol  - parametric solution of MPC problem in YALMIP
% option - various options
% K      - LQR policy (as u = K*x)

if option.solType
    % A) Update the feedback policy
    % Primal function: u =  Fi*x + gi
    % Tube MPC primal function: Utube = u_opt + K*( x - x_opt ), where:
    %   x_opt - computed initial state of the prediction (encoded in u)
    %   u_opt - computed control inputs
    %   K     - assumed LQR policy (from the TMPC setup)
    %   u     - u = [u_opt; x_opt]
    % ------------------------------------------------------
    % Primal function: Utube = (Fu_opt + K - K*Fx_opt)*x + (gu_opt-K*gx_opt)
    indx_x_opt = model.nu+1:model.nu+model.nx;
    indx_u_opt = 1:model.nu;
    T = mpsol{1};
    for k = 1:length(T.Fi)
        % F = (Fu_opt + K - K*Fx_opt)
        mpsol{1}.Fi{k} = mpsol{1}.Fi{k}(indx_u_opt,:) + ...
            K - K*mpsol{1}.Fi{k}(indx_x_opt,:);
        % g = (gu_opt - K*gx_opt)
        mpsol{1}.Gi{k} = mpsol{1}.Gi{k}(indx_u_opt,:) - ...
            K*mpsol{1}.Gi{k}(indx_x_opt,:);
    end
end

% B) Convert the solutino into MPT3 object
eMPC = EMPCController(mpt_mpsol2pu(mpsol));

% C) Recover lost information
% eMPC.N = N;
eMPC.N = 1; % XXX: Do we need this for: solType == 1
eMPC.model = model;

% D) Change the number of inputs (needed for '.evaluate')
if ~option.solType
    eMPC.nu = model.nu + model.nx;
end

% E) Create adjacency list 
if  option.adjList == 1
    adj_list = cell(1,eMPC.nr);
    for k = 1:eMPC.nr
        ShowProgress('Creating adjacency list ... ', k, eMPC.nr)
        t = 0;
        temp = cell(1,1);
        for kk = 1:eMPC.nr
            if eMPC.partition.Set(k).isAdjacent(eMPC.partition.Set(kk))
                intersection = eMPC.partition.Set(k) & eMPC.partition.Set(kk);
                if size(intersection.V,1) >= eMPC.partition.Set(k).Dim
                    t = t+1;
                    temp{t,1} = kk;
                end
            end
        end
        adj_list{k} = temp;
    end
    % adj_list = MPT_verify_graph(eMPC.partition.Set,adj_list);
    eMPC.optimizer.setInternal('adj_list', adj_list);
end
end


function ShowProgress(operation, k, kf)
% Inputs:
%          operation - sting
%          k         - current iteration
%          kf        - final number of iterations
%          ttime     - time of the previous operation

if k == 1, fprintf(1, strcat(operation,'       :')); end
    
kkf = fix(k/kf*100);
if kkf < 10
    fprintf(1, ' \b\b\b\b%d %%', kkf);
else
    fprintf(1, ' \b\b\b\b\b%d %%', kkf);
end

if k == kf
    fprintf('\n');
end
end



function F_alpha = approx_mRPIS(A,B,K,W,tol,MaxIter,Approx)
% Invariant Approximations of the Minimal Robust Positively Invariant Set
% by S. V. Rakovic, E. C. Kerrigan, K. I. Kouramas, D. Q. Mayne
% IEEE TRANSACTIONS ON AUTOMATIC CONTROL, VOL. 50, NO. 3, MARCH 2005
% URL: https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=1406138

if nargin < 4, error('Ups, not enough inputs...'); end
if nargin < 5, tol = 1e-4; end
if nargin < 6, MaxIter = 1000; end
if nargin < 7, Approx = 0; end

% Closed-loop/autonomous system matrix 
Acl = A + B*K;

% Problem size
nx = size(Acl, 1); 

% Initial values
s = 0; % "s" is integer going toward zero: s -> 0
alpha = 1e6; % Initialize "alpha"
Ms = 1000;
MSS = zeros(2*nx,MaxIter);

fprintf('Constructing approximated minimal robust invariant set ...\n')
tic,
while( alpha > tol/( tol + Ms ) ) && ( s < MaxIter )

    % Increment "s"
    s = s + 1;

    % Eq(10): (Acl^s)*W \subseteq alpha*W == h_W( (Acl^(s))'*f_i ) <= alpha*g_i
    % Eq(11): alpha(s) = max_{i \in I} ( ( h_W( (Acl^(s))'*f_i ) ) / g_i )
    alpha = max(W.support(Acl^(s)*(W.A)')./W.b);

    % Eq(13): M(s) = max_{j=1,...,n} ( sum_{i=0}^{s-1} ( h_W( (Acl^(i))'*e_j ) ), sum_{i=0}^{s-1} ( h_W( -(Acl^(i))'*e_j ) ) )
    MSS(:,s) = W.support([Acl^(s), -Acl^(s)]);
    mss = sum(MSS(:,1:s),2);
    Ms = max(mss);

    temp = alpha - tol/( tol + Ms );
    if mod(s,10) == 0
        fprintf('\tConvergance to 0: %.5f\n',max(0,temp))
    end
end

% Eq(2): Fs = Minkowski_Sum_{i=0}^{s-1}( Acl^(i)*W ), where F_0 = {0}
Fs = W; 
for i = 1 : (s - 1)
    ShowProgress('\tFinalizing ...         ',i,s-1)
    Fs = Fs + Acl^(i)*W;
    Fs.minVRep();
end
F_alpha = (1 - alpha)^(-1)*Fs;
compTime = toc;
fprintf('... done! (computation time: %.2f, #itterations: %.0f, #halfspaces = %.0f)\n',compTime,s,size(Fs.H,1))

end 