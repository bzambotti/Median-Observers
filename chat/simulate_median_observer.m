function [x_hat,e_norm,selected_delay,mode_history,...
          innovation_raw,innovation_filtered] = ...
    simulate_median_observer(A,B,C,L,x,u,y,q)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Robust Median Observer Simulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(A,1);
p = size(C,1);

sim_steps = size(y,2);

N = 2*q+1;

modes = generate_modes(q,p);
num_modes = size(modes,1);

% Gains

if ~iscell(L)

    Lcell = cell(num_modes,1);

    for i=1:num_modes
        Lcell{i}=L;
    end

else

    Lcell=L;

end

% Allocation

x_hat=zeros(n,sim_steps+1);

e_norm=zeros(1,sim_steps+1);

innovation_raw=zeros(p,sim_steps);

innovation_filtered=zeros(p,sim_steps);

selected_delay=zeros(p,sim_steps);

mode_history=zeros(1,sim_steps);

% Large initial estimation error

x_hat(:,1)=x(:,1)+50*randn(n,1);

e_norm(1)=norm(x(:,1)-x_hat(:,1));

% Innovation memory

sigma=zeros(p,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:sim_steps

    % Raw innovation

    raw=y(:,k)-C*x_hat(:,k);

    innovation_raw(:,k)=raw;

    % Update memory

    sigma=[raw sigma(:,1:end-1)];

    % Median computation

    med=zeros(p,1);

    active=zeros(1,p);

    for j=1:p

        s=sort(sigma(j,:));

        med(j)=s(q+1);

        idx=find(sigma(j,:)==med(j),1);

        active(j)=idx-1;

    end

    innovation_filtered(:,k)=med;

    selected_delay(:,k)=active(:);

    %Determine active mode

    mode=1;

    for i=1:num_modes

        if isequal(modes(i,:),active)

            mode=i;
            break

        end

    end

    mode_history(k)=mode;

    % Observer gain

    Lk=Lcell{mode};

   % Observer update

    x_hat(:,k+1)=A*x_hat(:,k)+ ...
                 B*u(:,k)+ ...
                 Lk*med;

    % Error

    e_norm(k+1)=norm(x(:,k+1)-x_hat(:,k+1));

end

end