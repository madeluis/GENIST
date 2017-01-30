function [d_long_scn,reg,set_max_reg,change_t,change_t_up,change_t_down,n_genes,n_levels] = GRN_preprocessing(long_scn,autoregulation)

    %% BUILD A BAYESIAN NETWORK OF THIS SUBSET OF GENES

    % Assumptions:
    % 1) The only genes that can regulate another gene are the ones that have a
    % change in the PREVIOUS or simultaneous time point. It makes sense in our
    % case!!
    % 2) At each gene, I look only at the
    % first change (to be a target), but the regulators of that gene can come
    % from changes that don't correspond to the first change in the other
    % genes.

    %% Estimate potential regulators
    % IN THIS BLOCK, I ESTIMATE THE POTENTIAL REG BASED ON THE ZOU & CONZEN
    % PAPER, BUT INSTEAD OF, FOR GENE i, SELECT ALL THE GENES THAT HAD A CHANGE B4 THE
    % 1ST CHANGE OF GENE i, I SELECT ALL THE GENES THAT HAD A CHANGE RIGHT B4 (PREVIOUS
    % TIME POINT) THE 1ST CHANGE OF GENE i 

    % I ESTIMATE THE POTENTIAL A BIT DIFFERENTLY: 
    % FOR GENE i, I SELECT ALL THE GENES THAT HAD A CHANGE RIGHT B4 (PREVIOUS
    % TIME POINT) ANY CHANGE OF GENE i

    T = size(long_scn,2);
    n_genes = size(long_scn,1);

    change_t = zeros(size(long_scn)); %var to store the positions where there is a 1.2 or -0.7 fold change 

    for i = 2:T
        change_t(:,i) = (long_scn(:,i) > 1.1*long_scn(:,i-1)) | (long_scn(:,i) < 0.9*long_scn(:,i-1));
        change_t_up(:,i) = (long_scn(:,i) > 1.1*long_scn(:,i-1));
        change_t_down(:,i) = (long_scn(:,i) < 0.9*long_scn(:,i-1));
    end

     % variable to store potential regulators and corregulators. Contains two 
     % 2-D matrices. Each mat: columns correspond to variables, so column i
     % contains all the (co)regulators of gene i. First mat stores regulators.
     % Second mat stores corregulators.
    reg = zeros(n_genes,n_genes,2);

    for i = 1:n_genes
        for j = 2:T
            if change_t(i,j) == 1 
                reg(:,i,1) = reg(:,i,1) == 1 | change_t(:,j-1) == 1; %store potential regulators
                reg(:,i,2) = reg(:,i,2) == 1 | change_t(:,j) == 1; %store potential coregulators
            end
        end
    end

    if autoregulation == 0
        % remove diagonal from reg matrix. i.e., dont let genes autoregulate
        for i = 1:2
            D = diag(reg(:,:,i));
            reg(:,:,i) = reg(:,:,i) - diag(D);
        end
    end
    


    %% Discretize original expression data in n_levels levels

    n_levels = 2;
    if n_levels == 2
        m = mean(long_scn,2) * ones(1,T);
        %d_long_scn = zeros(size(long_scn));
        d_long_scn = long_scn >= m;
    elseif n_levels == 3
        q = quantile(long_scn,2,2); % quantile(X,N,dim)
        for i = 1:T
            for j = 1:n_genes
                if long_scn(j,i) < q(j,1)
                    d_long_scn(j,i) = 0;
                elseif long_scn(j,i) > q(j,1) & long_scn(j,i) < q(j,2)
                    d_long_scn(j,i) = 1;
                else
                    d_long_scn(j,i) = 2;
                end
            end
        end
    elseif n_levels == 4
        q = quantile(long_scn,3,2); % quantile(X,N,dim)
        for i = 1:T
            for j = 1:n_genes
                if long_scn(j,i) < q(j,1)
                    d_long_scn(j,i) = 0;
                elseif long_scn(j,i) > q(j,1) & long_scn(j,i) < q(j,2)
                    d_long_scn(j,i) = 1;
                elseif long_scn(j,i) > q(j,2) & long_scn(j,i) < q(j,3)
                    d_long_scn(j,i) = 2;
                else
                    d_long_scn(j,i) = 3;
                end
            end
        end
    end


    max_reg = max(max(sum(reg)));% var to store the max # of potential regulators that have been found for one gene
    if max_reg > 3
        set_max_reg = 3; % I set the maximum number of possible regulators of a gene to 3.
    else
        set_max_reg = max_reg; 
    end

end