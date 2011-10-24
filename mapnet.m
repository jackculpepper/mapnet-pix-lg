

for t = 1:num_trials

    tic

    loaded_chunk = 0;

    X0r = zeros(J,Bmini);
    X1r = zeros(J,Bmini);

    for b = 1:Bmini
        vtest = 0;
        j = 0;

        while var(vtest) < lower_var_thresh
            %% choose a random batch element
            i = ceil(rand*B);

            %% choose a random point in time
            f = ceil(rand*(T-1));

            X0r(:,b) = Xr(:,f,i);
            X1r(:,b) = Xr(:,f+1,i);

            vtest = X1r(:,b)-X0r(:,b);

            j = j+1;
        end

        b_log(1, b, mod(update, b_log_length)+1) = i;
        b_log(2, b, mod(update, b_log_length)+1) = f;
    end


    time_choose = toc;

    tic



    c = zeros(N,Bmini);

    for Lscn = Lmin:5:L
        EI = zeros(Lscn,Bmini);
        E = zeros(Lscn,Bmini);
        snr_b = zeros(Bmini,1);

        %% pca project
        X0 = V(:,1:Lscn)'*X0r; X0 = diag(sqrt(1./diag(D(1:Lscn,1:Lscn)))) * X0;
        X1 = V(:,1:Lscn)'*X1r; X1 = diag(sqrt(1./diag(D(1:Lscn,1:Lscn)))) * X1;
        Xdiff = X1-X0;

        for b = 1:Bmini

            switch inf_type
                case 'lbfgsb'
        
                    %% no bounds
                    lb  = zeros(1,N); % lower bound
                    ub  = zeros(1,N); % upper bound
                    nb  = zeros(1,N); % bound type (none)

                    c0 = c(:,b);
 
                    [c1,fx,exitflag,userdata] = lbfgs(@objfun_c, c0(:), ...
                                                      lb, ub, nb, ...
                                                      opts_lbfgs_c, ...
                                                      psi(1:Lscn,1:Lscn,:), ...
                                                      X0(:,b), X1(:,b), ...
                                                      V(:,1:Lscn), lambda, mask);
                    c(:,b) = reshape(c1, N, 1); 
            end


            ExpA = expm(reshape(reshape(psi(1:Lscn,1:Lscn,:), Lscn^2, N)*c(:,b), Lscn, Lscn));
            if length(find(isnan(ExpA(:)))) > 0
                %% solver exploded; zero out so we can skip it in learning
                c(:,b) = 0;
                ExpA = eye(Lscn);
            end


            EI(:,b) = ExpA*X0(:,b);
            E(:,b) = X1(:,b) - EI(:,b);

            snr_b(b) = 10 * log10 ( sum(sum( (mask.*(V(:,1:Lscn)*X1(:,b))).^2)) ...
                / sum(sum((mask.*(V(:,1:Lscn)*(EI(:,b)-X1(:,b)))).^2)) );

            fprintf(' b %d/%d L %d/%d snr_b %.4f', b, Bmini, Lscn, L, snr_b(b));

        end



        if (display_every == 1 || mod(update,display_every) == 0) && ...
            mod(Lscn, 10) == 0

            sfigure(5); clf;
            sfigure(6); clf;
            sfigure(7); clf;
            %sfigure(8); clf;

            sc = max(max(abs([V(:,1:Lscn)*X0 ; V(:,1:Lscn)*X1 ; V(:,1:Lscn)*EI])));
            for b = 1:Bmini
                sfigure(5);
                subp(Bmini/5,5,b);
                    imagesc(reshape(mask.*(V(:,1:Lscn)*E(:,b)),Jsz,Jsz),[-sc sc]);
                    colormap(gray); axis image off; title('error');
                    colorbar;

                sfigure(6);
                subp(Bmini/5,5,b);
                    imagesc(reshape(V(:,1:Lscn)*X0(:,b),Jsz,Jsz),[-sc sc]), ...
                    colormap(gray); axis image off; title('frame 0');
                    colorbar;

                sfigure(7);
                subp(Bmini/5,5,b);
                    imagesc(reshape(V(:,1:Lscn)*EI(:,b),Jsz,Jsz),[-sc sc]), ...
                    colormap(gray); axis image off; title('frame 1 recon');
                    colorbar;

%                sfigure(8); colormap(gray);
%                psi_c = zeros(Lscn,Lscn);
%                for n = 1:N
%                    psi_c = psi_c + psi(1:Lscn,1:Lscn,n) * c(n,b);
%                end
%                V_psi_c_V = V(:,1:Lscn)*psi_c*V(:,1:Lscn)';
%                array_psi_c_b = render_2d_linkage(permute(V_psi_c_V, [2 1 3]));
%                subp(Bmini/5,5,b);
%                    imagesc(array_psi_c_b); axis image off;

            end
            drawnow
        end

        fprintf(' %d/%d', b, Bmini);
    end
    fprintf('\n');


    time_inf = toc;
 
    snr = mean(snr_b);
    
    if display_every == 1 || mod(update,display_every) == 0
        i = 1;

        sfigure(3);

        subplot(5,1,1),bar(X0(:,i)),title('X0 pixels'); axis tight;
        subplot(5,1,2); bar(X1(:,i)); title('X1 pixels'); axis tight;
        subplot(5,1,3); bar(EI(:,i)); title('X1 reconstructed'); axis tight;
        subplot(5,1,4); bar(E(:,i)); title('error'); axis tight;
        subplot(5,1,5); bar(c(:,i)); title('c'); axis tight;

%        sfigure(4);
%        sc = max(abs([V(:,1:L)*X0(:,i) ; V(:,1:L)*X1(:,i) ; V(:,1:L)*(X0(:,i)+EI(:,i))]));
%        subplot(2,2,1),imagesc(reshape(V(:,1:L)*X0(:,i),Jsz,Jsz),[-sc sc]), ...
%            colormap(gray),axis image off,title('frame 0');
%            colorbar;
%        subplot(2,2,2),imagesc(reshape(V(:,1:L)*X1(:,i),Jsz,Jsz),[-sc sc]), ...
%            colormap(gray),axis image off,title('frame 1');
%            colorbar;
%        subplot(2,2,3),imagesc(reshape(mask.*(V(:,1:L)*E(:,i)),Jsz,Jsz),[-sc sc]), ...
%            colormap(gray),axis image off,title('error');
%            colorbar;
%        subplot(2,2,4),imagesc(reshape(V(:,1:L)*(X0(:,i)+EI(:,i)),Jsz,Jsz),[-sc sc]), ...
%            colormap(gray),axis image off,title('frame 1 recon');
%            colorbar;

        drawnow;
    end

    tic

    % update bases

    switch lrn_type
        case 'gd'
            psi0 = psi(1:L,1:L,:);
            dpsi = zeros(L,L,N);

            fpsi0 = 0;
            for b = 1:Bmini
                if length(find(c(:,b))) > 0
                    [fpsi0_b, dpsi_b] = objfun_psi(psi0(:), c(:,b), X0(:,b), ...
                                                   X1(:,b), V(:,1:L), mask);
                    dpsi = dpsi + reshape(dpsi_b, L, L, N);
                    fpsi0 = fpsi0 + fpsi0_b;
                else
                    fprintf('warning: skipping %d due to nan\n', b);
                end
            end

            psi1 = psi0 - eta * dpsi;

            fpsi1 = 0;
            for b = 1:Bmini
                if length(find(c(:,b))) > 0
                    fpsi1_b = objfun_psi(psi1(:), c(:,b), X0(:,b), ...
                                         X1(:,b), V(:,1:L), mask);
                    fpsi1 = fpsi1 + fpsi1_b;
                end
            end
    end

    %% pursue a constant change in angle
    angle_psi = acos(psi1(:)' * psi0(:) / sqrt(sum(psi1(:).^2)) / sqrt(sum(psi0(:).^2)));
    if angle_psi < target_angle
        eta = eta*1.01;
    else
        eta = eta/1.01;
    end

    if angle_psi < angle_thresh
        psi(1:L,1:L,:) = reshape(psi1, L, L, N);
    else
        fprintf('warning: angle too big; not taking update\n');
    end


    %% normalize
    for n = 1:N
        psi(:,:,n) = sqrt(L)*psi(:,:,n)/sqrt(sum(sum(psi(:,:,n).^2)));
    end


    c_log(:, :, mod(update, c_log_length)+1) = c;


    time_updt = toc;

    % display
    
    tic
    if display_every == 1 || mod(update,display_every) == 0

        % Display the bfs
        V_psi_V = zeros(J,J,N);
        for n = 1:N
            V_psi_V(:,:,n) = V(:,1:L)*psi(1:L,1:L,n)*V(:,1:L)';
        end
        array = render_2d_linkage(permute(V_psi_V,[2 1 3]));
 
        sfigure(1); colormap(gray);
        imagesc(array, [-1 1]);
        axis image off;

        % histogram of coefficients   
        sfigure(8);
        for n = 1:N
            subplot(Nsz, Nsz, n);
            c_log_n = c_log(n, :, :);
            hist(c_log_n(:), 100);
        end


        if mod(update,save_every) == 0
            [sucess,msg,msgid] = mkdir(sprintf('state/%s', paramstr));

            array_frame = uint8(255*((array+1)/2)+1);
 
            imwrite(array_frame, ...
                sprintf('state/%s/psi_up=%06d.gif',paramstr,update), 'gif');
        end
        drawnow;
    end
    time_disp = toc;


    fprintf('%s update %d ch %.2f if %.2f ud %.2f dy %.2f', ...
        paramstr,update,time_choose,time_inf,time_updt,time_disp);
    fprintf(' ang %.4f', angle_psi);
    fprintf(' eta %.4f', eta);
    fprintf(' snr %.4f min_snr %.4f\n', mean(snr_b), min(snr_b));

    update = update + 1;
end

savestate

