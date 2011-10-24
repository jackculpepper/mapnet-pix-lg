

switch datasource
    case 'movies'
        [sucess,msg,msgid] = mkdir(sprintf('../cache'));

        Wsz = 128;
        T = 64;
        topmargin = 15;
        num_chunks = 56;
        load_interval = 5;
        buff = 4;

        data_root = sprintf('../../data/vid075-%s', datatype);

        filename_data = sprintf('../cache/Xr_vid075_%s_%dx%dx%d.mat', datatype, J, T, B);
        if exist( filename_data, 'file')
            cache = load(filename_data);
            Xr = cache.Xr;
            %Zr = cache.Zr;
            V = cache.V;
            D = cache.D;
        else
            Xr = zeros(J,T,B);
            Zr = zeros(J,2,T-1,B);
            fprintf('selecting image patches ..\n');

            for b = 1:B
                if ~exist('F','var') || mod(b,load_interval) == 0
                    % choose a movie for this batch
                    %% skip the first chunk, which is used for testing
                    j = 1 + ceil((num_chunks-1)*rand);
                    j = 24;
                    [F,G] = read_chunk(data_root,j,Wsz,T);
                    fprintf(' loading chunk %d', j);
                end

                r = topmargin+Jsz/2+buff+ceil((Wsz-Jsz-(topmargin+2*buff))*rand);
                c = Jsz/2+buff+ceil((Wsz-Jsz-2*buff)*rand);
                Y = reshape(F(r-Jsz/2:r+Jsz/2-1,c-Jsz/2:c+Jsz/2-1,:),J,T);
                Z = reshape(G(r-Jsz/2:r+Jsz/2-1,c-Jsz/2:c+Jsz/2-1,:,:),J,2,T-1);

                Y = Y - mean(Y(:));

                Xr(:,:,b) = Y;
                Zr(:,:,:,b) = Z;

                fprintf('\r%d / %d', b, B);
            end
            fprintf('\n');

            Xr = reshape(Xr, J, T*B) / std(Xr(:));

            % PCA
            C = Xr*Xr' / (B*T);
            [V,D] = eig(C);
            [val,idx] = sort(diag(D), 'descend');
            V = V(:,idx);
            D = diag(val);

            %% save to cache
            save(filename_data, 'Xr', 'V', 'D', 'Zr', '-v7.3');
        end

        %X = V(:,1:L)'*Xr;
        %X = diag(sqrt(1./diag(D(1:L,1:L)))) * X;

        Xr = reshape(Xr, J, T, B);

end


