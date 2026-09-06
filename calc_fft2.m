function [freq_vec,shifted_fft_vec,N,units]=calc_fft2(signal,fs,varargin)

    posFreqFlag=1;
    padFlag=0;
    dimFlag=1;
    unitFlag='amp';

    if nargin>2
        for karg=1:2:nargin-2
            switch lower(varargin{karg})
                case 'posfreqflag'
                    posFreqFlag=varargin{karg+1};
                case 'padflag'
                    padFlag=varargin{karg+1};
                case 'dimflag'
                    dimFlag=varargin{karg+1};
                case 'unitflag'
                    unitFlag=lower(varargin{karg+1});
                otherwise
                    error('Input argument not supported yet!')
            end
        end
    end
    

    %% pad signal if needed
    if padFlag
       signal=padarray(signal,padFlag*length(signal),0,'both'); 
    end
    

    %% calculate the physical (�real�) frequencies vector:
    N = 2^nextpow2(size(signal,dimFlag));
%     N=2^nextpow2(length(signal));	% number of time-domain samples
    dt=1/fs;		% time between samples = (sampling frequency)^-1
    df=1/(N*dt);      % the frequency resolution (df=1/max_T)

    
    %% generate proper physical frequency vector (after fftshift)
    if mod(N,2)==0
        freq_vec= df*((1:N)-1-N/2);     % frequency vector for EVEN length vectors: f =[-f_max,-f_max+df,...,0,...,f_max-df]
    else
        freq_vec= df*((1:N)-0.5-N/2);   % frequency vector for ODD length vectors f =[-f_max,-f_max+fw,...,0,...,f_max]
    end

    shifted_fft_vec=fftshift(fft(signal,N,dimFlag),dimFlag); % Fourier-transform the signal

    

    if posFreqFlag % Delete negative frequencies
        % Crop the negative frequencies
        freq_vec(1:N/2) = [];
        
        % Make a copy of the spectra matrix
        s_tmp = shifted_fft_vec;
        
        % Get the initial size and number of dimensions
        sizeInit = size(s_tmp);
        
        if dimFlag ==1 % to do: debug this section
            % Delete the first half of the first dimension of the spectra matrix
            s_tmp(1:N/2,:) = []; % Spectra matrix is now 2D

            % Reshape according to the new number of elements
            s_tmp = reshape(s_tmp,[size(s_tmp,1), sizeInit(2:end)]);
        else
            sdims = ndims(s_tmp);
            dims = 1:sdims;

            % Delete the dimension corresponding to dimFlag
            dims(dims==dimFlag) = [];

            % Permute the spectra matrix: dimFlag first
            s_tmp = permute(s_tmp,[dimFlag, dims]);

            % Compare the number of permuted dimensions with the initial number of dimensions
            % One dimension can drop if a singleton dimension is last after permutation
            pdims = ndims(s_tmp);
            if pdims ~= sdims
                sdims = pdims;
            end

            % Delete the first half of the first dimension of the spectra matrix
            s_tmp(1:N/2,:) = []; % Spectra matrix is now 2D

            % Reshape according to the new number of elements
            s_tmp = reshape(s_tmp,[size(s_tmp,1), sizeInit(2:end)]);

            % Permute back to the original dimensions arrangement
            s_tmp = permute(s_tmp,1:sdims);
        end
        
        % Save result
        shifted_fft_vec = s_tmp;
        shifted_fft_vec(2:end) = 2*shifted_fft_vec(2:end);
        

    end
    
    
    switch unitFlag
        case 'amp'
            shifted_fft_vec = shifted_fft_vec/N*sqrt(2);
            units = 'V';
        case 'raw'
            % Nothing to do
            units = 'raw';
        case 'rms'
            shifted_fft_vec = shifted_fft_vec/N;
            units = 'Vrms';
        case 'pow'
            shifted_fft_vec = abs(shifted_fft_vec/N).^2;
            units = 'Vrms^2';
        otherwise
            warning(['Unknown unitFlag string' unitFlag '. Defaulting to ''amp''.'])
            shifted_fft_vec = shifted_fft_vec/N*sqrt(2);
            units = 'V';
    end
            
    

end