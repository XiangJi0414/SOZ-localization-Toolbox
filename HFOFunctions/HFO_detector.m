% This is the main function to detect HFOs
function results=HFO_detector(signal,parameters)
% signal (N_channels,data_points): signal containing potential HFOs
% parameters of the detector
signal=single(signal); % converting to single for faster operation (necessary for gpu)
fs=parameters.fs; % sampling frequency
L=size(signal,2); % length of the signal
time=0:L-1;
time=time/fs; % creating time vector based on the length of the signal and fs
time=single(time);
%%
temp=mod(L,parameters.chunks); % making sure the user knows some data might be discarded.
% I can run the HFO detector for what remains but that does not make sense!
if temp~=0
    warning('%.0f data points are discarded because mod(L,chunk)~=0',temp);
end
l_chuncks=floor(L/parameters.chunks); % length of each data chunk
%
first=true;
for channel=parameters.channels % looping through the channels
    for chunk=1:parameters.chunks % looping through data chunks
        current_inds= (chunk-1)*l_chuncks+1:chunk*l_chuncks; % picking up the indices of the current chunk
        if parameters.verbose
            fprintf('Channel %.0f, chunk %.0f\n',channel,chunk); % displaying where we are if verbose is true.
        end
        %
        data=signal(channel,current_inds); % picking up the current data to process
        current_time=time(current_inds); % picking up current time
        if parameters.gpu==true
            data=gpuArray(data); % creating gpu array if chosen by the user
        end
        if first==true % is this the first time we are going through the data? We only need to create the filter bank once!
            first=false; % not anymore!
            % creating the filter bank for cwt analysis
            if strcmp(parameters.wavelet,'morse')
                fb = cwtfilterbank(SignalLength=length(data),SamplingFrequency=fs,Wavelet=parameters.wavelet,WaveletParameters=parameters.WaveletParameters,VoicesPerOctave=parameters.VoicesPerOctave,FrequencyLimits=parameters.FrequencyLimits);
            else
                fb = cwtfilterbank(SignalLength=length(data),SamplingFrequency=fs,Wavelet=parameters.wavelet,VoicesPerOctave=parameters.VoicesPerOctave,FrequencyLimits=parameters.FrequencyLimits);
            end
        end
        [B,f] = cwt(data,FilterBank=fb); % performing cwt analysis
        if parameters.gpu==true % sending the data back to cpu if we are doing gpu calculation
            B=gather(B);
            f=gather(f);
        end
        z=abs(B); % finding the magnitude.
        clear B; % clean up some big array
        f=f';
        %
        mean_magnitude=mean(z,2)'; % calculating mean of magnitude for the whole signal.
        if strcmp(parameters.compare,'power')
            z=z.^2; % getting read of z itself if it is chosen to use power for comparison
            mean_power=mean(z,2)'; % calculating mean of power for the whole signal.
            mean_z=mean_power;
        else
            temp=z.^2;
            mean_power=mean(temp,2)';
            clear temp;
            mean_z=mean_magnitude;
        end
        results.channel(channel,chunk).mean_magnitude=mean_magnitude; % saving mean magnitude in the resutls
        results.channel(channel,chunk).mean_power=mean_power; % saving mean power in the results.
        %
        if ~isempty(parameters.ripple) % did the user ask to find ripples?
            rng(1); % setting rng for reproducibility
            cutoff=(f>=parameters.ripple.range(1)).*(f<=parameters.ripple.range(2));cutoff=logical(cutoff); % finding the cut off frequency
            z_th=z(cutoff,:); % cutting off data
            r=randi(size(z_th,1)*size(z_th,2),floor(size(z,2)/5),1);% picking up 5% of the data for ecdf creation.
            [yy,xx]=ecdf(double(z_th(r))); % constructing ecdf. I convert the data back to double here because I have seen issues if I use signle here.
            [~,im]=min(abs(yy-parameters.ripple.ecdf)); % find where we are crossing threshold
            th=xx(im); % picking up threshold
            results.channel(channel,chunk).Prob.ripple.th=th; % saving the threshold value
            if parameters.ripple.save_ecdf==true
                results.channel(channel,chunk).Prob.ripple.ecdf_x=single(xx); % saving ecdf results
                results.channel(channel,chunk).Prob.ripple.ecdf_y=single(yy); % saving ecdf results
            end
            z_th(z_th<th)=0; % setting z (magnitude or power based on the user's choice) 0 for values leff than threshold
            f_cut=f(cutoff); % frequency range for the cutoff
            CC = bwconncomp(z_th,4); % finding connected blobs.
            % picking up blobs that pass criteria set with parameters
            [events_R,artifacts_R.frequency_range,artifacts_R.amplitude]=clean_results(CC,f_cut,z_th,f,mean_z,current_time,parameters.compare,parameters.center,parameters.ripple,parameters.unimodal);
            results.channel(channel,chunk).ripples=events_R; % saving resutls
            results.channel(channel,chunk).artifacts.ripple=artifacts_R; % saving events that are potentialy artifacts
        end
        %%
        if ~isempty(parameters.fast_ripple) % did the user ask to find fast-ripples?
            % comments are all the same as ripples. The parameters are the
            % only things that are channging to find fast-ripples.
            rng(1);
            cutoff=(f>=parameters.fast_ripple.range(1)).*(f<=parameters.fast_ripple.range(2));
            cutoff=logical(cutoff);
            z_th=z(cutoff,:);
            r=randi(size(z_th,1)*size(z_th,2),floor(size(z,2)/5),1);
            [yy,xx]=ecdf(double(z_th(r)));
            [~,im]=min(abs(yy-parameters.fast_ripple.ecdf));
            th=xx(im);
            results.channel(channel,chunk).Prob.fast_ripple.th=th;
            if parameters.fast_ripple.save_ecdf==true
                results.channel(channel,chunk).Prob.fast_ripple.ecdf_x=single(xx);
                results.channel(channel,chunk).Prob.fast_ripple.ecdf_y=single(yy);
            end
            z_th=z(cutoff,:);
            z_th(z_th<th)=0;
            f_cut=f(cutoff);
            CC = bwconncomp(z_th,4);
            [events_FR,artifacts_FR.frequency_range,artifacts_FR.amplitude]=clean_results(CC,f_cut,z_th,f,mean_z,current_time,parameters.compare,parameters.center,parameters.fast_ripple,parameters.unimodal);
            results.channel(channel,chunk).fast_ripples=events_FR;
            results.channel(channel,chunk).artifacts.fast_ripple=artifacts_FR;
        end
        %%
        if ~isempty(parameters.spike) % did the user ask to spikes?
            [events_S,artifacts_S]=find_spikes_basic(z,f,current_time,parameters.spike.range,parameters.spike.prom,parameters.spike.artifact);
            results.channel(channel,chunk).spikes=events_S;
            results.channel(channel,chunk).artifacts.spike=artifacts_S;
        end
    end
    %%
    results.frequency=f; % saving the frequncies of the filter bank
    results.parameters=parameters; % saving the parameters in the results to make sure that we know what parameters were used.
end
