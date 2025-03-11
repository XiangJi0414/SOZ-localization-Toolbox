% This function compares detected HFOs with ground truth. 
function [TP,FP,FN,inds_TP]=compare_resutls(detected,ground_truth,delta)
% detected: time of detected HFO
% ground_truth: time of ground truth event
% delta: length of time window to test detected events against ground truth
detected=sort(detected);
ground_truth=sort(ground_truth);
TP=0;
FN=0;
inds_TP=false(length(detected),1);
detected_og=detected;
% Looping through ground truth data
for i=1:length(ground_truth)
    % Creating a window with length of delta around each ground truth data
    current_time=ground_truth(i);
    CI_P=current_time+delta/2;
    CI_N=current_time-delta/2;

    inds=(detected>CI_N).*(detected<CI_P);
    inds=logical(inds);

    % Checking if any events were detected
    inds_og=(detected_og>CI_N).*(detected_og<CI_P);
    inds_og=logical(inds_og);

    if sum(inds)==0
        FN=FN+1; % Increase false negative by 1 if nothing was detected        
    else
        TP=TP+1; % Increase true postive by 1 if an event was detected
        detected(inds)=[]; % remove events that were counted as true positives
        inds_TP(inds_og)=true;
    end
end
FP=length(detected); % The length of detected events is now false postives beachse we looped through all the ground truth events and removed all true positives.