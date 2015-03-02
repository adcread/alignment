% perform alignment procedure multiple times to investigate effect on
% constellations at receivers

superman = 128;

constellationPri{1} = zeros(3,superman);
constellationPri{2} = zeros(2,superman);
dataTxPri{1} = zeros(3,superman);
dataTxPri{2} = zeros(2,superman);

dataRxPri{1} = zeros(3,superman);
dataRxPri{2} = zeros(2,superman);


constellationPub{1} = zeros(3,superman);
constellationPub{2} = zeros(2,superman);

dataTxPub{1} = zeros(3,superman);
dataTxPub{2} = zeros(2,superman);

dataRxPub{1} = zeros(3,superman);
dataRxPub{2} = zeros(2,superman);

for batman = 1:superman
    
    alignment;
    
    if ~isempty(decodedPri{1})
        constellationPri{1}(:,batman) = equalisedPri{1}.';
        dataTxPri{1}(:,batman) = privateCodeword(1,1:3).';
        dataRxPri{1}(:,batman) = dataPri{1}.';
    end
    
    if ~isempty(decodedPub{1})
        constellationPub{1}(:,batman) = equalisedPub{1}.';
        dataTxPub{1}(:,batman) = publicCodeword(1,1:3).';
        dataRxPub{1}(:,batman) = dataPub{1}.';
    end
    
    if ~isempty(decodedPri{2})
        constellationPri{2}(:,batman) = equalisedPri{2}.';
        dataTxPri{2}(:,batman) = privateCodeword(2,1:2).';
        dataRxPi{2}(:,batman) = dataPri{2}.';
    end
    
    if ~isempty(decodedPub{2})
        constellationPub{2}(:,batman) = equalisedPub{2}.';
        dataTxPri{2}(:,batman) = publicCodeword(2,1:2).';
        dataRxPri{2}(:,batman) = dataPri{2}.';
    end
    
end

%% Create plots of the constellations found

subplot(4,3,1);
    scatter(real(constellationPri{1}(1,:)),imag(constellationPri{1}(1,:)),'.');
    title('User 1 Private Stream 1');
subplot(4,3,2);
    scatter(real(constellationPri{1}(2,:)),imag(constellationPri{1}(2,:)),'.');
    title('User 1 Private Stream 2');  
subplot(4,3,3);
    scatter(real(constellationPri{1}(3,:)),imag(constellationPri{1}(3,:)),'.');
    title('User 1 Private Stream 3');
    
subplot(4,3,4);
    scatter(real(constellationPri{2}(1,:)),imag(constellationPri{2}(1,:)),'.');
    title('User 2 Private Stream 1');
subplot(4,3,5);
    scatter(real(constellationPri{2}(2,:)),imag(constellationPri{2}(2,:)),'.');
    title('User 2 Private Stream 2');
    
subplot(4,3,7);
    scatter(real(constellationPub{1}(1,:)),imag(constellationPub{1}(1,:)),'.');
    title('User 1 Public Stream 1');
subplot(4,3,8);
    scatter(real(constellationPub{1}(2,:)),imag(constellationPub{1}(2,:)),'.');
    title('User 1 Public Stream 2');    
subplot(4,3,9);
    scatter(real(constellationPub{1}(3,:)),imag(constellationPub{1}(3,:)),'.');
    title('User 1 Public Stream 3');
    
subplot(4,3,10);
    scatter(real(constellationPub{2}(1,:)),imag(constellationPub{2}(1,:)),'.');
    title('User 2 Public Stream 1');
subplot(4,3,11);
    scatter(real(constellationPub{2}(2,:)),imag(constellationPub{2}(2,:)),'.');
    title('User 2 Public Stream 2');  


%% Send an email to IFTTT to notify Mobile user.
    
matlabmail('trigger@recipe.ifttt.com','#matlab','#matlab','chris@celeritas.org.uk','aUgme2nt');


    