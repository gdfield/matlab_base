function NewElectrodes = NS512_519Inverse(Electrodes)

a = get512_519_inverseElecs(); 

% NewElectrodes = a(Electrodes);
NewElectrodes = zeros(1,length(Electrodes)); 
for ii=1:length(Electrodes)
    NewElectrodes(ii)=find(a==Electrodes(ii));
end
