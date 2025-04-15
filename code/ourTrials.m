function [NewData] = ourTrials(Data)

% Définition des valeurs à filtrer
val_col3 = ([20, 10, 0, -10, -20] / 180) * pi;
val_col456 = [0, 0.2, 0.4, 0.6];

% Création du masque de filtrage
mask = (Data(:,2) == 0) & ...
       ismember(Data(:,3), val_col3) & ...
       (ismember(Data(:,4), val_col456) | ismember(Data(:,5), val_col456) | ismember(Data(:,6), val_col456)) & ...
       ~( (Data(:,4) == 1) | (Data(:,5) == 1) | (Data(:,6) == 1) ); % Exclure les lignes contenant 1

% Extraction des lignes correspondant aux conditions
NewData = Data(mask, :);
NewData = NewData(:,1:12);

end