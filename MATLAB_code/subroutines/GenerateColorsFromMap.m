function  [colors,selectedcolormap] = GenerateColorsFromMap(NumberOfColors,Map)
%This function loads the colormap (also given as output) as specified by the
%string "Map". This colormap (which is effectively a 64x3 array) is then used
%to generate a NumberOfColorsx3 array, which can be used to color figures. 
%A single color for the generated sequence is colors(a,:).

%The available colormaps are:
% *RainbowReinier [Default]: Blue-Cyan-Green-Yellow-Red. Often found in MSc
%                            Thesis R.M.J. Janssen.
% *RBG_Miranda: Red-Black-Green. Created to make a chemical heat map for
%                            M.G.M. Kok
% *RWB_Miranda: Red-White-Blue. Created to make a chemical heat map for
%                            M.G.M. Kok
% *PGO_Miranda: Purple-Green-Orange. Created to make a chemical heat map for
%                            M.G.M. Kok
% *cool,hot: As matlab default colormaps.

%This routine requires the 'colormapstorage.mat' file to be present.

%Load the colormap
lastwarn('') %reset last warning
load('colormapstorage',Map); %attempt to load the map
%Catch unknown colormaps
if strcmp(lastwarn,['Variable ',Map,' not found.'])
    %Attempted to load unknown colormap. Switch to default.
    fprintf(['Warning GenerateColorsFromMap: No colormap found matching the specified name. ', ...
        'The default colormap is chosen. \n'])
    load('colormapstorage','RainbowReinier')
end
selectedcolormap = eval(genvarname(Map));

%Use the colormap and interpolation to generate the requested list of
%colors.
ML = size(selectedcolormap,1); %length of the colormap.
if NumberOfColors == 1
    %Problem with regular method due to division by NumberOfColors-1
    map_locations = ML/2;
else
    %Regular operation
    map_locations = 1:(ML-1)/(NumberOfColors-1):ML;
end

colors = zeros(NumberOfColors,3);
for p=1:NumberOfColors
    for q=1:3 %Loop over 3 RGB values
        colors(p,q) = interp1(1:ML,selectedcolormap(:,q),map_locations(p),'linear');
    end
end

%==========================================================================
%==========================================================================
% HOW TO ADD YOUR OWN COLORMAP
%==========================================================================
% You can create your own colormap and add it to this function by the
% following procedure. Any line starting with %> is actual code you should
% type in matlab. Change "yourmapname" to the desired colormap name.

% Create your own colormap using the matlab colormap editor
%> colormapeditor
% Mess around until you like the colormap. Then press OK.
%> yourmapname = colormap;
% Save it inside the 'colormapstorage'
%> save('colormapstorage','yourmapname','-append')
% You can now use your own colormap with this program
% colors = GenerateColorsFromMap(NumberOfColors,'yourmapname')

%==========================================================================
%==========================================================================