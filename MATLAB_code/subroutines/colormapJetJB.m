function CM = colormapJetJB(npts)
% creates colormap jet with the yellow bnot too bright.
CM = colormap(jet(npts));
CM(:,2)=CM(:,2)*0.8;

end