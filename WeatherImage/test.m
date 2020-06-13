nasa = wmsfind('NASA','SearchField','serverurl');
layer = nasa.refine('bluemarbleng', ...
    'SearchFields','layername', ...
    'MatchType', 'exact');

[A,R] = wmsread(layer(1));

axesm globe
axis off
geoshow(A,R)
title('Blue Marble')