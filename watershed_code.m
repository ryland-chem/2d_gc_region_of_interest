%calculate the gradient, find regions (edges) with extreme changes
gmag = imgradient(pValCutOff);

%need to mark objects in the foreground
se = strel('disk', 4000);
Io = imopen(pValCutOff, se);

%opening by reconstruction
Ie = imerode(pValCutOff, se);
Iobr = imreconstruct(Ie, pValCutOff);

%followed by closing, helps to fill in regions that look like dougnuts
Ioc = imclose(Io, se);

Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);

%regional maxima calc
%this affects the watershed here
fgm = imregionalmax(Iobrcbr);

se2 = strel(ones(5,5));
fgm2 = imclose(fgm,se2);
fgm3 = imerode(fgm2,se2);

%again affects the watershed
fgm4 = bwareaopen(fgm3, 5);

%compute background markers
bw = imbinarize(Iobrcbr);

%watershed
D = bwdist(bw);
%DL has our different regions of interest in it
DL = watershed(D);

%overlay probability cutoff with the watershed
Lrgb = label2rgb(DL, 'jet','w','shuffle');
image(pValCutOff, 'CDataMapping', 'scaled');
hold on
himage = imshow(Lrgb);
himage.AlphaData = 0.3;





