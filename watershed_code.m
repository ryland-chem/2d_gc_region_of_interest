%calculate the gradient, find regions (edges) with extreme changes
x = pValCutOff;
gmag = imgradient(x);

%need to mark objects in the foreground
se = strel('disk', 1);
Io = imopen(x, se);

%opening by reconstruction
Ie = imerode(x, se);
Iobr = imreconstruct(Ie, x);

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
clims = [0 1];
imagesc((flipud(pValCutOff)), clims);
colormap jet;
set(gca,'YDir','normal');
ylabel("2nd Dimension Acquisitions"); xlabel("1st Dimension Acquisitions");
hold on
himage = imagesc(flipud(Lrgb));
himage.AlphaData = 0.5;





