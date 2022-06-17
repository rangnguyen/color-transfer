% Test color_transfer.m
pd = '..\data\example10\';
Is = im2double(imread([pd 'source.jpg']));
It = im2double(imread([pd 'target.jpg']));
Io = color_transfer(Is, It);
imwrite(Io, [pd 'color_transfer.jpg']);

% Compute the evaluation metric
Io = im2double(imread([pd 'color_transfer.jpg']));
d = evaluate_metric(It, Io);
disp(['Gamut distance: ' num2str(d)]);