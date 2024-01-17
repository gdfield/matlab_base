function fitmovie = concat_movie(fitmovie_cell)
% NB 2016-06-03
% takes in movies in cells and puts them back to back.
n_blocks = length(fitmovie_cell);
fitmovie = [];
if isfield(fitmovie_cell{1}, 'matrix')
    for block = 1:n_blocks
        fitmovie = cat(3, fitmovie, fitmovie_cell{block}.matrix);
    end
else
    for block = 1:n_blocks
        fitmovie = cat(3, fitmovie, fitmovie_cell{block});
    end
end
end