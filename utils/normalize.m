function [nX] = normalize(X)

nX = (X-min(X(:)))./(max(X(:))-min(X(:)));