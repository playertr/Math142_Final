function plot_ellipsoid(Cov, mu)
    % plot 3D ellipsoid
    % developed from the original demo by Rajiv Singh
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/42966
    % 5 Dec, 2002 13:44:34
    % Example data (Cov=covariance,mu=mean) is included.

    [U,L] = eig(Cov);
    % L: eigenvalue diagonal matrix
    % U: eigen vector matrix, each column is an eigenvector

    % For N standard deviations spread of data, the radii of the eliipsoid will
    % be given by N*SQRT(eigenvalues).

    N = 1; % choose your own N
    radii = N*sqrt(diag(L));

    % generate data for "unrotated" ellipsoid
    [xc,yc,zc] = ellipsoid(0,0,0,radii(1),radii(2),radii(3));

    % rotate data with orientation matrix U and center mu
    a = kron(U(:,1),xc); 
    b = kron(U(:,2),yc); 
    c = kron(U(:,3),zc);

    data = a+b+c; n = size(data,2);

    x = data(1:n,:)+mu(1); 
    y = data(n+1:2*n,:)+mu(2); 
    z = data(2*n+1:end,:)+mu(3);

    % now plot the rotated ellipse
    % sc = surf(x,y,z); shading interp; colormap copper
    surfl(x, y, z); colormap copper

end