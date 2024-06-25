%% PARAMETRIC ANALYSIS OF A UNIT CELL WITH A SPHERICAL INCLUSION
clc, clear all

% Parameters
L = [1 1] ; dx = 1/20 ; % unit cell size & element length (mm)
D = linspace(0,.99,25) ; .6 ; % inclusion diameter
p = 2 ; logspace(-1,1,25) ; % inclusion power function (2==sphere,inf==cube,1==diamond)
xc = L.*cat(3,[.5 .5]) ; % position of inclusion(s)
Ep = [4e3 2e3] ; % phase young modulus (MPa)
nuP = [.3 .3] ; % phase Poisson ratio
rhoP = [1220e-12 1215e-12] ; % phase material densities (tons/mm3)
k = [1 0].*pi/norm(L(1)).*linspace(0.01,1,100)' ; % wavevectors (rad/mm) [nK nCoord]
nModes = 20 ; % extract the first nModes
redux = true ; % apply basis reduction ?

% Parametric analysis
IMG = {} ;
for par = 1:numel(D+p)
% UNIT CELL DEFINITION
    unitcell = [] ;
    unitcell.mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
    unitcell.mesh.setElementTypes(pkg.geometry.mesh.elements.LagrangeElement('quad',2)) ; % quadratic elements
% Two-phase geometry (defined by a levelset function at nodes)
% lvlst<0 for phase 1 and lvlst>0 for phase 2 (lvlst==0 is the interface)
%     rxP = min((sum(abs(unitcell.mesh.Nodes-xc).^p(par),2)).^(1/p(par)),[],3) ;
%     unitcell.lvlst = ((D/2)-rxP) ; % phase 1 outside
    rxP = min((sum(abs(unitcell.mesh.Nodes-xc).^p,2)).^(1/p),[],3) ;
    unitcell.lvlst = ((D(par)/2)-rxP) ; % phase 1 outside
% Material phase properties
    unitcell.E_p = Ep ; % young moduli of phases (MPa)
    unitcell.nu_p = nuP ; % poisson ratio of phases
    unitcell.rho_p = rhoP ; % density of phases (tons/mm3)
    unitcell.G_p = unitcell.E_p./2./(1+unitcell.nu_p) ; 
    unitcell.C_p = pkg.fem.bloch.stiffness(unitcell.E_p,unitcell.G_p) ;
    unitcell.S_p = pkg.math.inv(unitcell.C_p) ;
% Distribute on the mesh
    [ee,we,ie] = unitcell.mesh.integration() ; 
    Ne = unitcell.mesh.interpMat(ee,ie) ; % interpolation at quadrature points
    unitcell.phase = 1+double(Ne*unitcell.lvlst>0) ;
    unitcell.C = unitcell.C_p(:,:,unitcell.phase) ;
    unitcell.rho = unitcell.rho_p(unitcell.phase) ;
% Estimate mean properties
    unitcell.volFrac = full(sum(sparse((1:numel(ie))',unitcell.phase,we,numel(ie),numel(unitcell.E_p)),1))/sum(we(:)) ;
    unitcell.C_voigt = sum(unitcell.C_p.*reshape(unitcell.volFrac,1,1,[]),3) ;
    unitcell.C_reuss = inv(sum(unitcell.S_p.*reshape(unitcell.volFrac,1,1,[]),3)) ;
    unitcell.rho_mean = sum(unitcell.rho_p.*unitcell.volFrac,2) ;
% Estimate isotropic phase velocities
    unitcell.cl_voigt = sqrt(unitcell.C_voigt(1,1)/unitcell.rho_mean) ;
    unitcell.ct_voigt = sqrt(unitcell.C_voigt(end,end)/unitcell.rho_mean) ;
    unitcell.cl_reuss = sqrt(unitcell.C_reuss(1,1)/unitcell.rho_mean) ;
    unitcell.ct_reuss = sqrt(unitcell.C_reuss(end,end)/unitcell.rho_mean) ;
    
% WIDE BAND WAVE PROPAGATION IN THE 3D SOLID
% Reduction parameters
    [k,w,U,Q,lq] = dispersion(unitcell,k,NaN,nModes,[],redux) ;
    
% DISPLAY 
clf reset ; 
% dispersion branches
    knorm = sqrt(sum(abs(k).^2,2)) ;
    clrs = get(gca,'colororder') ;
    pl = plot3(repelem(knorm,nModes),real(w(:)),imag(w(:)),'.k','displayname','SAFE') ;
    xlabel 'Wavenumber $k$ (rad/mm)'
    ylabel 'Frequency $\omega$ (rad/s)'
    title("\boldmath $D = "+num2str(D(par),2)+"$") ;
%     title("\boldmath $p = "+num2str(p(par),2)+"$") ;
% Rescaling
    ct = sqrt(unitcell.ct_voigt.*unitcell.ct_reuss) ; % average shear velocity
    fc = .5*ct/L(1) ; % cutoff frequency
    kc = pi/L(1) ; % cutoff wavenumber
    pl.XData = pl.XData./kc ;
    pl.YData = pl.YData./2./pi./fc ;
    xlabel 'Normalized Wavenumber $\frac{k \times L}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2 L}{c_t}$'
    set(gca,'ylim',[0 4])
% the unitcell
    axes('position',[0.17412      0.24167      0.13294      0.13452]) ; 
    axis equal tight off ; 
    plot(unitcell.mesh,'VisibleEdges','none','FaceAlpha',0) ; 
    cut = unitcell.mesh.simplex.cut(unitcell.lvlst) ; % the cut is performed on the simplex mesh for speed
    plot(cut.IN,'VisibleEdges','none') ; 
    light ;
% Save Image
    IMG{end+1} = getframe(gcf) ;
end

%% EXPORT THE IMAGES AS GIF
pkg.export.gif(IMG) ;
