%% SIMULATION ASSOCIATED TO THE WIDEBAND CHARACTERIZATION OF AN ARCHITECTURED MATERIAL
clc, clear all

%% UNIT CELL DEFINITION
% Define the unit cell as a structure containing at least the fields:
%   - "mesh" as the 3D unit cell mesh (must be periodic)
%   - the material properties ("C" & "rho") at the mesh quadrature points

unitcell = [] ;
% Unit cell geometry & discretization
    L = [1 1] ; dx = 1/20 ; % unit cell size & element length (mm)
    unitcell.mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
% Use quadratic elements improves a lot the convergence
    unitcell.mesh.setElementTypes(pkg.geometry.mesh.elements.LagrangeElement('quad',2)) ;
% Two-phase geometry (defined by a levelset function at nodes)
% lvlst<0 for phase 1 and lvlst>0 for phase 2 (lvlst==0 is the interface)
    switch 1
        case 0 % uniform phase=phase1
            unitcell.lvlst = -ones(unitcell.mesh.nNodes,1) ;
        case 1 % spherical inclusion(s)
            p = 4/2 ; % power (inf for cubes, 2 for spheres, 1 for diamonds
            D = 15/20 ; % radius 
            xc = L.*cat(3,[.5 .5] ...
                       ...,[0 0] ,[0 1] ,[1 1] ,[1 0] ...
                        ) ;
            rxP = min(sum(abs(unitcell.mesh.Nodes-xc).^p,2),[],3) ;
            unitcell.lvlst = (((D/2).^p)-rxP) ; % phase 1 outside
        case 2 % 2-phases gyroid
            x = 2*pi.*unitcell.mesh.Nodes./L ;
            unitcell.lvlst = sin(2*x(:,1)).*cos(3*x(:,2)) + sin(3*x(:,2)) ;
    end
    cut = unitcell.mesh.simplex.cut(unitcell.lvlst) ; % the cut is performed on the simplex mesh for speed
% Display
    figure(1)
    clf reset ; axis equal tight ; 
    plot(unitcell.mesh,'VisibleEdges','none','FaceAlpha',0) ; 
    plot(cut.IN,'VisibleEdges','none') ; 
    light ;
% Material phase properties
    unitcell.E_p = [4e3 2e3] ; % young moduli of phases (MPa)
    unitcell.nu_p = [.3 .45] ; % poisson ratio of phases
    unitcell.rho_p = [1220e-12 1215e-12] ; % density of phases (tons/mm3)
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
% Display
    unitcell
    
%% WIDE BAND WAVE PROPAGATION IN THE 3D SOLID
    k = [1 0].*pi/norm(L(1)).*linspace(0.01,1,100)' ; % wavevectors (rad/mm) [nK nCoord]
    nModes = 20 ; % extract the first nModes
% Reduction parameters
    redux = true ; uint8(20) ; sqrt(eps) ; % apply basis reduction ?
    [k,w,U,Q,lq] = dispersion(unitcell,k,NaN,nModes,[],redux) ;
% Display
    clf(figure(2)) ;
    knorm = sqrt(sum(abs(k).^2,2)) ;
    clrs = get(gca,'colororder') ;
    pl = plot3(repelem(knorm,nModes),real(w(:)),imag(w(:)),'.k','displayname','SAFE') ;
    xlabel 'Wavenumber $k$ (rad/mm)'
    ylabel 'Frequency $\omega$ (rad/s)'
    title('\bf Dispersion curves') ;
% Rescaling
    ct = sqrt(unitcell.ct_voigt.*unitcell.ct_reuss) ; % average shear velocity
    fc = .5*ct/L(1) ; % cutoff frequency
    kc = pi/L(1) ; % cutoff wavenumber
    pl.XData = pl.XData./kc ;
    pl.YData = pl.YData./2./pi./fc ;
    xlabel 'Normalized Wavenumber $\frac{k \times L}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2 L}{c_t}$'
    
%% Animate modes
    Kdir = repelem(k./sqrt(sum(k.^2,2)),nModes,1) ;
    extrude = 3 ;
    gifOnClick = true ;
    plw = pkg.fem.bloch.waveModeAnimation(unitcell.mesh,Kdir,U,pl,extrude,gifOnClick) ;
    plw.AlphaData = -repmat(unitcell.phase(:),extrude^2,1) ;
    set(gca,'alim',get(gca,'alim')+[-1 1]*1) ;
    
%% Animate the reduced basis
    clf reset  ; pl = plot(lq,'o') ; set(gca,'yscale','log')
    plw = pkg.fem.bloch.waveModeAnimation(unitcell.mesh,zeros(size(Q,2),2),reshape(Q,unitcell.mesh.nNodes,2,[]),pl,extrude,gifOnClick) ;
    plw.AlphaData = -repmat(unitcell.phase(:),extrude^2,1) ;
    set(gca,'alim',get(gca,'alim')+[-1 1]*1) ;
    
%% TRY OTHER REDUCTION BASES
    nModesR = 100
    %[~,w0,U0] = dispersion(unitcell,[0 0],NaN,nModesR,[],false) ; % cutoff frequencies
    [k0,~,U0] = dispersion(unitcell,[1 0],0,nModesR,[],false) ; % cutoff wavenumbers
    Q = reshape(U0,prod(size(U0,[1 2])),[]) ;
    
%% EFFECT OF THE REDUCED BASIS
    clf ;
    plot(repelem(knorm,nModes)/kc,real(w(:))/2/pi/fc,'.k','markersize',10) ;
    pl = plot(NaN,NaN,'or','markersize',3,'linewidth',1) ;
    xlabel 'Normalized Wavenumber $\frac{k \times L}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2 L}{c_t}$'
    set(gca,'ylimmode','manual','xlimmode','manual') ;
    IMG = {} ;
    for rr = 1:size(Q,2)
        [~,wr] = dispersion(unitcell,k,NaN,min(nModes,rr),[],Q(:,1:rr)) ;
        pl.XData = repelem(knorm,min(nModes,rr))/kc ;
        pl.YData = real(wr(:))/2/pi/fc ;
        title("\boldmath $R="+string(rr)+"$")
        drawnow ;
        IMG{end+1} = getframe(gcf) ;
    end
%% VISUALIZE THE REDUCED BASIS
    
    Pn = Q./sqrt(sum(abs(Q).^2,1)) ;
    Un = U./sqrt(sum(abs(U).^2,[1 2])) ;
    Un = reshape(Un,size(Q,1),[]) ;
    an = Pn\Un ; 
    
    clr = (abs(an(1,:))).' ;
    clrs = linspecer(5)*1 ; get(gca,'colororder') ; jet(4) ;
    clr = 1-((abs(an(1:size(clrs,1),:))).'*(1-clrs)) ;
    
    clf ;
    scatter3(repelem(knorm,nModes)/kc,real(w(:))/2/pi/fc,imag(w(:))/2/pi/fc,15,clr,'filled') ;
    xlabel 'Normalized Wavenumber $\frac{k \times L}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2 L}{c_t}$'
    %%
    
    
    
    sig = sqrt(sum(abs(an).^2,2)) ;
    [sig,is] = sort(sig,'descend') ;
    clf ; pl = plot(sig,'o') ;
    plw = pkg.fem.bloch.waveModeAnimation(unitcell.mesh,zeros(size(Pn,2),2),reshape(Pn,unitcell.mesh.nNodes,2,[]),pl,extrude,gifOnClick) ;
    plw.AlphaData = -repmat(unitcell.phase(:),extrude^2,1) ;
    set(gca,'alim',get(gca,'alim')+[-1 1]*1) ;
    
%% PRINCIPAL COMPONENT ANALYSIS
% Basis of modes to reduce
    U_pca = reshape(U,prod(size(U,[1 2])),[]) ;
% Normalization
    [K0,K0i,Kij,M] = pkg.fem.bloch.FEM(unitcell.mesh,unitcell.C,unitcell.rho,[],[]) ;
    M = M(1:end*2/3,1:end*2/3) ;
    K0 = K0(1:end*2/3,1:end*2/3) ;
    K1 = K0i{1}(1:end*2/3,1:end*2/3) ;
    K2 = Kij{1,1}(1:end*2/3,1:end*2/3) ;
    kk = repmat(k(:,1),[nModes 1]) ;
    Ku = U_pca'*((K0*U_pca) + (K1*U_pca).*abs(kk(:).') + (K2*U_pca).*abs(kk(:).').^2) ;
    w2Mu = U_pca'*((M*U_pca).*abs(w(:).').^2) ;
    Ku = Ku+Ku' ;
    w2Mu = w2Mu+w2Mu' ;
%% The PCA
    [Qm,s2] = eig(Ku,w2Mu,'vector') ;
    sigma = real(sqrt(s2)) ;
    [sigma,is] = sort(sigma,'ascend') ;
    Qm = Qm(:,is) ;
    Um = reshape(U_pca*Qm,size(U,1),size(U,2),[]) ;
% Plot eigenvalues
    clf(figure(3)) ; 
%     plot(sigma./max(abs(sigma)),'+') ;
    pls = plot(1-cumsum(sigma)./sum(sigma),'+') ;
    set(gca,'yscale','log') ;
    %set(gca,'ylim',[1e-6 1])
    set(gca,'xlim',[1 numel(sigma)])
%% Animate modes
    Kdir = zeros(numel(sigma),2) ;
    gifOnClick = true ;
    plw = pkg.fem.bloch.waveModeAnimation(unitcell.mesh,Kdir,Um,pls,false,gifOnClick) ;
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE PROPAGATION SOLVING
function [k,w,U,Pr,lq] = dispersion(unitcell,k,w,nModes,perVec,redux)
% Compute the dispersion diagram, either:
    mode = '' ;
% - the frequencies w associated to a given list of wavevectors k [nK nCoord] and w==NaN
    if all(~isnan(k(:))) || any(isnan(w)) ; mode = 'w' ; end
% - the wavenumbers k associated to a given list of wavevectors directions k [nK nCoord] and frequencies w [nK 1]
    if all(~isnan(k(:))) || all(~isnan(w(:))) ; mode = 'knorm' ; end
% - the wavenumbers component ki associated to a given list of wavevectors components k [nK nCoord] with one NaN and frequencies w [nK 1]
    if all(sum(isnan(k),1)==1) ; mode = 'ki' ; end
    if isempty(mode) ; error('Dispersion mode not detectable') ; end
% perVec gives the periodicity vector
    if nargin<5 || isempty(perVec) ; perVec = diag(range(unitcell.mesh.boundingBox,1)) ; end
% apply basis reduction ?
    if nargin<6 ; redux = true ; end
    nCoordK = size(k,2) ;
% Reduction parameters
    if isscalar(redux) % the reduced basis has to be determined
        kr = k(round(linspace(1,end,3)),:) ; % number of 'key' wavevectors
        nModesR = min(nModes+8,2*nModes) ; % extended basis for reduction
        if islogical(redux)  
            tolR = sqrt(eps) ; % default tolerance to cull basis vectors
        else 
            tolR = redux ; 
            redux = true ;  
        end % prescribed tolerance
    else % the reduced basis is given from a previous computation
        Pr = redux ;
        perVec = inf*perVec ; % the reduced basis already includes peridicity...
        redux = false ;
    end
% Build FEM Matrices
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(unitcell.mesh,unitcell.C,unitcell.rho,[],perVec) ;
% Keep only needed peridicity coordinates
    if nCoordK<3
        K0i = K0i(1:nCoordK) ;
        Kij = Kij(1:nCoordK,1:nCoordK) ;
    end
% For plane problems, remove the 3rd displacement component...
    if unitcell.mesh.nCoord<3 || range(unitcell.mesh.Nodes(:,3))==0
        M = M(1:2*end/3,1:2*end/3) ;
        P = P(1:2*end/3,1:2*end/3) ;
        K00 = K00(1:2*end/3,1:2*end/3) ;
        K0i = cellfun(@(Ki)Ki(1:2*end/3,1:2*end/3),K0i,'uni',false) ;
        Kij = cellfun(@(Ki)Ki(1:2*end/3,1:2*end/3),Kij,'uni',false) ;
    end
% Vectorize
    K00v = K00(:) ;
    K0iv = cellfun(@(Ki)reshape(Ki,[],1),K0i,'uni',false) ; K0iv = cat(2,K0iv{:}) ;
    Kijv = cellfun(@(Ki)reshape(Ki,[],1),Kij,'uni',false) ; Kijv = cat(2,Kijv{:}) ;
% Symmetrize
    M = .5*(M+M') ;
% Compute the reduced basis
    nDOF = size(M,1) ;
    if redux
        nKr = size(kr,1) ;
        Vr = zeros(nDOF,nModesR,nKr) ;
    % Compute the modes for the key wavevectors
        wtbr = waitbar(0,'Reduced basis...') ;
        for kk = 1:nKr
            k1 = kr(kk,:) ;
            k2 = reshape(kr(kk,:).'.*kr(kk,:),1,[]) ;
            K = reshape(K00v + sum(K0iv.*k1,2) + sum(Kijv.*k2,2),[nDOF nDOF]) ;
            [Vr(:,:,kk),~] = eigs(K,M,nModesR,'sm') ;
            wtbr = waitbar(kk/nKr,wtbr) ;
        end
        delete(wtbr) ;
    % Orthogonalize the basis
        Vr = Vr(:,:) ;
        Mo = Vr'*(M*Vr) ; Mo = .5*(Mo+Mo') ;
        [Q,lq] = eig(Mo,'vector') ;
    % Sort eigenvectors
        [lq,is] = sort(lq,'descend') ;
        Q = Q(:,is) ;
    % Truncate the basis ...
        if isinteger(tolR) % ... by a prescribed order
            overtol = (1:numel(lq))<=tolR ;
        else % ... by tolerance
            overtol = abs(lq)>tolR*max(abs(lq)) ;
        end
        disp("Ratio of modes kept: "+string(sum(overtol)/numel(overtol))) ;
        lq = lq(overtol) ;
        Q = Q(:,overtol) ;
        Vr = Vr*Q ;
    % Apply the projection
        K00v = reshape(Vr'*(K00*Vr),[],1) ;
        K0iv = cellfun(@(Ki)reshape(Vr'*(Ki*Vr),[],1),K0i,'uni',false) ; K0iv = cat(2,K0iv{:}) ;
        Kijv = cellfun(@(Ki)reshape(Vr'*(Ki*Vr),[],1),Kij,'uni',false) ; Kijv = cat(2,Kijv{:}) ;
        M = Vr'*(M*Vr) ; M = .5*(M+M') ;
        Pr = P*Vr ;
    else
        if ~exist('Pr','var') % no reduction applied (only periodic BC) 
            Pr = P ; 
        else
            K00v = reshape(Pr'*(K00*Pr),[],1) ;
            K0iv = cellfun(@(Ki)reshape(Pr'*(Ki*Pr),[],1),K0i,'uni',false) ; K0iv = cat(2,K0iv{:}) ;
            Kijv = cellfun(@(Ki)reshape(Pr'*(Ki*Pr),[],1),Kij,'uni',false) ; Kijv = cat(2,Kijv{:}) ;
            M = Pr'*(M*Pr) ; M = .5*(M+M') ;
        end
    end
% Compute the dispersion diagram
    nDOF = size(Pr,2) ;
    nK = size(k,1) ;
    w = zeros(nModes,nK) ;
    U = zeros(nDOF,nModes,nK) ;
    %wtbr = waitbar(0,'Dispersion...') ;
    for kk = 1:nK
        k1 = k(kk,:) ;
        k2 = reshape(k(kk,:).'.*k(kk,:),1,[]) ;
        K = reshape(K00v + sum(K0iv.*k1,2) + sum(Kijv.*k2,2),[nDOF nDOF]) ;
        if issparse(M)
            [U(:,:,kk),w2] = eigs(K,M,nModes,'sm') ;
            w2 = diag(w2) ;
        else
            [Q,w2] = eig(K,M,'vector','chol') ;
            [w2,is] = sort(w2,'ascend','comparisonmethod','abs') ;
            U(:,:,kk) = Q(:,is(1:nModes)) ;
            w2 = w2(1:nModes) ;
        end
        w(:,kk) = sqrt(w2) ;
        %wtbr = waitbar(kk/nK,wtbr) ;
    end
    %delete(wtbr) ;
    U = reshape(Pr*U(:,:),unitcell.mesh.nNodes,[],nModes,nK) ;
end





