%% SIMULATION ASSOCIATED TO THE WIDEBAND CHARACTERIZATION OF AN ARCHITECTURED MATERIAL
%clc 
clear all

%% UNIT CELL DEFINITION
% Define the unit cell as a structure containing at least the fields:
%   - "mesh" as the 3D unit cell mesh (must be periodic)
%   - the material properties ("C" & "rho") at the mesh quadrature points

unitcell = [] ;
% Unit cell geometry & discretization
    L = .5*[1 1 1] ; dx = 1/50 ; % unit cell size & element length (mm)
    %L = [dx dx L(3)] ;
    unitcell.mesh = pkg.geometry.mesh.GridMesh(L,dx) ;
% Use quadratic elements improves a lot the convergence
    unitcell.mesh.setElementTypes(pkg.geometry.mesh.elements.LagrangeElement('hex',2)) ;
% Two-phase geometry (defined by a levelset function at nodes)
% lvlst<0 for phase 1 and lvlst>0 for phase 2 (lvlst==0 is the interface)
    switch 1
        case 0 % uniform phase=phase1
            unitcell.lvlst = -ones(unitcell.mesh.nNodes,1) ;
        case 1 % spherical inclusion(s)
            p = 2 ; 4/2 ; % power (inf for cubes, 2 for spheres, 1 for diamonds
            R = 50/100*min(L)/2 ; % radius 
            xc = L.*cat(3,[.5 .5 .5] ...
                        ...,[0 0 0] ,[0 1 0] ,[1 1 0] ,[1 0 0] ,[0 0 1] ,[0 1 1] ,[1 1 1] ,[1 0 1] ...
                        ) ;
            rx2 = min(sum(abs(unitcell.mesh.Nodes-xc).^p,2),[],3) ;
            unitcell.lvlst = ((R.^p)-rx2) ; % phase 1 inside
        case 2 % 2-phases gyroid
            t = 20/100 ; % gyroid thickness
            x = 2*pi.*unitcell.mesh.Nodes./L ;
            gyrFcn = sin(x(:,1)).*cos(x(:,2)) + sin(x(:,2)).*cos(x(:,3)) + sin(x(:,3)).*cos(x(:,1)) ;
            unitcell.lvlst = abs(gyrFcn)-pi*t ; % phase1 is the gyroid wall
        case 3 % 2-phases schwartz's surface
            t = 2*dx ; % gyroid thickness
            x = 2*pi.*unitcell.mesh.Nodes./L ;
            schwFcn = sum(cos(x),2) ;
            unitcell.lvlst = abs(schwFcn)-pi*t ; % phase1 is the wall
    end
    cut = unitcell.mesh.simplex.cut(unitcell.lvlst) ; % the cut is performed on the simplex mesh for speed
    clf reset ; view(60,30) ;  axis equal ; 
    plot(unitcell.mesh,'VisibleEdges','none','FaceAlpha',0) ; 
    plot(cut.ON,'VisibleEdges','none') ; 
    light ;
% Material phase properties
    unitcell.E_p = [4e3*(1+0*.01i) 2e3*(1+0*.05i)] ; % young moduli of phases (MPa)
    unitcell.nu_p = [.3 .45] ; % poisson ratio of phases
    unitcell.rho_p = [1220e-12 1215e-12] ; % density of phases (tons/mm3)
    unitcell.G_p = unitcell.E_p./2./(1+unitcell.nu_p) ;
    unitcell.C_p = pkg.fem.bloch.stiffness(unitcell.E_p,unitcell.G_p) ;
    unitcell.S_p = pkg.math.inv(unitcell.C_p) ;
% Distribute on the mesh
    [ee,we,ie] = unitcell.mesh.integration() ; 
    Ne = unitcell.mesh.interpMat(ee,ie) ; % interpolation at quadrature points
    unitcell.phase = .5*(sign(Ne*sign(unitcell.lvlst))+1) + 1 ;
    unitcell.C = unitcell.C_p(:,:,unitcell.phase) ;
    unitcell.rho = unitcell.rho_p(unitcell.phase) ;
% Estimate mean properties
    unitcell.volFrac = full(sum(sparse((1:numel(ie))',unitcell.phase,we,numel(ie),numel(unitcell.E_p)),1))/sum(we(:)) ;
    unitcell.C_voigt = sum(unitcell.C_p.*reshape(unitcell.volFrac,1,1,[]),3) ;
    unitcell.C_reuss = inv(sum(unitcell.S_p.*reshape(unitcell.volFrac,1,1,[]),3)) ;
    unitcell.rho_mean = sum(unitcell.rho_p.*unitcell.volFrac,2) ;
% Estimate isotropic phase velocities
    unitcell.cl_voigt = sqrt(unitcell.C_voigt(1,1)/unitcell.rho_mean) ;
    unitcell.ct_voigt = sqrt(unitcell.C_voigt(6,6)/unitcell.rho_mean) ;
    unitcell.cl_reuss = sqrt(unitcell.C_reuss(1,1)/unitcell.rho_mean) ;
    unitcell.ct_reuss = sqrt(unitcell.C_reuss(6,6)/unitcell.rho_mean) ;
% Display
    unitcell
    
    
%% WIDE BAND WAVE PROPAGATION IN THE 3D SOLID
    k = [1 0 0].*pi/norm(L(1)).*linspace(0.01,1,100)' ; % wavevectors (rad/mm) [nK 3]
    nModes = 20 ; % extract the first nModes
% Reduction parameters
    redux = true ; % apply basis reduction ?
    [k,w,U] = dispersion(unitcell,k,NaN,nModes,[],redux) ;
% Display
    clf ;
    knorm = sqrt(sum(abs(k).^2,2)) ;
    clrs = get(gca,'colororder') ;
    fill([knorm ; flip(knorm)],[unitcell.cl_voigt.*knorm ; unitcell.cl_reuss.*flip(knorm)],clrs(1,:),'linestyle','none','facealpha',.3,'displayname','$c_\ell$')
    fill([knorm ; flip(knorm)],[unitcell.ct_voigt.*knorm ; unitcell.ct_reuss.*flip(knorm)],clrs(2,:),'linestyle','none','facealpha',.3,'displayname','$c_t$')
    plot(repmat(knorm,[nModes 1]),real(w(:)),'.k','displayname','SAFE') ;
    xlabel 'Wavenumber (rad/mm)'
    ylabel 'Frequency (rad/s)'
    legend('location','southeast')
    title('\bf Waves in the 3D solid') ;
% Detect Longi & shear modes
    Uu = kron(eye(3),ones(unitcell.mesh.nNodes,1))/sqrt(unitcell.mesh.nNodes) ;
    Uv = reshape(U,[],nModes*size(k,1)) ;
    MAC = reshape((Uu'*Uv)./sqrt(sum(abs(Uv).^2,1)),[size(Uu,2) nModes size(k,1)]) ;
    [~,imax] = max(abs(MAC),[],2) ; 
    wu = w(sub2ind(size(w),repmat((1:size(k,1))',[1 size(Uu,2)]),permute(imax,[3 1 2]))) ;
    pl = plot(knorm,real(wu),'o','markersize',4,'linewidth',1) ;
    set(pl,{'displayname'},{'$\ell$-like';'$t_2$-like';'$t_3$-like'}) ;
    cl_eq = wu(:,1)./knorm ; % equivalent longitudinal velocity
    ct_eq = mean(wu(:,2:3),2)./knorm ; % equivalent longitudinal velocity
% Rescaling
    fc = .5*sqrt(unitcell.ct_voigt.*unitcell.ct_reuss)/L(1) ; % cutoff frequency
    kc = pi/L(1) ; % cutoff wavenumber
    obj = allchild(gca) ;
    set(obj,{'xdata'},cellfun(@(x)x./kc,get(obj,'XData'),'uni',false))
    set(obj,{'ydata'},cellfun(@(x)x./2/pi/fc,get(obj,'YData'),'uni',false))
    xlabel 'Normalized Wavenumber $\frac{k \times L}{\pi}$' ; 
    ylabel 'Normalized Frequency $\frac{f \times 2 L}{c_t}$'
    
%% Equivalent moduli
    G_eq = unitcell.rho_mean.*ct_eq.^2 ;
    C11_eq = unitcell.rho_mean.*cl_eq.^2 ;
    clf ;
    plot(knorm,real(G_eq),'.','displayname','$\bar{G}$')
    plot(knorm,real(C11_eq),'.','displayname','$\bar{C}_{11}$')
    plot(knorm,unitcell.G_p+knorm.*0,'.','displayname','$G_\phi$')
    plot(knorm,reshape(unitcell.C_p(1,1,:),1,[])+knorm.*0,'.','displayname','$C_{\phi,11}$')
    xlabel 'Wavenumber (rad/mm)'
    ylabel 'Moduli (MPa)'
    legend
    
    
%% RIBBON-SHAPED WAVEGUIDE
    ribbon = unitcell ; 
    N = [10 2] ; % number of unit cells in the ribbon [horizontal vertical]
    ribbon.b = N(1)*L(2) ; % ribbon width
    ribbon.h = N(2)*L(3) ; % ribbon height
% Replicate the mesh
    ribbon.mesh = copy(unitcell.mesh) ;
    repMesh = pkg.geometry.mesh.GridMesh(uint64([0 N-1])) ;
    ribbon.mesh.replicate(ribbon.mesh.Nodes+permute(repMesh.Nodes.*L,[3 2 1]),false) ;
    [ribbon.mesh,~,~,meanMat] = ribbon.mesh.cullDuplicates ;
% Replicate the properties
    ribbon.C = repmat(unitcell.C,[1 1 prod(N)]) ;
    ribbon.rho = repmat(unitcell.rho,[1 1 prod(N)]) ;
    ribbon.lvlst = meanMat*repmat(unitcell.lvlst,[prod(N) 1]) ;
% Display
    clf ; view(60,30) ; axis equal tight ; %plot(ribbon.mesh) ; return
    cut = ribbon.mesh.simplex.cut(ribbon.lvlst) ; % the cut is performed on the simplex mesh for speed
    clf reset ; view(60,30) ;  axis equal ; 
    plot(ribbon.mesh,'VisibleEdges','none','FaceAlpha',0) ; 
    plot(cut.IN,'VisibleEdges','none') ; 
    light ;
    
%% COMPUTE THE RIBBON DISPERSION DIAGRAM
    k = [1 0 0].*pi/ribbon.h.*linspace(0.01,1,100)' ; % wavevectors (rad/mm) [nK 3]
    nModes = 20 ; % extract the first nModes
% Reduction parameters
    redux = true ; % apply basis reduction ?
    [k,w,U] = dispersion(ribbon,k,NaN,nModes,L.*[1 0 0],redux) ;
% Display
    clf ;
    knorm = sqrt(sum(abs(k).^2,2)) ;
    plot(repmat(knorm,[nModes 1]),real(w(:)),'.k','displayname','SAFE') ;
    xlabel 'Wavenumber (rad/mm)'
    ylabel 'Frequency (rad/s)'
    legend('location','southeast')
    title('\bf Waves in the ribbon') ;
    
    
%% ULTRASONIC THETA-SCAN
    Nc = 8 ; % number of unit cells in the specimen thickness
    w = 2*pi*2.25e6 ; linspace(1e6,10e6,100) ; % frequency (rad/s)
    theta1 = (0:64)*pi/180 ; % first indicence angle (rad)
    theta2 = 0*pi/180 ; % second indicence angle (rad)
    rho0 = 1000e-12 ; % fluid density (tons/mm3)
    c0 = 1500e3 ; % fluid wave speed (mm/s) ;
    aI = 1 ; % incident wave amplitude
    plate = unitcell ;
% Replicate the mesh
    plate.mesh = copy(unitcell.mesh) ;
    plate.mesh.replicate(plate.mesh.Nodes+[0 0 L(3)].*reshape(0:Nc-1,1,1,Nc),false) ;
    [plate.mesh,~,~,meanMat] = plate.mesh.cullDuplicates ;
    plate.h = Nc*L(3) ;
% Replicate the properties
    plate.C = repmat(unitcell.C,[1 1 Nc]) ;
    plate.rho = repmat(unitcell.rho,[Nc 1]) ;
    plate.lvlst = meanMat*repmat(unitcell.lvlst,[Nc 1]) ;
% Display
    clf ; view(60,30) ; axis equal tight ;
    cut = plate.mesh.simplex.cut(plate.lvlst) ; % the cut is performed on the simplex mesh for speed
    clf reset ; view(60,30) ;  axis equal ; 
    plot(plate.mesh,'VisibleEdges','none','FaceAlpha',0) ; 
    plot(cut.ON,'VisibleEdges','none') ; 
    light ;
%% Compute the ultrasonic theta-scan..
    [aT,aR,u] = us(plate,aI,w,theta1,theta2,c0,rho0) ;
    aBalanceErr = abs(aI).^2 - abs(aT).^2 - abs(aR).^2 ;
            
%% Homogeneous equivalent slab validation
    [aTiv,aRiv] = us_hom_iso(plate.cl_voigt,plate.ct_voigt,plate.rho_mean,plate.h,aI,w,theta1,theta2,c0,rho0) ;
    [aTir,aRir] = us_hom_iso(plate.cl_reuss,plate.ct_reuss,plate.rho_mean,plate.h,aI,w,theta1,theta2,c0,rho0) ;
%     [aTav,aRav] = us_hom_aniso(plate.C_voigt,plate.rho_mean,plate.h,aI,w,theta1,theta2,c0,rho0) ;
%     [aTar,aRar] = us_hom_aniso(plate.C_reuss,plate.rho_mean,plate.h,aI,w,theta1,theta2,c0,rho0) ;
    
% Display
    clf ; 
    if numel(unique(theta1+theta2))==1
        p = w/2e6/pi ; xlabel 'Frequency (MHz)' ;
    elseif numel(unique(w))==1
        p = (theta1+theta2)*180/pi ; xlabel 'Incidence Angle (deg)' ;
    end
    clrs = get(gca,'colororder') ;
    fill([p flip(p)]',[abs(aTiv) flip(abs(aTir))]',clrs(1,:),'linestyle','-','edgecolor',clrs(1,:),'facealpha',.3,'displayname','$|a_T|$')
    fill([p flip(p)]',[abs(aRiv) flip(abs(aRir))]',clrs(2,:),'linestyle','-','edgecolor',clrs(2,:),'facealpha',.3,'displayname','$|a_R|$')
    plot(p,abs(aT),'.','color',clrs(1,:)) ; plot(p,abs(aR),'.','color',clrs(2,:)) ;  
    set(findobj(gca,'type','line'),'markersize',15,'linewidth',1.5)
    legend  
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAVE PROPAGATION SOLVING
function [k,w,U] = dispersion(unitcell,k,w,nModes,perVec,redux)
% Compute the dispersion diagram, either:
    mode = '' ;
% - the frequencies w associated to a given list of wavevectors k [nK nCoord] and w==NaN
    if all(~isnan(k(:))) || any(isnan(w)) ; mode = 'w' ; end
% - the wavenumbers k associated to a given list of wavevectors directions k [nK nCoord] and frequencies w [nK 1]
    if all(~isnan(k(:))) || all(~isnan(w(:))) ; mode = 'knorm' ; end
% - the wavenumbers component ki associated to a given list of wavevectors components k [nK nCoord] with one NaN and frequencies w [nK 1]
    if all(sum(isnan(k),1)==1) || all(~isnan(w(:))) ; mode = 'ki' ; end
    if isempty(mode) ; error('Dispersion mode not detectable') ; end
% perVec gives the periodicity vector
    if nargin<5 || isempty(perVec) ; perVec = diag(range(unitcell.mesh.boundingBox,1)) ; end
% apply basis reduction ?
    if nargin<6 ; redux = true ; end
% Reduction parameters
    kr = k(round(linspace(1,end,3)),:) ; % number of 'key' wavevectors
    nModesR = min(nModes+8,2*nModes) ; % extended basis for reduction
    tolR = sqrt(eps) ; % tolerance to cull basis vectors
% Prepare the FEM Matrices
    [K00,K0i,Kij,M,P] = pkg.fem.bloch.FEM(unitcell.mesh,unitcell.C,unitcell.rho,[],perVec) ;
    M = .5*(M+M') ;
    K00v = K00(:) ;
    K0iv = cellfun(@(Ki)reshape(Ki,[],1),K0i,'uni',false) ; K0iv = cat(2,K0iv{:}) ;
    Kijv = cellfun(@(Ki)reshape(Ki,[],1),Kij,'uni',false) ; Kijv = cat(2,Kijv{:}) ;
    nDOF = size(M,1) ;
% Compute the reduced basis
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
        overtol = abs(lq)>tolR*max(abs(lq)) ;
        disp("Ratio of modes kept: "+string(sum(overtol)/numel(overtol))) ;
        Q = Q(:,overtol) ;
        Vr = Vr*Q ;
    % Apply the projection
        nDOF = size(Vr,2) ;
        K00v = reshape(Vr'*(K00*Vr),[],1) ;
        K0iv = cellfun(@(Ki)reshape(Vr'*(Ki*Vr),[],1),K0i,'uni',false) ; K0iv = cat(2,K0iv{:}) ;
        Kijv = cellfun(@(Ki)reshape(Vr'*(Ki*Vr),[],1),Kij,'uni',false) ; Kijv = cat(2,Kijv{:}) ;
        M = Vr'*(M*Vr) ; M = .5*(M+M') ;
        P = P*Vr ;
    end
% Compute the dispersion diagram
    nK = size(k,1) ;
    w = zeros(nK,nModes) ;
    U = zeros(nDOF,nModes,nK) ;
    wtbr = waitbar(0,'Dispersion...') ;
    for kk = 1:nK
        k1 = k(kk,:) ;
        k2 = reshape(k(kk,:).'.*k(kk,:),1,[]) ;
        K = reshape(K00v + sum(K0iv.*k1,2) + sum(Kijv.*k2,2),[nDOF nDOF]) ;
        if ~redux
            [U(:,:,kk),w2] = eigs(K,M,nModes,'sm') ;
            w2 = diag(w2) ;
        else
            [Q,w2] = eig(K,M,'vector','chol') ;
            [w2,is] = sort(w2,'ascend','comparisonmethod','abs') ;
            U(:,:,kk) = Q(:,is(1:nModes)) ;
            w2 = w2(1:nModes) ;
        end
        w(kk,:) = sqrt(w2) ;
        wtbr = waitbar(kk/nK,wtbr) ;
    end
    delete(wtbr) ;
    U = reshape(P*U(:,:),unitcell.mesh.nNodes,[],nModes,nK) ;
end


function [aT,aR,u] = us(plate,aI,w,theta1,theta2,c0,rho0)
% Each case...
    w = w + 0*theta1 + 0*theta2 + 0*aI ;
    theta1 = 0*w + theta1 + 0*theta2 + 0*aI ;
    theta2 = 0*w + 0*theta1 + theta2 + 0*aI ;
    aI = 0*w + 0*theta1 + 0*theta2 + aI ;
% Wavenumbers
    kappa1 = (w./c0).*sin(theta1) ;
    kappa2 = (w./c0).*sin(theta2) ;
    k0 = sqrt((w./c0).^2-kappa1.^2-kappa2.^2) ;
% Interpolation Matrices
    % Full mesh bulk
        [ee,we,ie] = plate.mesh.integration() ;
        nQP = numel(we) ;
        O = sparse(nQP,plate.mesh.nNodes) ;
        N = plate.mesh.interpMat(ee,ie) ;
        D = plate.mesh.diffMat(ee,ie) ; 
        [D{end+1:3}] = deal(O) ;
        x = plate.mesh.Nodes ;
    % Input face
        fi = plate.mesh.Faces.subpart(plate.mesh.near([NaN NaN 0],plate.mesh.Faces)) ;
        [ei,wi,ii] = plate.mesh.integration(fi) ;
        Ni = plate.mesh.interpMat(ei,ii,fi) ;
        xi = Ni*plate.mesh.Nodes(:,1:2) ;
    % Output face
        fo = plate.mesh.Faces.subpart(plate.mesh.near([NaN NaN plate.h],plate.mesh.Faces)) ;
        [eo,wo,io] = plate.mesh.integration(fo) ;
        No = plate.mesh.interpMat(eo,io,fo) ;
        xo = No*plate.mesh.Nodes(:,1:2) ;
% Constant system matrices
    % Strains E = [E11;E22;E33;2E13;2E23;2E12] = B*[u1;u2;u3]
        B = [D{1} O O ; O D{2} O ; O O D{3} ; D{3} O D{1} ; O D{3} D{2} ; D{2} D{1} O] ;
    % Stiffness matrix
        Cqp = plate.C.*reshape(we,1,1,[]) ; % times weights
        iii = (0:5)'*nQP + zeros(1,6) + reshape(1:nQP,1,1,[]) ;
        jjj = zeros(6,1) + (0:5)*nQP + reshape(1:nQP,1,1,[]) ;
        CW = sparse(iii(:),jjj(:),Cqp(:),6*nQP,6*nQP) ;
        K = B'*CW*B ;
    % Mass Matrix
        rhoW = plate.rho(:).*we(:) ;
        M = (N'*diag(sparse(rhoW))*N) ;
        M = blkdiag(M,M,M) ;
% Solve for each case..
    u = zeros([plate.mesh.nNodes 3 size(aI)]) ;
    u3om = zeros(size(aI)) ; % mean output normal displacement
    u3im = zeros(size(aI)) ; % mean input normal displacement
    wtbr = waitbar(0,'US SCAN...') ; tttt = tic ;
    for pp = 1:numel(w)
    % Acoustic radiation
        % Impedance
            phii = exp(-1i.*(kappa1(pp).*xi(:,1)+kappa2(pp).*xi(:,2))) ;
            phio = exp(-1i.*(kappa1(pp).*xo(:,1)+kappa2(pp).*xo(:,2))) ;
            Z = w(pp).*rho0./k0(pp).*kron(diag(sparse([0 0 1])),Ni'*diag(sparse(wi))*Ni + No'*diag(sparse(wo))*No) ;
            f = 2*aI(pp)*kron([0 0 1]',Ni'*(wi.*phii)) ;
    % Periodicity
        P = plate.mesh.perNodeMat([range(x(:,1)) 0 0 ; 0 range(x(:,2)) 0]) ;
        phi = exp(-1i*(kappa1(pp).*x(:,1)+kappa2(pp).*x(:,2))) ;
        P = diag(sparse(phi))*P ;
        P(:,sum(P,1)==0) = [] ; % delete slave DOFs
        P = blkdiag(P,P,P) ;
        Kp = P'*K*P ;
        Mp = P'*M*P ;
        Zp = P'*Z*P ;
        fp = P'*f ;
    % Solve
        D = Kp + 1i*w(pp)*Zp - w(pp)^2*Mp ;
        us = D\fp ;
        u(:,:,pp) = reshape(P*us,plate.mesh.nNodes,3) ;
%         clf ; axis equal tight
%         pl = plot(plate.mesh...
%                 ,'Deformation',real(u)...
%                 ,'CData',real(u(:,3))...
%                 ,'VisibleFaces','all'...
%                 ,'VisibleEdges','none'...
%                 ,'FaceAlpha',.15...
%                 ) ; 
%         colorbar
%         view([30 30])
    % Input/Output normal displacements
        u3om(pp) = ((wo./phio).'*No*u(:,3,pp))/sum(wo) ;
        u3im(pp) = ((wi./phii).'*Ni*u(:,3,pp))/sum(wi) ;
    % Update waitbar
        if toc(tttt)>.1 ; wtbr = waitbar(pp/numel(w),wtbr) ; tttt = tic ; end
    end
    delete(wtbr) ;
% Evaluate transmitted & reflected waves
    aT = (1i*rho0*w.^2./k0.*u3om.*exp(1i*plate.h*k0)) ;
    aR = (aI - 1i*rho0*w.^2./k0.*u3im) ;
end

function [aT,aR] = us_hom_iso(cl,ct,rho,h,aI,w,theta1,theta2,c0,rho0)
% Analytic model for an isotropic homogeneous slab
% Wavenumbers
    kappa1 = (w./c0).*sin(theta1) ;
    kappa2 = (w./c0).*sin(theta2) ;
    k0 = sqrt((w./c0).^2-kappa1.^2-kappa2.^2) ;
    kappa = sqrt(kappa1.^2 + kappa2.^2) ;
    kl = sqrt((w./cl).^2-kappa.^2) ;
    kt = sqrt((w./ct).^2-kappa.^2) ;
% Cosines & Sines
    cl = cos(kl.*h/2) ; sl = sin(kl.*h/2) ;
    ct = cos(kt.*h/2) ; st = sin(kt.*h/2) ;
    pi0 = exp(1i*h*k0) ;
% Terms
    alpha = (kt.^2-kappa.^2).^2 ;
    beta = 4*kappa.^2.*kl.*kt ;
    gamma = (rho0./rho).*(kl./k0).*(kt.^2+kappa.^2).^2 ;
% Sym/Antisym
    aS = .5*pi0.*aI.*((alpha.*cl.*st+beta.*sl.*ct)-1i*gamma.*sl.*st)./((alpha.*cl.*st+beta.*sl.*ct)+1i*gamma.*sl.*st) ;
    aA = .5*pi0.*aI.*(gamma.*cl.*ct-1i*(alpha.*sl.*ct+beta.*cl.*st))./(gamma.*cl.*ct+1i*(alpha.*sl.*ct+beta.*cl.*st)) ;
% Transmitted/reflected
    aT = (aS+aA) ;
    aR = (aS-aA).*exp(-1i*h*k0) ;
end

function [aT,aR] = us_hom_aniso(Cij,rho,h,aI,w,theta1,theta2,c0,rho0)
% Analytic model for an ANisotropic homogeneous slab
    %Cij = plate.C_voigt ; rho = plate.rho_mean ; h = plate.h ;
% Fluid Wavenumbers
    kappa1 = (w./c0).*sin(theta1) ;
    kappa2 = (w./c0).*sin(theta2) ;
    k0 = sqrt((w./c0).^2-kappa1.^2-kappa2.^2) ;
    sz = size(k0) ; nsz = numel(sz) ;
% Reshaping by letting the 4 first dimensions free
    sz4 = [1 1 1 1 sz] ;
    let4 = @(x)reshape(x,[1 1 1 1 size(x)]) ; % [1 1 1 1 sz]
    AI = let4(aI) ; 
    W = let4(w) ; 
    K1 = let4(kappa1) ; 
    K2 = let4(kappa2) ; 
    K0 = let4(k0) ;
% Solid wavenumbers: first eigenvalue problem
% governing equations: S_ij,j = (C_ijkl u_k,l),j = C_ijkl u_k,lj = - k_l k_j C_ijkl u_k = -rho w^2 u_i
    % Convert C to 4-th order tensor
        % Indices in Voigt notation
        iv = [1 1 ; 2 2 ; 3 3 ; 2 3 ; 1 3 ; 1 2] ; % [6 2]
        iv = reshape(iv,6,1,2) + zeros(1,6) ; % [6 6 2]
        iv = cat(3,iv,permute(iv,[2 1 3])) ; % [6 6 4]
        % Include symetries
        iv_s = [iv ; iv(:,:,[2 1 3 4]) ; iv(:,:,[1 2 4 3]) ; iv(:,:,[2 1 4 3])] ; % [24 6 4]
        Cij_s = repmat(Cij,[4 1]) ; % [24 6 ?]
        % Build the 4th order tensor
        Cijkl = zeros([81 size(Cij,3:max(ndims(Cij),3))]) ; % [81 ?]
        ind = sub2ind([3 3 3 3],iv_s(:,:,1),iv_s(:,:,2),iv_s(:,:,3),iv_s(:,:,4)) ; % [24*6 1]
        Cijkl(ind,:) = reshape(Cij_s(:,:,:),24*6,[]) ; % [81 ?]
        Cijkl = reshape(Cijkl,[3 3 3 3 size(Cijkl,2:ndims(Cijkl))]) ; % [3 3 3 3 ?]
    % Acoustic tensor decomposition 
    % Tik = k_l k_j C_ijkl 
    %     = kappa_a kappa_b C_ibka + k_3 kappa_b C_ibk3 +  kappa_a k_3 C_i3ka + k_3 k_3 C_i3k3
    %     = T0(kappa) + k_3 T1(kappa) + k_3^2 T2(kappa)
        Cikjl = permute(Cijkl,[1 3 2 4 5:ndims(Cijkl)]) ;  % [3 3 3 3 ?] easier contractions
        kappa = cat(1,K1+zeros(sz4),K2+zeros(sz4),zeros(sz4)) ; % [3 1 1 1 sz]
        ik3 = cat(1,zeros(sz4),zeros(sz4),ones(sz4)) ;
        T0 = sum(reshape(kappa,[1 1 3 1 sz]).*reshape(kappa,[1 1 1 3 sz]).*Cikjl,[3 4]) ; % [3 3 1 1 sz] 
        T1 = sum(reshape(ik3,[1 1 3 1 sz]).*reshape(kappa,[1 1 1 3 sz]).*Cikjl,[3 4]) ...
             + sum(reshape(kappa,[1 1 3 1 sz]).*reshape(ik3,[1 1 1 3 sz]).*Cikjl,[3 4]) ; % [sz 3 3]  
        T2 = sum(reshape(ik3,[1 1 3 1 sz]).*reshape(ik3,[1 1 1 3 sz]).*Cikjl,[3 4]) ; % [sz 3 3] 
    % Solve the quadratic eigenvalue problem in k3
    % (T0 + k_3 T1 + k_3^2 T2 - w^2.rho.eye(3)).uk = 0
        w2rI = W.^2.*rho*eye(3) ; % [3 3]
        array2matrices = @(A)num2cell(A,[1 2]) ;
        matrices2array = @(M)reshape(cat(4,M{:}),[size(M{1}) size(M,3:max(ndims(M),3))])  ;
        A0 = array2matrices(T0-w2rI) ; % {sz4} of [3 3]
        A1 = array2matrices(T1) ; % {sz4} of [3 3]
        A2 = array2matrices(T2) ; % {sz4} of [3 3]
        [U,k3] = cellfun(@polyeig,A0,A1,A2,'uni',false) ; % {sz4} of U:[3 6],k3:[6 1]
        U = matrices2array(U) ; % [3 6 1 1 sz]
        k3 = matrices2array(k3) ; % [6 1 1 1 sz]
        k3 = permute(k3,[2 1 3:ndims(k3)]) ; % [1 6 1 1 sz] eases what follows
% Apply boundary conditions to find the amplitudes
% two interfaces at x3 = +/- h/2, outgoing normals n = -/+ e_3
    x3 = .5*h*cat(4,-1,1) ; % [1 1 1 2]
% Every field is modulated by exp(-1i*(kappa_a.x_a) thus this term is discarded
% Fluid pressure: p = a_+ exp(-1i*k0.x_3) + a_- exp(1i*k0.x_3) = Pe*[aR;aT;aI]
    % x3 <= -h/2: a^+ = aI, a^- = aR
    % x3 >= +h/2: a^+ = aT, a^- = 0
    e0p = exp(-1i.*K0.*x3) ; % [1 1 1 2 sz]
    e0m = exp(1i.*K0.*x3) ; % [1 1 1 2 sz]
    Pe = e0p.*cat(4,[0 0 1],[0 1 0]) + e0m.*cat(4,[1 0 0],[0 0 0]) ; % [1 3 1 2 sz]
    % dp_dx3 = -1i*k0 [ a_+ exp(-1i*k0.x_3) - a_- exp(1i*k0.x_3)] = d3Pe*[aR;aT;aI]
    d3Pe = -1i*K0.*(e0p.*cat(4,[0 0 1],[0 1 0]) - e0m.*cat(4,[1 0 0],[0 0 0])) ; % [1 3 1 2 sz]
% Solid fields: 
%   u = sum_r{ a_r U^r exp(-1i*k3^r*x3) } = Ue*a
    Ue = U.*exp(-1i.*k3.*x3) ; % [3 6 1 2 sz]
%   grad(u) = u_i,j = -1i * sum_r{ a_r U^r_i k^r_j exp(-1i*k3^r*x3) } = Ge*a
    K = kappa + ik3.*k3 ; % [3 6 1 1 sz]
    Ge = -1i*Ue.*permute(K,[3 2 1 4 5:ndims(K)]) ; % [3 6 3 2 sz]
%   Eij = sym(grad(u)) = Ee*a
    Ee = .5*(Ge+permute(Ge,[3 2 1 4 5:ndims(Ge)])) ; % [3 6 3 2 sz] 
%   Sij = Cijkl Ekl = Se*a
    Se = sum(...
            permute(Cijkl,[1 nsz+5 2 nsz+6 4+(1:nsz) 3 4]) ... [3 1 3 1 sz 3 3]
            .*permute(Ee,[nsz+5 2 nsz+6 4 4+(1:nsz) 1 3]) ... [1 6 1 2 sz 3 3]
          ,nsz+4+[1 2]) ; % [3 6 3 2 sz]
% Assemble the system M*[a_r;aR;aT;aI] = 0
    Ms = [] ; % for solid wave amplitudes a_r
    M0 = [] ; % for fluid wave amplitudes [aR;aT;aI]
    ss = prod(sz) ; 
% stress continuity: S_ij.n_j = -p n_i 
    % S13 = 0 at both interfaces
    Ms = [Ms ; Se(1,:,3,1,:) ; Se(1,:,3,2,:)] ; 
    M0 = [M0 ; zeros([2 3 1 1 ss])] ;
    % S23 = 0 at both interfaces
    Ms = [Ms ; Se(2,:,3,1,:) ; Se(2,:,3,2,:)] ;
    M0 = [M0 ; zeros([2 3 1 1 ss])] ;
%   S_33 + p = 0 at both interfaces
    Ms = [Ms ; Se(3,:,3,1,:) ; Se(3,:,3,2,:)] ; 
    M0 = [M0 ; Pe(:,:,:,1,:) ; Pe(:,:,:,2,:)] ;
% displacement continuity: w^2.rho0.u_i.n_i = dp/dx_i.n_i
%   w^2.rho0.u_3 - dp/dx_3 = 0
    w2rU = rho0.*w.^2.*Ue ;
    Ms = [Ms ; w2rU(3,:,:,1,:) ; w2rU(3,:,:,2,:)] ; 
    M0 = [M0 ; -d3Pe(:,:,:,1,:) ; -d3Pe(:,:,:,2,:)] ;
% Full system
    M = [Ms M0] ;
% Solve for given aI
    A = M(:,1:8,:) ;
    b = -M(:,9,:).*reshape(aI,1,1,[]) ;
    a = cellfun(@mldivide,array2matrices(A),array2matrices(b),'uni',false) ;
    a = reshape(cat(3,a{:}),[8 1 sz]) ;
% Transmitted & reflected
    aT = reshape(a(8,:),sz) ;
    aR = reshape(a(7,:),sz) ;
end




