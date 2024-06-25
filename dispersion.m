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