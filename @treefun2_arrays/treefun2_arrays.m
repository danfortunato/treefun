classdef treefun2_arrays  %#ok<*PROP,*PROPLC>

    properties

        boxes
        n = 16
        domain
        level
        height
        id
        parent
        children
        coeffs
        col
        row
        flatNeighbors
        leafNeighbors

    end

    methods

        function f = treefun2_arrays(varargin)

            if ( nargin < 1 )
                return
            end

            [f.domain, f.level, f.height, f.id, f.parent, f.children, ...
                f.coeffs, f.col, f.row] = deal(varargin{:});
            f.n = size(f.coeffs(end), 1);
            
            % Now do level restriction
            f = balance2(f);
            %[f.flatNeighbors, f.leafNeighbors] = generateNeighbors(f);

        end

    end

    methods ( Access = private )

        function f = refineBox(f, id)

            % Split into four child boxes
            dom = f.domain(:,id);
            xmid = mean(dom(1:2));
            ymid = mean(dom(3:4));
            
            cid1 = length(f.id)+1;
            f.domain(:,cid1)   = [dom(1) xmid dom(3) ymid];
            f.id(cid1)         = cid1;
            f.parent(cid1)     = id;
            f.children(:,cid1) = 0;
            f.level(cid1)      = f.level(id)+1;
            f.height(cid1)     = 0;
            f.coeffs{cid1}     = [];
            f.col(cid1)        = 2*(f.col(id)-1) + 1;
            f.row(cid1)        = 2*(f.row(id)-1) + 1;

            cid2 = length(f.id)+1;
            f.domain(:,cid2)   = [xmid dom(2) dom(3) ymid];
            f.id(cid2)         = cid2;
            f.parent(cid2)     = id;
            f.children(:,cid2) = 0;
            f.level(cid2)      = f.level(id)+1;
            f.height(cid2)     = 0;
            f.coeffs{cid2}     = [];
            f.col(cid2)        = 2*(f.col(id)-1) + 2;
            f.row(cid2)        = 2*(f.row(id)-1) + 1;

            cid3 = length(f.id)+1;
            f.domain(:,cid3)   = [dom(1) xmid ymid dom(4)];
            f.id(cid3)         = cid3;
            f.parent(cid3)     = id;
            f.children(:,cid3) = 0;
            f.level(cid3)      = f.level(id)+1;
            f.height(cid3)     = 0;
            f.coeffs{cid3}     = [];
            f.col(cid3)        = 2*(f.col(id)-1) + 1;
            f.row(cid3)        = 2*(f.row(id)-1) + 2;

            cid4 = length(f.id)+1;
            f.domain(:,cid4)   = [xmid dom(2) ymid dom(4)];
            f.id(cid4)         = cid4;
            f.parent(cid4)     = id;
            f.children(:,cid4) = 0;
            f.level(cid4)      = f.level(id)+1;
            f.height(cid4)     = 0;
            f.coeffs{cid4}     = [];
            f.col(cid4)        = 2*(f.col(id)-1) + 2;
            f.row(cid4)        = 2*(f.row(id)-1) + 2;

            f.children(:,id) = [cid1 cid2 cid3 cid4];
            f.height(id) = 1;
            f.coeffs{id} = [];

        end

        function f = balance2(f)
            L = leaves(f);
            while ( ~isempty(L) )
                id = L(1);
                neighborIDs = neighbors(f, id);
                split = false;

                if ( ~split && ~isnan(neighborIDs(1)) ) % Left neighbor
                    nbr = neighborIDs(1);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children([2 4], nbr); % SE, NE
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(2)) ) % Right neighbor
                    nbr = neighborIDs(2);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children([1 3], nbr); % SW, NW
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(3)) ) % Down neighbor
                    nbr = neighborIDs(3);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children([3 4], nbr); % NW, NE
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(4)) ) % Up neighbor
                    nbr = neighborIDs(4);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children([1 2], nbr); % SW, SE
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(5)) ) % Left down corner
                    nbr = neighborIDs(5);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children(4, nbr); % Right up
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end
                
                if ( ~split && ~isnan(neighborIDs(6)) ) % Right down corner
                    nbr = neighborIDs(6);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children(3, nbr); % Left up
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end
                
                if ( ~split && ~isnan(neighborIDs(7)) ) % Left up corner
                    nbr = neighborIDs(7);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children(2, nbr); % Right down
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end

                if ( ~split && ~isnan(neighborIDs(8)) ) % Right up corner
                    nbr = neighborIDs(8);
                    if ( ~isLeaf(f, nbr) )
                        childrenIDs = f.children(1, nbr); % Left down
                        split = split | any( ~isLeaf(f, childrenIDs) );
                    end
                end

                if ( split )
                    % This was a leaf, so we'll use its coeffs to evaluate
                    % on the new children
                    coeffs = f.coeffs{id};
                    % Split into four child boxes
                    f = refineBox(f, id);
                    children = f.children(:,id);
                    [LL, LR, UL, UR] = coeffs2children(coeffs);
                    f.coeffs{children(1)} = LL; % Lower left
                    f.coeffs{children(2)} = LR; % Lower right
                    f.coeffs{children(3)} = UL; % Upper left
                    f.coeffs{children(4)} = UR; % Upper right
                    L = [L children(:).']; %#ok<AGROW>
                    neighborIDs = neighborIDs( ~isnan(neighborIDs) );
                    idx = f.level(neighborIDs) < f.level(id);
                    L = [L neighborIDs(idx)]; %#ok<AGROW>
                    L = unique(L, 'stable');
                end

                L(1) = [];
            end
        end
        
    end

end
