classdef SOLN < handle
    %Holds solution and its properties
    
    properties
        C = [];
        S = [];
        D = [];
        mat = [];
        vec = [];
        mesh = [];
        x;
        v;
        k;
        per = [];
        iper = [];
    end
    
    methods
        function obj = SOLN(x,v,k,per,iper)
            obj.x = x;
            obj.v = v;
            obj.k = k;
            if nargin > 3
                obj.per = per;
                obj.iper = iper;
            end
        end
        
        
        function [] = getMatFromLR(obj)
            %Create full rank matrix from low rank factors
            obj.mat = obj.C*obj.S*obj.D';
        end
        
        
        function [] = getLRFactors(obj,r)
            %Create low rank factors from full matrix          
            assert(~isempty(obj.mat),'Error: obj.mat Needed for SVD');
            [obj.C,obj.S,obj.D] = svd(obj.mat,'econ');
            if nargin < 2
                r = sum(diag(obj.S) > 1e-14);
            end
            obj.C = obj.C(:,1:r);
            obj.S = obj.S(1:r,1:r);
            obj.D = obj.D(:,1:r);
        end
        
        
        function [] = mat2Vec(obj)
            assert(~isempty(obj.mat));    
            if ~isempty(obj.iper)
                tmp = reshape(obj.mat,[],1);
                obj.vec = tmp(obj.iper);
            else
                obj.vec = convertMattoVec(obj.x.nodes,obj.v.nodes,obj.k,obj.mat);
            end
        end
        
        
        function [] = vec2Mat(obj)
            assert(~isempty(obj.vec));
            if ~isempty(obj.per)
                tmp = obj.vec(obj.per);
                obj.mat = reshape(tmp,obj.x.N*(obj.k+1),obj.v.N*(obj.k+1));
            else
                obj.mat = convertVectoMat(obj.x.nodes,obj.v.nodes,obj.k,obj.vec);
            end
        end
        
 
        function r = rank(obj)
            assert(~isempty(obj.S));
            r = size(obj.S,1);
        end
        
        
        function plot(obj)
            assert(~isempty(obj.vec),'Run obj.getVec() to create vector before calling obj.plot()');
            plotVec(obj.x.nodes,obj.v.nodes,obj.k,obj.vec);
        end
    end
end



